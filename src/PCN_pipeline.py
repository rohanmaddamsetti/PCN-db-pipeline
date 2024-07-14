#!/usr/bin/env python

"""
PCN_pipeline.py by Maggie Wilson and Rohan Maddamsetti.

For this pipeline to work, ncbi datasets, pysradb, and kallisto must be in the $PATH.
On the Duke Compute Cluster (DCC), run the following to get these programs into the path:
conda activate PCNdb-env

IMPORTANT: this script assumes it is being run on DCC if sys.platform == "linux".
This means that users on a linux machine will need to modify a couple functions if they
are running this code locally, and cannot use slurm to submit many jobs in parallel.

Currently, this pipeline only analyzes Illumina short-read data.
Empirically, long-read data do not give accurate PCN information.

"""

import subprocess
import json
import os
import sys
import gzip
import re
import csv
import math
from Bio import SeqIO
from os.path import basename, exists
import urllib.request
import time
import logging
import glob
import pprint
import polars as pl
from tqdm import tqdm
import HTSeq ## for filtering fastq multireads.
import numpy as np ## for matrix multiplications for running PIRA.


################################################################################
## Functions.

def get_SRA_ID_from_RefSeqID(refseq_id):
    ## datasets must be in $PATH.
    bash_command = f'datasets summary genome accession {refseq_id}'
    cmd_output = subprocess.check_output(bash_command, shell=True)
    json_output = cmd_output.decode('utf-8')
    json_data = json.loads(json_output)
    sra_id = "NA"
    reports = json_data.get("reports")
    if reports:
        sample_ids = reports[0].get("assembly_info", {}).get("biosample", {}).get("sample_ids")
        if sample_ids:
            for sample_id in sample_ids:
                if sample_id.get("db") == "SRA":
                    sra_id = sample_id.get("value")
                    break
    return(sra_id)


def fetch_Run_IDs_with_pysradb(sra_id):
    
    ## pysradb must be in $PATH.
    pysradb_command = f'pysradb metadata {sra_id}'
    pysradb_attempts = 5
    pysra_command_worked = False
    while pysradb_attempts:
        try:
            pysradb_output = subprocess.check_output(pysradb_command, shell=True)
            ## assuming the previous line worked--
            pysradb_attempts = 0
            pysra_command_worked = True
        except subprocess.CalledProcessError:
            pysradb_attempts -= 1
    ## initialize an empty list of run_ids for this sra_id.
    run_accessions = list()
    ## check to see if pysradb_output is meaningful.
    if pysra_command_worked:
        pysradb_output_str = pysradb_output.decode('utf-8')
        # Splits the metadata of the SRA ID into respective rows. 
        # And isolates the rows that use the Illumina instrument.
        rows = pysradb_output_str.strip().split('\n')
        run_accessions = list()
        for i, row in enumerate(rows):
            if i == 0: continue ## skip the header
            fields = row.split("\t")
            try:
                study_accession, study_title, experiment_accession, experiment_title, experiment_desc, organism_taxid, organism_name, library_name, library_strategy, library_source, library_selection, library_layout, sample_accession, sample_title, instrument, instrument_model, instrument_model_desc, total_spots, total_size, run_accession, run_total_spots, run_total_bases = fields
                int_total_size = int(total_size)
            except ValueError: ## if either the number of fields is wrong, or if the typecast fails (say if total_size == "<NA>"),  then skip this run_accession.
                continue
            ## if there is data associated with this accession (total_size > 0),
            ## this is Illumina WGS data, and the run_accession is valid, then add to the list of run_accessions.
            if int_total_size > 0 and library_strategy == "WGS" and instrument_model_desc == "ILLUMINA" and run_accession != "nan":
                run_accessions.append(run_accession)
    return(run_accessions)


def create_RefSeq_SRA_RunID_table(prokaryotes_with_plasmids_file, RunID_table_outfile):
    ## first, get all RefSeq IDs in the prokaryotes-with-plasmids.txt file.
    with open(prokaryotes_with_plasmids_file, "r") as prok_with_plasmids_file_obj:
        prok_with_plasmids_lines = prok_with_plasmids_file_obj.read().splitlines()
    ## skip the header.
    prok_with_plasmids_data = prok_with_plasmids_lines[1:]
    ## get the right column (5th from end) and turn GCA Genbank IDs into GCF RefSeq IDs.
    refseq_id_column = [line.split("\t")[-5].replace("GCA", "GCF") for line in prok_with_plasmids_data]
    ## filter for valid IDs (some rows have a '-' as a blank placeholder).
    refseq_ids = [x for x in refseq_id_column if x.startswith("GCF")]
    ## now make the RunID csv file.
    with open(RunID_table_outfile, "w") as RunID_table_outfile_obj:
        header = "RefSeq_ID,SRA_ID,Run_ID\n"
        RunID_table_outfile_obj.write(header) 
        for RefSeq_accession in refseq_ids:
            my_SRA_ID = get_SRA_ID_from_RefSeqID(RefSeq_accession)
            if my_SRA_ID == "NA": continue ## skip genomes without SRA data.
            Run_IDs = fetch_Run_IDs_with_pysradb(my_SRA_ID)
            for my_Run_ID in Run_IDs:
                if my_Run_ID == "0" or my_Run_ID == "nan": continue ## skip bad Run_IDs
                row = f"{RefSeq_accession},{my_SRA_ID},{my_Run_ID}\n"
                print(row) ## just to show that the program is running properly.
                RunID_table_outfile_obj.write(row)
    return


def create_refseq_accession_to_ftp_path_dict(prokaryotes_with_plasmids_file):
    refseq_accession_to_ftp_path_dict = dict()
    ## first, get all RefSeq IDs in the prokaryotes-with-plasmids.txt file.
    with open(prokaryotes_with_plasmids_file, "r") as prok_with_plasmids_file_obj:
        for i, line in enumerate(prok_with_plasmids_file_obj):
            if i == 0: continue ## skip the header.
            ## get the accession field (5th from end) and turn GCA Genbank IDs into GCF RefSeq IDs.
            refseq_id = line.split("\t")[-5].replace("GCA", "GCF")        
            ## get the ftp_url field (3rd from end) and make sure that we turn the GCA Genbank URL
            ## into the GCF RefSeq FTP URL.
            ftp_url = line.split("\t")[-3].replace("GCA", "GCF")
            ## check for for valid IDs and URLs (some rows have a '-' as a blank placeholder).
            if refseq_id.startswith("GCF") and refseq_id in ftp_url:
                refseq_accession_to_ftp_path_dict[refseq_id] = ftp_url
    return refseq_accession_to_ftp_path_dict


def reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
    with open(md5_file, "r") as checksum_fh:
        target_string = "_genomic.gbff.gz"
        for line in checksum_fh:
            if target_string in line:          
                my_target_checksum, my_target_filename = line.strip().split()
                break
    if sys.platform == "darwin":
        my_md5_cmd = "md5" ## run md5 on my mac,
    elif sys.platform == "linux":
        my_md5_cmd = "md5sum" ## but run md5sum on DCC (linux)
    else:
        raise AssertionError("UNKNOWN PLATFORM")

    ## run md5 on the local file and get the output.
    md5_call = subprocess.run([my_md5_cmd, gbff_gz_file], capture_output=True, text=True)

    if sys.platform == "darwin": ## then the checksum is the LAST 'word' in the output.
        my_md5_checksum = md5_call.stdout.split()[-1].strip()
    elif sys.platform == "linux": ## then the checksum is the FIRST 'word' in the output.
        my_md5_checksum = md5_call.stdout.split()[0].strip()
    else:
        raise AssertionError("UNKNOWN PLATFORM")

    ## verify that the checksums match.
    return my_md5_checksum == my_target_checksum


def fetch_reference_genomes(RunID_table_file, refseq_accession_to_ftp_path_dict, reference_genome_dir, log_file):
    ## we get RefSeq IDs from the RunID table because this file *only* contains those RefSeq IDs 
    ## for which we could download raw Illumina short reads from the NCBI Short Read Archive.

    with open(RunID_table_file, "r") as RunID_file_obj:
        RunID_table_lines = RunID_file_obj.read().splitlines()

    ## remove the header from the imported data.
    RunID_table_data = RunID_table_lines[1:]
    ## get the first column to get all refseq_ids of interest.
    ## set comprehension to remove duplicates (there can be multiple SRA datasets per reference genome).
    refseq_ids = {line.split(",")[0] for line in RunID_table_data}
    ## now look up the FTP URLs for each refseq id.
    ftp_paths = [refseq_accession_to_ftp_path_dict[x] for x in refseq_ids]

    with open(log_file, 'w') as log_fh: ## for tracking which genomes are downloaded and which failed.
        for ftp_path in tqdm(ftp_paths):
            ## note that the format of this accession is {refseqid}_{assemblyid}.
            my_full_accession = basename(ftp_path)
            my_base_filename = my_full_accession + "_genomic.gbff.gz"
            ## files on the NCBI FTP site to download
            gbff_ftp_path = os.path.join(ftp_path, my_base_filename)
            md5_ftp_path = os.path.join(ftp_path, "md5checksums.txt")
            ## local paths to download these files
            gbff_gz_file = os.path.join(reference_genome_dir, my_base_filename)
            md5_file = os.path.join(reference_genome_dir, my_full_accession + "_md5checksums.txt")

            if exists(gbff_gz_file) and exists(md5_file): ## then check whether the reference genome is OK.
                if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
                    print(f"{gbff_gz_file} SUCCEEDED.\n", file=log_fh) ## print to the log file,
                    print(f"{gbff_gz_file} SUCCEEDED.\n") ## and print to stdout as well.
                    continue
                else:
                    os.remove(gbff_gz_file)
                    os.remove(md5_file)

            gbff_fetch_attempts = 5
            gbff_fetched = False

            while not gbff_fetched and gbff_fetch_attempts:
                try:
                    urllib.request.urlretrieve(gbff_ftp_path, filename=gbff_gz_file)
                    urllib.request.urlretrieve(md5_ftp_path, filename=md5_file)
                except urllib.error.URLError:
                    ## if some problem happens, try again.
                    gbff_fetch_attempts -= 1
                    if gbff_fetch_attempts == 0:
                        print(f"{gbff_gz_file} FAILED.", file=log_fh) ## print to the log file,
                        print(f"{gbff_gz_file} FAILED.") ## and print to stdout as well.
                    ## delete the corrupted files if they exist.
                    if exists(gbff_gz_file):
                        os.remove(gbff_gz_file)
                    if exists(md5_file):
                        os.remove(md5_file)
                ## if we are here, then assume the try block worked.
                if exists(gbff_gz_file) and exists(md5_file): ## then check whether the reference genome is OK.
                    if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
                        print(f"{gbff_gz_file} SUCCEEDED.", file=log_fh) ## print to the log file,
                        print(f"{gbff_gz_file} SUCCEEDED.") ## and print to stdout as well.
                        gbff_fetched = True  ## assume success if the checksum matches,
                        gbff_fetch_attempts = 0  ## and don't try again.
                    else:
                        os.remove(gbff_gz_file)
                        os.remove(md5_file)
    return
 

def get_Run_IDs_from_RunID_table(RunID_table_file):
    Run_IDs = list()
    with open(RunID_table_file, "r") as RunID_table_fh:
        table_csv = csv.DictReader(RunID_table_fh)
        Run_IDs = [row["Run_ID"] for row in table_csv]
    return Run_IDs


def download_fastq_reads(SRA_data_dir, Run_IDs):
        """
        the Run_ID has to be the last part of the directory.
        see documentation here:
        https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
        """
        for Run_ID in Run_IDs:
            prefetch_dir_path = os.path.join(SRA_data_dir, Run_ID)
            if os.path.exists(prefetch_dir_path): ## skip if we have already prefetched the read data.
                continue
            ## prefetch will create the prefetch_dir_path automatically-- give it the SRA_data_dir.
            ## default max-size is 20G, but some datasets are larger. So kick up the max-size.
            prefetch_args = ["prefetch", "--max-size", "100G", Run_ID, "-O", SRA_data_dir]
            print (" ".join(prefetch_args))
            subprocess.run(prefetch_args)
        print("prefetch step completed.")
        my_cwd = os.getcwd()
        os.chdir(SRA_data_dir)
        for Run_ID in Run_IDs:
            sra_fastq_file_1 = Run_ID + "_1.fastq"
            sra_fastq_file_2 = Run_ID + "_2.fastq"
            ## since we ran os.chdir(SRA_data_dir), this next line should work right.
            if os.path.exists(sra_fastq_file_1) and os.path.exists(sra_fastq_file_2):
                continue
            else:
                print ("Generating fastq for: " + Run_ID)
                fasterq_dump_args = ["fasterq-dump", "--threads", "10", Run_ID]
                print(" ".join(fasterq_dump_args))
                subprocess.run(fasterq_dump_args)
        ## now change back to original working directory.
        os.chdir(my_cwd)
        return


def all_fastq_data_exist(Run_IDs, SRA_data_dir):
    ## check to see if all the expected files exist on disk (does not check for corrupted data).
    SRA_file_list = os.listdir(SRA_data_dir)
    prefetch_dirs = [f for f in SRA_file_list if os.path.isdir(os.path.join(SRA_data_dir, f))]
    fastq_files = [f for f in SRA_file_list if f.endswith(".fastq")]
    for Run_ID in Run_IDs:
        ## does the prefetch directory exist?
        if Run_ID not in SRA_file_list:
            return False
        ## does the first fastq file exist?
        sra_fastq_path_1 = os.path.join(SRA_data_dir, Run_ID + "_1.fastq")
        if not os.path.exists(sra_fastq_path_1):
            return False
        ## does the second fastq file exist?
        sra_fastq_path_2 = os.path.join(SRA_data_dir, Run_ID + "_2.fastq")
        if  not os.path.exists(sra_fastq_path_2):
            return False
    return True
        

def generate_gene_level_fasta_reference_for_kallisto(gbk_gz_path, outfile):
    print("making as output: ", outfile)
    print("reading in as input:", gbk_gz_path)
    with open(outfile, "w") as outfh:
        with gzip.open(gbk_gz_path, 'rt') as gbk_gz_fh:
            SeqID = None
            SeqType = None
            for i, record in enumerate(SeqIO.parse(gbk_gz_fh, "genbank")):
                SeqID = record.id
                if "chromosome" in record.description or i == 0:
                    ## IMPORTANT: we assume here that the first record is a chromosome.
                    SeqType = "chromosome"
                elif "plasmid" in record.description:
                    SeqType = "plasmid"
                else:
                    continue
                for feature in record.features:
                    ## only analyze protein-coding genes.
                    if feature.type != "CDS": continue
                    locus_tag = feature.qualifiers["locus_tag"][0]
                    ## Important: for kallisto, we need to replace spaces with underscores in the product annotation field.
                    product = feature.qualifiers["product"][0].replace(" ","_")
                    DNAseq = feature.extract(record.seq)
                    header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"locus_tag="+locus_tag,"product="+product])
                    outfh.write(header + "\n")
                    outfh.write(str(DNAseq) + "\n")
    return


def make_NCBI_gene_fasta_refs_for_kallisto(refgenomes_dir, kallisto_ref_outdir):
    ## this function makes fasta sequences for every gene in every genome.
    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff.gz")]
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        genome_id = gzfile.split(".gbff.gz")[0]
        fasta_outfile = os.path.join(kallisto_ref_outdir, genome_id+".fna")
        generate_gene_level_fasta_reference_for_kallisto(gzpath, fasta_outfile)
    return


def generate_replicon_level_fasta_reference_for_kallisto(gbk_gz_path, outfile):
    print("making as output: ", outfile)
    print("reading in as input:", gbk_gz_path)
    with open(outfile, "w") as outfh:
        with gzip.open(gbk_gz_path, 'rt') as gbk_gz_fh:
            SeqID = None
            SeqType = None
            for i, record in enumerate(SeqIO.parse(gbk_gz_fh, "genbank")):
                SeqID = record.id
                if "chromosome" in record.description or i == 0:
                    ## IMPORTANT: we assume here that the first record is a chromosome.
                    SeqType = "chromosome"
                elif "plasmid" in record.description:
                    SeqType = "plasmid"
                else:
                    continue
                ## Important: for kallisto, we need to replace spaces with underscores in the replicon annotation field.
                replicon_description = record.description.replace(" ","_")
                header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"replicon="+replicon_description])
                outfh.write(header + "\n")
                outfh.write(str(record.seq) + "\n")
    return


def make_NCBI_replicon_fasta_refs_for_kallisto(refgenomes_dir, kallisto_ref_outdir):
    ## this function makes a genome fasta file for each genome.
    ## each genome fasta file contains fasta sequences for every replicon.
    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff.gz")]
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        genome_id = gzfile.split(".gbff.gz")[0]
        fasta_outfile = os.path.join(kallisto_ref_outdir, genome_id+".fna")
        generate_replicon_level_fasta_reference_for_kallisto(gzpath, fasta_outfile)
    return


def make_NCBI_kallisto_indices(kallisto_ref_dir, kallisto_index_dir):
    ref_fasta_filelist = [x for x in os.listdir(kallisto_ref_dir) if x.endswith(".fna")]
    for ref_fasta_file in ref_fasta_filelist:
        ref_fasta_path = os.path.join(kallisto_ref_dir, ref_fasta_file)
        genome_id = ref_fasta_file.split(".fna")[0]
        index_file = genome_id + ".idx"
        index_path = os.path.join(kallisto_index_dir, index_file)
        kallisto_index_args = ["kallisto", "index", "-i", index_path, ref_fasta_path]
        subprocess.run(kallisto_index_args)
    return


def run_kallisto_quant(RefSeq_to_SRA_RunList_dict, kallisto_index_dir, SRA_data_dir, results_dir):
    ## IMPORTANT: kallisto needs -l -s parameters supplied when run on single-end data.
    ## to avoid this complexity, I only process paired-end Illumina data, and skip single-end data altogether.
    index_list = [x for x in os.listdir(kallisto_index_dir) if x.endswith(".idx")]
    for index_file in index_list:
        index_path = os.path.join(kallisto_index_dir, index_file)
        genome_id = index_file.split(".idx")[0]
        refseq_id = "_".join(genome_id.split("_")[:2])
        Run_ID_list = RefSeq_to_SRA_RunList_dict[refseq_id]
        ## make read_path_arg_list.
        read_path_arg_list = list()
        for Run_ID in Run_ID_list:
            SRA_file_pattern = f"{SRA_data_dir}/{Run_ID}*.fastq"
            matched_fastq_list = sorted(glob.glob(SRA_file_pattern))
            if len(matched_fastq_list) != 2: ## skip if we didn't find paired fastq reads.
                continue
            read_path_arg_list += matched_fastq_list
        ## run kallisto quant with 10 threads by default.
        output_path = os.path.join(results_dir, genome_id)
        if len(read_path_arg_list): ## if we found paired-end fastq reads for this genome, then run kallisto.
            kallisto_quant_args = ["kallisto", "quant", "-t", "10", "-i", index_path, "-o", output_path, "-b", "100"] + read_path_arg_list
            kallisto_quant_string = " ".join(kallisto_quant_args)
            slurm_string = "sbatch -p scavenger --mem=16G --wrap=" + "\"" + kallisto_quant_string + "\""
            print(slurm_string)
            subprocess.run(slurm_string, shell=True)
    return


def make_RefSeq_to_SRA_RunList_dict(RunID_table_csv):
    RefSeq_to_SRA_RunList_dict = dict()
    with open(RunID_table_csv, "r") as csv_fh:
        for i, line in enumerate(csv_fh):
            if i == 0: continue ## skip the header.
            line = line.strip() 
            RefSeqID, SRA_ID, RunID = line.split(',')
            if RefSeqID in RefSeq_to_SRA_RunList_dict:
                RefSeq_to_SRA_RunList_dict[RefSeqID].append(RunID)
            else:
                RefSeq_to_SRA_RunList_dict[RefSeqID] = [RunID]
    return RefSeq_to_SRA_RunList_dict


def parse_gene_metadata_in_header(target_id):
    fields = target_id.split("|")
    SeqID = fields[0].split("=")[-1]
    SeqType = fields[1].split("=")[-1]
    locus_tag = fields[2].split("=")[-1]
    ## convert underscores back into spaces.
    product = fields[3].split("=")[-1].replace("_", " ")
    metadata_tuple = (SeqID, SeqType, locus_tag, product)
    return(metadata_tuple)


def estimate_chr_plasmid_copy_numbers_from_genes(genecount_tsv_path):
    genome_dict = dict()
    ## keys are SeqIDs.
    ## values are a dict: {SeqType: "chromosome", total_length: 10000, total_est_counts: 100}
    with open(genecount_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, locus_tag, product = parse_gene_metadata_in_header(target_id)
            if SeqID in genome_dict:
                genome_dict[SeqID]["total_length"] += float(length)
                genome_dict[SeqID]["total_est_counts"] += float(est_counts)
            else: ## Initialize the dictionary.
                genome_dict[SeqID] = {"SeqType" : SeqType, "total_length" : float(length), "total_est_counts": float(est_counts)}
    coverage_dict = dict()
    ##keys are seq_ids, value is (SeqType, coverage) pair.
    ## we set the default value to -1 so that we can catch error cases
    ## where the chromosome is not found in the genome.
    chromosome_coverage = -1
    for SeqID, replicon_dict in genome_dict.items():
        coverage = replicon_dict["total_est_counts"]/replicon_dict["total_length"]
        coverage_dict[SeqID] = (replicon_dict["SeqType"], coverage)
        if replicon_dict["SeqType"] == "chromosome":
            chromosome_coverage = coverage
            
    ## now normalize by chromosome coverage to get copy number estimates.
    copy_number_dict = dict()
    for SeqID, value_tuple in coverage_dict.items():
        seqtype, coverage = value_tuple
        copy_number_dict[SeqID] = (seqtype, coverage/chromosome_coverage)
    return(copy_number_dict)


def calculate_NCBI_replicon_copy_numbers_from_genes(kallisto_quant_results_dir, copy_number_csv_file):
    """
    define lists to encode the following columns of the table.
    AnnotationAccession, SeqID, SeqType, CopyNumber
    """
    AnnotationAccessionVec = []
    SeqIDVec = []
    SeqTypeVec = []
    CopyNumberVec = []
    ## skip .DS_Store and any other weird files.
    genomedirectories = [x for x in os.listdir(kallisto_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        ## I probably should have trimmed the '_genomic' suffix in an earlier step.
        annotation_accession = genomedir.split("_genomic")[0]
        genome_quantfile_path = os.path.join(kallisto_quant_results_dir, genomedir, "abundance.tsv")
        copy_number_dict = estimate_chr_plasmid_copy_numbers_from_genes(genome_quantfile_path)
        for SeqID, value_tuple in copy_number_dict.items():
            seqtype, coverage = value_tuple
            AnnotationAccessionVec.append(annotation_accession)
            SeqIDVec.append(SeqID)
            SeqTypeVec.append(seqtype)
            CopyNumberVec.append(coverage)

    assert len(AnnotationAccessionVec) == len(SeqIDVec) == len(SeqTypeVec) == len(CopyNumberVec)
    ## now write the copy number data to file.
    with open(copy_number_csv_file, "w") as outfh:
        header = "AnnotationAccession,SeqID,SeqType,CopyNumber"
        outfh.write(header + "\n")
        for i in range(len(AnnotationAccessionVec)):
            outfh.write(AnnotationAccessionVec[i] + "," + SeqIDVec[i] + "," + SeqTypeVec[i] + "," + str(CopyNumberVec[i]) + "\n")
    return


def estimate_gene_copy_numbers(genecount_tsv_path):

    chromosomal_gene_length = 0.0
    chromosomal_gene_est_counts = 0.0

    gene_coverage_dict = dict()
    ## get the chromosomal gene coverage, and get the coverage for all genes
    with open(genecount_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, locus_tag, product = parse_gene_metadata_in_header(target_id)
            coverage = float(est_counts) / float(length)
            gene_coverage_dict[locus_tag] = (SeqID, SeqType, product, coverage)
            if SeqType == "chromosome":
                chromosomal_gene_length += float(length)
                chromosomal_gene_est_counts += float(est_counts)
    ## NOTE: GCF_026154285.1_ASM2615428v1 did not have any reads pseudoalign.
    ## Return an empty dict() when nothing aligns to the chromosome.
    if chromosomal_gene_length == 0:
        print("WARNING: no reads pseudoaligned to chromosome in file: ", genecount_tsv_path)
        print("estimate_gene_copy_numbers is returning an empty dict.")
        return(dict())
    chromosome_coverage = chromosomal_gene_est_counts / chromosomal_gene_length
    ## now normalize by chromosome coverage to get copy number estimates.
    gene_copy_number_dict = dict()
    for locus_tag, value_tuple in gene_coverage_dict.items():
        my_SeqID, my_SeqType, my_product, my_coverage = value_tuple
        my_gene_copy_number = my_coverage / chromosome_coverage
        gene_copy_number_dict[locus_tag] = (my_SeqID, my_SeqType, my_product, my_gene_copy_number)
    return(gene_copy_number_dict)


def measure_NCBI_gene_copy_numbers(kallisto_gene_quant_results_dir, gene_copy_number_csv_file):
    """
    define lists to encode the following columns of the table.
    RefSeqID, SeqID, SeqType, locus_tag, product, CopyNumber
    """
    RefSeqIDVec = []
    SeqIDVec = [] ## this is for the replicon.
    SeqTypeVec = []
    LocusTagVec = []
    ProductVec = []
    CopyNumberVec = []
    
    ## skip .DS_Store and any other weird files.
    genomedirectories = [x for x in os.listdir(kallisto_gene_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        refseq_id = "_".join(genomedir.split("_")[:2])
        genome_quantfile_path = os.path.join(kallisto_gene_quant_results_dir, genomedir, "abundance.tsv")
        gene_copy_number_dict = estimate_gene_copy_numbers(genome_quantfile_path)
        for locus_tag, value_tuple in gene_copy_number_dict.items():
            SeqID, seqtype, product, copy_number = value_tuple
            RefSeqIDVec.append(refseq_id)
            SeqIDVec.append(SeqID)
            SeqTypeVec.append(seqtype)
            LocusTagVec.append(locus_tag)
            ProductVec.append(product)
            CopyNumberVec.append(str(copy_number))

    assert len(RefSeqIDVec) == len(SeqIDVec) == len(SeqTypeVec) == len(LocusTagVec) == len(ProductVec) == len(CopyNumberVec)

    ## we have to double-quote all columns-- some fields in the product column contain commas!
    RefSeqIDVec = ["\"" + x + "\"" for x in RefSeqIDVec]
    SeqIDVec = ["\"" + x + "\"" for x in SeqIDVec]
    SeqTypeVec = ["\"" + x + "\"" for x in SeqTypeVec]
    LocusTagVec = ["\"" + x + "\"" for x in LocusTagVec]
    ProductVec = ["\"" + x + "\"" for x in ProductVec]
    CopyNumberVec = ["\"" + x + "\"" for x in CopyNumberVec]
    
    ## now write the gene copy number data to file.
    with open(gene_copy_number_csv_file, "w") as outfh:
        ## double-quote each column name in the header for consistency.
        header = "\"RefSeqID\",\"SeqID\",\"SeqType\",\"locus_tag\",\"product\",\"CopyNumber\""
        outfh.write(header + "\n")
        for i in range(len(RefSeqIDVec)):
            outfh.write(RefSeqIDVec[i] + "," + SeqIDVec[i] + "," + SeqTypeVec[i] + "," + LocusTagVec[i] + "," + ProductVec[i] + "," + CopyNumberVec[i] + "\n")
    return


def parse_replicon_metadata_in_header(target_id):
    fields = target_id.split("|")
    SeqID = fields[0].split("=")[-1]
    SeqType = fields[1].split("=")[-1]
    ## convert underscores back into spaces.
    replicon_description = fields[2].split("=")[-1].replace("_", " ")
    metadata_tuple = (SeqID, SeqType, replicon_description)
    return(metadata_tuple)


def estimate_replicon_copy_numbers(kallisto_replicon_count_tsv_path):

    chromosomal_length = 0.0
    chromosomal_est_counts = 0.0

    replicon_coverage_dict = dict()
    ## get the chromosomal gene coverage, and get the coverage for all genes
    with open(kallisto_replicon_count_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, replicon_description = parse_replicon_metadata_in_header(target_id)
            coverage = float(est_counts) / float(length)
            replicon_coverage_dict[SeqID] = (SeqType, replicon_description, coverage)
            if SeqType == "chromosome":
                chromosomal_length += float(length)
                chromosomal_est_counts += float(est_counts)
    ## Return an empty dict() when nothing aligns to the chromosome.
    if chromosomal_length == 0:
        print("WARNING: no reads pseudoaligned to chromosome in file: ", kallisto_replicon_count_tsv_path)
        print("estimate_replicon_copy_numbers is returning an empty dict.")
        return(dict())
    chromosome_coverage = chromosomal_est_counts / chromosomal_length
    ## now normalize by chromosome coverage to get copy number estimates.
    replicon_copy_number_dict = dict()
    for my_SeqID, value_tuple in replicon_coverage_dict.items():
        my_SeqType, replicon_description, my_coverage = value_tuple
        my_replicon_copy_number = my_coverage / chromosome_coverage
        replicon_copy_number_dict[my_SeqID] = (my_SeqType, replicon_description, my_replicon_copy_number)
    return(replicon_copy_number_dict)


def measure_NCBI_replicon_copy_numbers(kallisto_replicon_quant_results_dir, replicon_copy_number_csv_file):
    """
    define lists to encode the following columns of the table.
    RefSeqID, SeqID, SeqType, CopyNumber
    """
    RefSeqIDVec = []
    SeqIDVec = [] ## this is for the replicon.
    SeqTypeVec = []
    RepliconDescriptionVec = []
    CopyNumberVec = []
    
    ## skip .DS_Store and any other weird files.
    genomedirectories = [x for x in os.listdir(kallisto_replicon_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        refseq_id = "_".join(genomedir.split("_")[:2])
        genome_quantfile_path = os.path.join(kallisto_replicon_quant_results_dir, genomedir, "abundance.tsv")
        replicon_copy_number_dict = estimate_replicon_copy_numbers(genome_quantfile_path)
        for SeqID, value_tuple in replicon_copy_number_dict.items():
            seqtype, replicon_description, copy_number = value_tuple
            RefSeqIDVec.append(refseq_id)
            SeqIDVec.append(SeqID)
            SeqTypeVec.append(seqtype)
            RepliconDescriptionVec.append(replicon_description)
            CopyNumberVec.append(str(copy_number))

    assert len(RefSeqIDVec) == len(SeqIDVec) == len(SeqTypeVec) == len(RepliconDescriptionVec) == len(CopyNumberVec)

    ## we have to double-quote all columns-- some fields in the product column contain commas!
    RefSeqIDVec = ["\"" + x + "\"" for x in RefSeqIDVec]
    SeqIDVec = ["\"" + x + "\"" for x in SeqIDVec]
    SeqTypeVec = ["\"" + x + "\"" for x in SeqTypeVec]
    RepliconDescriptionVec = ["\"" + x + "\"" for x in RepliconDescriptionVec]
    CopyNumberVec = ["\"" + x + "\"" for x in CopyNumberVec]
    
    ## now write the replicon copy number data to file.
    with open(replicon_copy_number_csv_file, "w") as outfh:
        ## double-quote each column name in the header for consistency.
        header = "\"RefSeqID\",\"SeqID\",\"SeqType\",\"RepliconDescription\",\"CopyNumber\""
        outfh.write(header + "\n")
        for i in range(len(RefSeqIDVec)):
            outfh.write(RefSeqIDVec[i] + "," + SeqIDVec[i] + "," + SeqTypeVec[i] + "," + RepliconDescriptionVec[i] + "," + CopyNumberVec[i] + "\n")
    return


def tabulate_NCBI_replicon_lengths(refgenomes_dir, replicon_length_csv_file):
    with open(replicon_length_csv_file, 'w') as outfh:
        header = "AnnotationAccession,SeqID,SeqType,replicon_length\n"
        outfh.write(header)
        for gbk_gz in os.listdir(refgenomes_dir):
            if not gbk_gz.endswith(".gbff.gz"): continue
            annotation_accession = gbk_gz.split("_genomic.gbff")[0]
            infile = os.path.join(refgenomes_dir, gbk_gz)
            with gzip.open(infile, "rt") as genome_fh:
                for i, replicon in enumerate(SeqIO.parse(genome_fh, "gb")):
                    SeqID = replicon.id
                    if "chromosome" in replicon.description or i == 0:
                        ## IMPORTANT: we assume here that the first record is a chromosome.
                        SeqType = "chromosome"
                    elif "plasmid" in replicon.description:
                        SeqType = "plasmid"
                    else:
                        continue
                    replicon_length = str(len(replicon))
                    ## now write out the data for the replicon.
                    row = ','.join([annotation_accession, SeqID, SeqType, replicon_length])
                    outfh.write(row + "\n")
    return


def filter_gene_copy_number_file_for_ARGs(gene_copy_number_csv_file, ARG_copy_number_csv_file):
    ## define ARG keywords for pattern matching.
    chloramphenicol_keywords = "chloramphenicol|Chloramphenicol"
    tetracycline_keywords = "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
    MLS_keywords = "macrolide|lincosamide|streptogramin"
    multidrug_keywords = "Multidrug resistance|multidrug resistance|antibiotic resistance"
    beta_lactam_keywords = "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\\S*"
    glycopeptide_keywords = "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
    polypeptide_keywords = "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
    diaminopyrimidine_keywords = "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
    sulfonamide_keywords = "sulfonamide|Sul1|sul1|sulphonamide"
    quinolone_keywords = "quinolone|Quinolone|oxacin|qnr|Qnr"
    aminoglycoside_keywords = "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
    macrolide_keywords = "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
    antimicrobial_keywords = "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\\S*"

    antibiotic_keywords = "|".join([chloramphenicol_keywords, tetracycline_keywords, MLS_keywords, multidrug_keywords,
                                    beta_lactam_keywords, glycopeptide_keywords, polypeptide_keywords, diaminopyrimidine_keywords,
                                    sulfonamide_keywords, quinolone_keywords, aminoglycoside_keywords, macrolide_keywords, antimicrobial_keywords])

    with open(gene_copy_number_csv_file, "r") as gene_fh, open(ARG_copy_number_csv_file, "w") as ARG_fh:
        for i, line in enumerate(gene_fh):
            if i == 0: ## dealing with the header
                ARG_fh.write(line)
            else:
                if re.search(antibiotic_keywords, line):
                    ARG_fh.write(line)
    return


def generate_replicon_fasta_references_for_themisto(gbk_gz_path, fasta_outdir):
    print("reading in as input:", gbk_gz_path)
    ## open the input reference genome file.
    with gzip.open(gbk_gz_path, 'rt') as gbk_gz_fh:
        SeqID = None
        SeqType = None
        for i, record in enumerate(SeqIO.parse(gbk_gz_fh, "genbank")):
            SeqID = record.id
            if "chromosome" in record.description or i == 0:
                ## IMPORTANT: we assume here that the first record is a chromosome.
                SeqType = "chromosome"
            elif "plasmid" in record.description:
                SeqType = "plasmid"
            else:
                continue
            ## replace spaces with underscores in the replicon annotation field.
            replicon_description = record.description.replace(" ","_")
            header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"replicon="+replicon_description])
            my_replicon_fastafile = SeqID + ".fna"
            my_replicon_outfilepath = os.path.join(fasta_outdir, my_replicon_fastafile)
            with open(my_replicon_outfilepath, "w") as outfh:
                outfh.write(header + "\n")
                outfh.write(str(record.seq) + "\n")
    return


def generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir):
    """
    IMPORTANT: the replicons in the fasta reference list file need to be sorted by length.
    Directly measure the length of each of the replicons here, and then use this to sort
    the list of replicon FASTA paths, before writing to file.
    """
    genome_id = os.path.basename(fasta_outdir)
    ## IMPORTANT: exclude any FASTA file that has the genome_id as a name.
    ## (this is a file containing all the replicons, created and used in the PIRA stages downstream).
    replicon_fasta_filelist = [x for x in os.listdir(fasta_outdir) if x.endswith(".fna") and not x.startswith(genome_id)]
    replicon_fasta_pathlist = [os.path.join(fasta_outdir, x) for x in replicon_fasta_filelist]

    ## Decorate-sort-undecorate by fasta sequence length.
    decorated_replicon_fasta_pathlist = list()
    for fastapath in replicon_fasta_pathlist:
        my_replicon = SeqIO.read(fastapath, "fasta")
        ## Get the length of the replicon
        replicon_length = len(my_replicon.seq)
        ## append a tuple of (fastapath, replicon_length)
        my_tuple = (fastapath, replicon_length)
        decorated_replicon_fasta_pathlist.append(my_tuple)
    ## sort the decorated list in place by replicon_length, in descending order from
    ## from largest to smallest replicon.
    decorated_replicon_fasta_pathlist.sort(key=lambda x: x[1], reverse=True)
    ## undecorate the path list, which is now in sorted order.
    sorted_replicon_fasta_pathlist = [x[0] for x in decorated_replicon_fasta_pathlist]

    ## write the path list to file for themisto.
    replicon_listfile = os.path.join(fasta_outdir, genome_id + ".txt")
    with open(replicon_listfile, "w") as fastatxtfile_fh:
        for fasta_path in sorted_replicon_fasta_pathlist:
            fastatxtfile_fh.write(fasta_path + "\n")
    return


def make_NCBI_replicon_fasta_refs_for_themisto(refgenomes_dir, themisto_fasta_ref_outdir):
    ## this function makes a genome directory for each genome.
    ## each directory contains separate fasta files for each replicon.

    ## make the output directory if it does not exist.
    if not exists(themisto_fasta_ref_outdir):
        os.mkdir(themisto_fasta_ref_outdir)

    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff.gz")]
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        genome_id = gzfile.split(".gbff.gz")[0]
        fasta_outdir = os.path.join(themisto_fasta_ref_outdir, genome_id)
        ## make the fasta output directory if it does not exist.
        if not exists(fasta_outdir):
            os.mkdir(fasta_outdir)
        generate_replicon_fasta_references_for_themisto(gzpath, fasta_outdir)
        generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir)
    return


def make_NCBI_themisto_indices(themisto_ref_dir, themisto_index_dir):
    ## make the output directory if it does not exist.
    if not exists(themisto_index_dir):
        os.mkdir(themisto_index_dir)

    ## each directory is named after the genome_id of the given genome.
    for genome_id in os.listdir(themisto_ref_dir):
        ## get the actual path for this directory
        ref_fasta_dir = os.path.join(themisto_ref_dir, genome_id)
        ## make sure that this path is real, and not an artifact of some weird file like .DS_Store in this directory
        if not os.path.isdir(ref_fasta_dir):
            continue
        index_input_filelist = os.path.join(ref_fasta_dir, genome_id + ".txt")

        ## make the directory for the index files, if it does not exist.
        genome_index_dir = os.path.join(themisto_index_dir, genome_id)
        if not exists(genome_index_dir):
            os.mkdir(genome_index_dir)
        ## set the index_prefix to write index files into the genome_index_dir.
        index_prefix = os.path.join(genome_index_dir, genome_id)
        
        ## make the temp directory if it doesn't exist.
        tempdir = os.path.join(genome_index_dir, "temp")
        if not exists(tempdir):
            os.mkdir(tempdir)
        
        themisto_build_args = ["themisto", "build", "-k","31", "-i", index_input_filelist, "--index-prefix", index_prefix, "--temp-dir", tempdir, "--mem-gigas", "4", "--n-threads", "4", "--file-colors"]
        themisto_build_string = " ".join(themisto_build_args)
        slurm_string = "sbatch -p scavenger --mem=4G --cpus-per-task=4 --wrap=" + "\"" + themisto_build_string + "\""
        print(slurm_string)
        subprocess.run(slurm_string, shell=True)
    return


def run_themisto_pseudoalign(RefSeq_to_SRA_RunList_dict, themisto_index_dir, SRA_data_dir, themisto_pseudoalignment_dir):
    ## make the output directory if it does not exist.
    if not exists(themisto_pseudoalignment_dir):
        os.mkdir(themisto_pseudoalignment_dir)

    ## each directory in themisto_index_dir is named after the genome_id of the given genome.
    genome_id_list = [x for x in os.listdir(themisto_index_dir)]
    for genome_id in genome_id_list:
        my_index_dir = os.path.join(themisto_index_dir, genome_id)
        ## make sure that this path is real, and not an artifact of some weird file in this directory
        if not os.path.isdir(my_index_dir):
            continue
        ## make the arguments relating to the paths to the index files for this genome.
        my_index_prefix = os.path.join(my_index_dir, genome_id)

        ## make the file containing the paths to the sequencing read data for this genome.
        refseq_id = "_".join(genome_id.split("_")[:2])

        ## check for inconsistencies between Run_ID_list and the genomes with themisto indices
        if refseq_id not in RefSeq_to_SRA_RunList_dict.keys():
            print(f"INCONSISTENCY ERROR: {refseq_id} not found in RunID_table.csv!")
            continue
        
        Run_ID_list = RefSeq_to_SRA_RunList_dict[refseq_id]
        ## make read_path_arg_list.
        readpath_list = list()
        for Run_ID in Run_ID_list:
            SRA_file_pattern = f"{SRA_data_dir}/{Run_ID}*.fastq"
            matched_fastq_list = sorted(glob.glob(SRA_file_pattern))
            readpath_list += matched_fastq_list

        ## if we didn't find any fastq data for this genome, then skip to the next genome.
        if not len(readpath_list):
            continue
        
        ## now, write the paths for the SRA read data to disk for themisto pseudoalign.
        SRAdata_listfile = os.path.join(themisto_pseudoalignment_dir, genome_id + "_SRAdata.txt")
        with open(SRAdata_listfile, "w") as SRAtxtfile_fh:
            for readpath in readpath_list:
                SRAtxtfile_fh.write(readpath + "\n")

        ## make the directory for the pseudoalignments.
        my_pseudoalignment_output_dir = os.path.join(themisto_pseudoalignment_dir, genome_id)
        if not exists(my_pseudoalignment_output_dir):
            os.mkdir(my_pseudoalignment_output_dir)

        ## make the temp directory if it doesn't exist.
        tempdir = os.path.join(my_pseudoalignment_output_dir, "temp")
        if not exists(tempdir):
            os.mkdir(tempdir)
            
        ## we make corresponding pseudoalignment output files for each SRA dataset.
        ## This list goes into the output listfile.
        output_listfile = os.path.join(themisto_pseudoalignment_dir, genome_id + "_pseudoalignments.txt")
        with open(output_listfile, "w") as output_listfile_fh:
            for readpath in readpath_list:
                read_filename = os.path.basename(readpath).split(".fastq")[0]
                output_filename = os.path.join(my_pseudoalignment_output_dir, read_filename + "_pseudoalignment.txt")
                output_listfile_fh.write(output_filename + "\n")
            
        ## now run themisto pseudoalign.
        themisto_pseudoalign_args = ["themisto", "pseudoalign", "--query-file-list", SRAdata_listfile, "--index-prefix", my_index_prefix, "--temp-dir", tempdir, "--out-file-list", output_listfile, "--n-threads", "4", "--threshold", "0.7"]
        themisto_pseudoalign_string = " ".join(themisto_pseudoalign_args)
        slurm_string = "sbatch -p scavenger --mem=4G --cpus-per-task=4 --wrap=" + "\"" + themisto_pseudoalign_string + "\""
        print(slurm_string)
        subprocess.run(slurm_string, shell=True)
    return


def map_themisto_IDs_to_replicon_metadata(themisto_replicon_ref_dir, my_genome_dirname):
    ## now, let's map the themisto replicon ID numbers to a (SeqID, SeqType) tuple.
    themisto_ID_to_seq_metadata_dict = dict()
    my_themisto_replicon_ID_mapping_file = os.path.join(themisto_replicon_ref_dir, my_genome_dirname, my_genome_dirname + ".txt")
    with open(my_themisto_replicon_ID_mapping_file, "r") as replicon_ID_mapping_fh:
        for i,fasta_file_path in enumerate(replicon_ID_mapping_fh):
            fasta_file_path = fasta_file_path.strip()
            ## get the header from this fasta file.
            with open(fasta_file_path, "r") as fasta_fh:
                my_header = fasta_fh.readline().strip()
                my_seq = fasta_fh.readline().strip()
            fields = my_header.strip(">").split("|")
            SeqID = fields[0].split("=")[-1]
            SeqType = fields[1].split("=")[-1]
            SeqLength = len(my_seq)
            themisto_ID_to_seq_metadata_dict[i] = (SeqID, SeqType, SeqLength)
    return themisto_ID_to_seq_metadata_dict


def summarize_themisto_pseudoalignment_results(themisto_replicon_ref_dir, themisto_pseudoalignment_dir, themisto_results_csvfile_path):
    with open(themisto_results_csvfile_path, "w") as output_csv_fh:
        ## first, write the output header.
        output_header = "AnnotationAccession,SeqID,SeqType,ReadCount"
        output_csv_fh.write(output_header + "\n")
        print(output_header)
        ## Iterate over the directories in themisto_pseudoalignment_dir.
        ## These contain the pseudoalignments for each genome.
        themisto_pseudoalignment_result_dirs = [x for x in os.listdir(themisto_pseudoalignment_dir) if x.startswith("GCF") and os.path.isdir(os.path.join(themisto_pseudoalignment_dir, x))]
        for my_genome_dirname in themisto_pseudoalignment_result_dirs:
            if not my_genome_dirname.startswith("GCF"): continue ## just an additional check to remove the temp directory.
            my_cur_pseudoalignment_dir_path = os.path.join(themisto_pseudoalignment_dir, my_genome_dirname)
            my_cur_AnnotationAccession = my_genome_dirname.strip("_genomic")
            ## initialize a dictionary to store the pseudoalignment counts.
            pseudoalignment_read_count_dict = dict()
            pseudoalignment_filepaths = [os.path.join(my_cur_pseudoalignment_dir_path, x) for x in os.listdir(my_cur_pseudoalignment_dir_path) if x.endswith("_pseudoalignment.txt")]
            ## now, summarize the read counts for this genome.
            for my_pseudoalignment_filepath in pseudoalignment_filepaths:
                with open(my_pseudoalignment_filepath) as my_pseudoalignment_fh:
                    for line in my_pseudoalignment_fh:
                        ## handle odd behavior in themisto: we need to sort the themisto replicon ID numbers ourselves.
                        ## IMPORTANT: sort the themisto ids (strings) by their numerical value.
                        replicon_set_string = " ".join(sorted(line.strip().split()[1:], key=lambda x: int(x)))
                        if replicon_set_string in pseudoalignment_read_count_dict:
                            pseudoalignment_read_count_dict[replicon_set_string] += 1
                        else:
                            pseudoalignment_read_count_dict[replicon_set_string] = 1

            ## now, let's map the themisto replicon ID numbers to a (SeqID, SeqType) tuple.
            themisto_ID_to_seq_metadata_dict = map_themisto_IDs_to_replicon_metadata(themisto_replicon_ref_dir, my_genome_dirname)
            
            ## now write the pseudoalignment counts to file.
            for replicon_set_string in sorted(pseudoalignment_read_count_dict.keys()):
                read_count = pseudoalignment_read_count_dict[replicon_set_string]
                replicon_ID_list = replicon_set_string.split()
                if len(replicon_ID_list) == 0:
                    SeqID = "NA"
                    SeqType = "NA"
                elif len(replicon_ID_list) == 1:
                    my_replicon_ID = replicon_ID_list[0]
                    ## the SeqLength parameter in themisto_ID_to_seq_metadata_dict is not used here.
                    SeqID, SeqType, _ = themisto_ID_to_seq_metadata_dict[my_replicon_ID]
                else:
                    SeqID = "&".join([themisto_ID_to_seq_metadata_dict[replicon_ID][0] for replicon_ID in replicon_ID_list])
                    SeqType = "multireplicon_sequence"
                ## now write to file.
                rowdata = ",".join([my_cur_AnnotationAccession, SeqID, SeqType, str(read_count)])
                print(rowdata)
                output_csv_fh.write(rowdata + "\n")
    return


def naive_themisto_PCN_estimation(themisto_results_csv_file, replicon_length_csv_file, naive_themisto_PCN_csv_file):
    ## This function simply ignores multireplicon reads when estimating PCN.
    ##Also note that this function omits replicons with zero mapped reads.
    print("running naive themisto PCN estimation (ignoring multireplicon reads)")
    ## import the data as polars dataframes.
    replicon_length_df = pl.read_csv(replicon_length_csv_file)
    naive_themisto_read_count_df = pl.read_csv(themisto_results_csv_file).filter(
        (pl.col("SeqType") == "chromosome") | (pl.col("SeqType") == "plasmid")).join(
            replicon_length_df, on = "SeqID").with_columns(
                (pl.col("ReadCount") / pl.col("replicon_length")).alias("SequencingCoverage"))
    
    ## make a second dataframe containing just the sequencing coverage for the longest replicon for each genome.
    ## to do so, first group by AnnotationAccession and compute maximum replicon_length within each group.
    longest_replicon_df = naive_themisto_read_count_df.group_by(
        "AnnotationAccession").agg(pl.col("replicon_length").max()).join(
            ## now join with the original DataFrame to filter for rows with the maximum replicon_length
            naive_themisto_read_count_df, on=["AnnotationAccession", "replicon_length"], how="inner").select(
                pl.col("AnnotationAccession", "SequencingCoverage")).with_columns(
                    pl.col("SequencingCoverage").alias("LongestRepliconCoverage")).select(
                        pl.col("AnnotationAccession", "LongestRepliconCoverage"))

    ## now normalize SequencingCoverage by LongestRepliconCoverage for each genome to calculate PCN.
    naive_themisto_PCN_df = naive_themisto_read_count_df.join(
        ## WEIRD BEHAVIOR in polars: the AnnotationAccession key for both dataframes is preserved,
        ## I don't know why. So remove the newly made AnnotationAccession_right column.
        longest_replicon_df, on = "AnnotationAccession").select(pl.col("*").exclude("AnnotationAccession_right")).with_columns(
            (pl.col("SequencingCoverage") / pl.col("LongestRepliconCoverage")).alias("CopyNumber"))

    ## now write the naive PCN estimates to file.
    naive_themisto_PCN_df.write_csv(naive_themisto_PCN_csv_file)    
    return


def assign_multireplicon_reads(genome_df):
    ## create a new data frame with just chromosome and plasmid data.
    updated_df = genome_df.filter(
        (pl.col("SeqType") == "chromosome") | (pl.col("SeqType") == "plasmid")).with_columns(
            ## cast the ReadCount to floats, and
            ## replace null values in the ReadCount columns with floating-point zeros.
            pl.col("ReadCount").cast(pl.Float64).fill_null(pl.lit(0.0)))

    ## get the multireplicon data.
    multireplicon_df = genome_df.filter(pl.col("SeqType") == "multireplicon_sequence")
    ## iterate over each multireplicon read set.
    for row_dict in multireplicon_df.iter_rows(named=True):
        seq_id_list = row_dict["SeqID"].split("&")
        read_count = row_dict["ReadCount"]
        num_reads_to_assign_to_each_replicon = float(read_count) / float(len(seq_id_list))
        ## iterate over each replicon in this multireplicon read set.
        for seq_id in seq_id_list:
            ## update the relevant values in the dataframe
            old_seqID_readcount = updated_df.filter(updated_df["SeqID"] == seq_id)["ReadCount"][0]
            new_seqID_readcount = old_seqID_readcount + num_reads_to_assign_to_each_replicon
            ## a kludgy hacky solution, but works:
            ## create a new column called temp, that has the new_seqID_readcount for the given SeqID,
            ## and zeros elsewhere in the column.
            ## then make a new column called updatedReadCount that takes the maximum of ReadCount and temp.
            ## then drop ReadCount and temp and rename updatedReadCount as ReadCount.
            updated_df = updated_df.with_columns(
                "ReadCount",
                pl.when(pl.col("SeqID") == seq_id).then(new_seqID_readcount).otherwise(pl.lit(0.0)).alias("temp")
            ).with_columns(pl.max_horizontal("ReadCount", "temp").alias("updatedReadCount")).select(
                pl.col("*").exclude(["ReadCount", "temp"])).rename({"updatedReadCount":"ReadCount"})
    return updated_df


def simple_themisto_PCN_estimation(themisto_results_csv_file, replicon_length_csv_file, simple_themisto_PCN_csv_file):
    ## This function divides multireplicon reads equally among the relevant replicons.

    ## Use a split-apply-combine strategy to break the big dataframe into smaller dataframes for each genome.
    ## Then pass each data frame to a function that assigns multi-replicon reads to each of the replicons,
    ## and produces a data frame that only has plasmids and chromosomes with the updated read counts.
    ## then merge all these pieces together into a big dataframe.
    ## finally, calculate PCN with this updated dataframe.
    print("running simple themisto PCN estimation (partitioning multireplicon reads equally among replicons)")
    ## import the data as polars dataframes.
    replicon_length_df = pl.read_csv(replicon_length_csv_file)
    themisto_read_count_df = pl.read_csv(themisto_results_csv_file)

    ## initialize the list of dataframes for each genome with updated read counts.
    list_of_updated_dfs = []
    ## iterate over the data for each genome
    for annotation_accession_group, genome_reads_df in themisto_read_count_df.group_by(["AnnotationAccession"]):
        ## awkward notation to avoid deprecated behavior for iterating over groups in polars.
        annotation_accession = annotation_accession_group[0]        
        ## IMPORTANT: some replicons may not have any reads assigned to them,
        ## say if they are completely contained in another replicon.
        ## To handle this case, we need to initialize genome_df to have all the replicons,
        ## with zeros assigned to ReadCount for replicons that are not present in genome_reads_df.
        ## We can get the missing replicons from replicon_length_df.
        genome_df = replicon_length_df.filter(pl.col("AnnotationAccession") == annotation_accession).join(
            genome_reads_df, on = ["AnnotationAccession","SeqID","SeqType"], how="outer_coalesce")
        
        ## now assign multireplicon reads equally among the replicons that they match.
        updated_df = assign_multireplicon_reads(genome_df)
        list_of_updated_dfs.append(updated_df)

    ## now merge the list of updated data frames.
    simple_themisto_read_count_df =  pl.concat(list_of_updated_dfs).join(
        ## and calculate SequencingCoverage
        replicon_length_df, on = "SeqID").with_columns(
            (pl.col("ReadCount") / pl.col("replicon_length")).alias("SequencingCoverage"))

    ## make a second dataframe containing just the sequencing coverage for the longest replicon for each genome.
    ## to do so, first group by AnnotationAccession and compute maximum replicon_length within each group.
    longest_replicon_df = simple_themisto_read_count_df.group_by(
        "AnnotationAccession").agg(pl.col("replicon_length").max()).join(
            ## now join with the original DataFrame to filter for rows with the maximum replicon_length
           simple_themisto_read_count_df, on=["AnnotationAccession", "replicon_length"], how="inner").select(
                pl.col("AnnotationAccession", "SequencingCoverage")).with_columns(
                    pl.col("SequencingCoverage").alias("LongestRepliconCoverage")).select(
                        pl.col("AnnotationAccession", "LongestRepliconCoverage"))

    ## now normalize SequencingCoverage by LongestRepliconCoverage for each genome to calculate PCN.
    simple_themisto_PCN_df = simple_themisto_read_count_df.join(
        ## WEIRD BEHAVIOR in polars: the AnnotationAccession key for both dataframes is preserved,
        ## I don't know why. So remove the newly made AnnotationAccession_right column.
        longest_replicon_df, on = "AnnotationAccession").select(pl.col("*").exclude("AnnotationAccession_right")).with_columns(
            (pl.col("SequencingCoverage") / pl.col("LongestRepliconCoverage")).alias("CopyNumber"))

    ## now write the simple PCN estimates to file.
    print(simple_themisto_PCN_df)
    simple_themisto_PCN_df.write_csv(simple_themisto_PCN_csv_file)
    return


def make_gbk_annotation_table(reference_genome_dir, gbk_annotation_file):
    with open(gbk_annotation_file,"w") as out_fh:
        header = "Annotation_Accession,host,isolation_source\n"
        out_fh.write(header)

        gbk_files = [x for x in os.listdir(reference_genome_dir) if x.endswith("_genomic.gbff.gz")]

        for x in tqdm(gbk_files):
            gbk_path = reference_genome_dir + x
            annotation_accession = x.split("_genomic.gbff.gz")[0]
            with gzip.open(gbk_path,'rt') as gbk_fh:
                host = "NA"
                isolation_source = "NA"
                ## use a buffer to store the line, to handle cases where
                ## the annotation spans multiple lines.
                in_host_field = False
                in_isolation_source_field = False
                line_buffer = []
                for line in gbk_fh:
                    line = line.strip()
                    ## We're going to delete all double-quote characters,
                    ## and replace all commas with semicolons so that they
                    ## don't clash with the csv format.
                    ## BUT: make sure that the line terminates with a double-quote--
                    ## otherwise, we need to slurp up data from multiple lines.
                    if line.startswith("/host"):
                        line_annot = line.split('=')[-1].replace('\"','').replace(',',';')
                        if line.endswith('"'):
                            host = line_annot
                        else: ## have to look at the next line too.
                            line_buffer.append(line_annot)
                            in_host_field = True
                    elif in_host_field and len(line_buffer):
                        line_annot = line.replace('\"','').replace(',',';')
                        line_buffer.append(line_annot)
                        if line.endswith('"'): ## flush the line buffer and reset the flag.
                            host = ' '.join(line_buffer)
                            line_buffer = []
                            in_host_field = False
                    elif line.startswith("/isolation_source"):
                        line_annot = line.split('=')[-1].replace('\"','').replace(',',';')
                        if line.endswith('"'):
                            isolation_source = line_annot
                        else: ## then have to look at next line too.
                            line_buffer.append(line_annot)
                            in_isolation_source_field = True
                    elif in_isolation_source_field and len(line_buffer):
                        line_annot = line.replace('\"','').replace(',',';')
                        line_buffer.append(line_annot)
                        if line.endswith('"'): ## flush the line buffer and reset flag.
                            isolation_source = ' '.join(line_buffer)
                            line_buffer = []
                            in_isolation_source_field = False
                    ## break if we have the annotation.
                    if (host != "NA") and (isolation_source != "NA"): break
                    ## also break if we're looking at gene annotation since
                    ## there's insufficient isolate annotation in this file.
                    if line.startswith("/gene"): break ## NOTE: has to be '/gene'
                    ## because we don't want to match the Genome Annotation Data in the
                    ## COMMENT field of the Genbank metadata.
                    if line.startswith("ORIGIN"): break ## break if we got to the first fasta seq
                row = ','.join([annotation_accession, host, isolation_source]) + '\n'
                out_fh.write(row)
    return


def filter_fastq_files_for_multireads(multiread_data_dir, themisto_pseudoalignment_dir, SRA_data_dir):

    ## make the directory for filtered multireads  if it does not exist.
    if not exists(multiread_data_dir):
        os.mkdir(multiread_data_dir)

    files_in_pseudoalignment_dir = [x for x in os.listdir(themisto_pseudoalignment_dir)]
    paths_in_pseudoalignment_dir = [os.path.join(themisto_pseudoalignment_dir, x) for x in files_in_pseudoalignment_dir]
    pseudoalignment_dirpaths = [x for x in paths_in_pseudoalignment_dir if os.path.isdir(x)]

    for my_pseudoalignment_results_dir in pseudoalignment_dirpaths:
        ## make sure to ignore the temp directory by filtering by filename suffixes!
        my_pseudoalignment_files = [x for x in os.listdir(my_pseudoalignment_results_dir) if x.endswith("_pseudoalignment.txt")]
        my_pseudoalignment_paths = [os.path.join(my_pseudoalignment_results_dir, x) for x in my_pseudoalignment_files]

        for cur_pseudoalignment_path in my_pseudoalignment_paths:
            ## before sorting the pseudoalignment results, check whether the filtered multireads
            ## for this pseudoalignment already exist on disk. If so, then skip--
            ## don't repeat work that has already been done.

            ## construct the path to the original fastq file.
            my_fastq_file = basename(cur_pseudoalignment_path).split("_pseudoalignment.txt")[0] + ".fastq"
            my_fastq_path = os.path.join(SRA_data_dir, my_fastq_file)

            ## construct the path to the filtered fastq file.
            my_genome = basename(my_pseudoalignment_results_dir)
            ## note: multiread_genome_dir is only created if filtered multireads exist.
            ## this means that this path only exists on disk if filtered multireads for this genome
            ## were already written to disk in a previous call of this function.
            multiread_genome_dir = os.path.join(multiread_data_dir, my_genome)
            my_filtered_fastq_file = "multireads_" + my_fastq_file
            my_filtered_fastq_path = os.path.join(multiread_genome_dir, my_filtered_fastq_file)
            ## now check to see if the filtered fastq file already exists.
            if exists(my_filtered_fastq_path): continue ## if so, then don't repeat the work.
            
            list_of_multiread_tuples = list()
            with open(cur_pseudoalignment_path, "r") as pseudoalignment_fh:
                for line in pseudoalignment_fh:
                    line = line.strip() ## remove whitespace
                    fields = line.split()
                    if len(fields) <= 2: continue ## skip reads that map to zero or one reference sequences.
                    read_index = int(fields[0])
                    matched_reference_tuple = tuple(fields[1:])
                    multiread_tuple = (read_index, matched_reference_tuple)
                    list_of_multiread_tuples.append(multiread_tuple)
            ## sort the list of multiread tuples in-place by the read_index.
            list_of_multiread_tuples.sort(key=lambda x: x[0])
            ## if there are multireads, then filter the original fastq file for multireads.
            if len(list_of_multiread_tuples):

                ## if there are multireads, then make multiread_genome_dir for this genome.
                if not exists(multiread_genome_dir):
                    os.mkdir(multiread_genome_dir)
                
                ## IMPORTANT: read indices in the themisto pseudoalignments are zero-based (first index is 0).
                multiread_indices = {x[0] for x in list_of_multiread_tuples} ## this is a set
                print(f"filtering {my_fastq_file} for multireads. multireads written into {my_filtered_fastq_path}")
                ## write out the filtered reads
                with open(my_filtered_fastq_path, "w") as filtered_fastq_fh:
                    my_fastq_reader = HTSeq.FastqReader(my_fastq_path)
                    for i, read in enumerate(my_fastq_reader):
                        ## skip reads that aren't in the set of multireads
                        if i not in multiread_indices: continue
                        read.write_to_fastq_file(filtered_fastq_fh)          
    return


def make_fasta_reference_genomes_for_minimap2(themisto_replicon_ref_dir):
    paths_in_themisto_replicon_ref_dir = [os.path.join(themisto_replicon_ref_dir, x) for x in os.listdir(themisto_replicon_ref_dir)]
    ## extra checking to make sure we are handling directories.
    dirpaths_in_themisto_replicon_ref_dir = [x for x in paths_in_themisto_replicon_ref_dir if os.path.isdir(x)]
    for my_dir in dirpaths_in_themisto_replicon_ref_dir:
        ## get the paths to each replicon FASTA file.
        ## IMPORTANT: the order of replicons in this file MATTERS.
        ## This corresponds to the zero-indexed replicon assignment done by themisto pseudoalignment.
        replicon_fasta_pathlist = []
        my_replicon_fasta_listfilepath = [os.path.join(my_dir, x) for x in os.listdir(my_dir) if x.endswith(".txt")].pop()
        with open(my_replicon_fasta_listfilepath, "r") as replicon_listfile_fh:
            for line in replicon_listfile_fh:
                line = line.strip()
                replicon_fasta_pathlist.append(line)

        my_fasta_reference_genome_outfile = basename(my_dir) + ".fna"
        my_fasta_reference_genome_outpath = os.path.join(my_dir, my_fasta_reference_genome_outfile)

        with open(my_fasta_reference_genome_outpath, "w") as fasta_outfh:
            for themisto_replicon_num, replicon_fasta_path in enumerate(replicon_fasta_pathlist):
                with open(replicon_fasta_path, "r") as my_fasta_infh:
                    for i, line in enumerate(my_fasta_infh):
                        ## add the replicon ID number assigned by themisto to the FASTA header
                        if i == 0:
                            header_string = line.lstrip(">")
                            replicon_id_string = "ThemistoRepliconID=" + str(themisto_replicon_num)
                            updated_header = ">" +  replicon_id_string + "|" + header_string
                            fasta_outfh.write(updated_header)
                        else:
                            fasta_outfh.write(line)
    return


def align_multireads_with_minimap2(themisto_replicon_ref_dir, multiread_data_dir, multiread_alignment_dir):
    ## Example: minimap2 -x sr ref.fa read1.fq read2.fq > aln.paf
    ## Details of PAF format are here: https://github.com/lh3/miniasm/blob/master/PAF.md

    ## make the directory for multiread alignments  if it does not exist.
    if not exists(multiread_alignment_dir):
        os.mkdir(multiread_alignment_dir)
    
    genomes_with_multireads = [x for x in os.listdir(multiread_data_dir) if x.startswith("GCF")]
    for my_genome in genomes_with_multireads:
        
        ## make a subdirectory for the multiread alignments.
        multiread_genome_alignment_dir = os.path.join(multiread_alignment_dir, my_genome)
        if not exists(multiread_genome_alignment_dir):
            os.mkdir(multiread_genome_alignment_dir)
        
        ref_genome_fasta_file = my_genome + ".fna"
        reference_genome_path = os.path.join(themisto_replicon_ref_dir, my_genome, ref_genome_fasta_file)
        my_multiread_data_dir = os.path.join(multiread_data_dir, my_genome)
        my_multiread_data_pathlist = [os.path.join(my_multiread_data_dir, x) for x in os.listdir(my_multiread_data_dir)]
        my_multiread_data_pathlist.sort() ## sort the filtered data.
        
        ## run single-end read alignment on each fastq file separately.
        for my_index, cur_multiread_data_path in enumerate(my_multiread_data_pathlist):
            my_alignment_file = basename(cur_multiread_data_path).split(".fastq")[0] + ".paf"
            my_multiread_alignment_outpath = os.path.join(multiread_genome_alignment_dir, my_alignment_file)

            minimap2_cmd_string = " ".join(["minimap2 -x sr", reference_genome_path, cur_multiread_data_path, ">", my_multiread_alignment_outpath])
            print(minimap2_cmd_string)
            subprocess.run(minimap2_cmd_string, shell=True)
    return


def parse_multiread_alignments(genome_dir):
    ## make a dictionary from reads to the multiset of replicons that the read maps to.
    paf_alignment_files = [x for x in os.listdir(genome_dir) if x.endswith(".paf")]
    paf_alignment_paths = [os.path.join(genome_dir, x) for x in paf_alignment_files]
    
    multiread_mapping_dict = dict()
    for paf_path in paf_alignment_paths:
        with open(paf_path, "r") as paf_fh:
            for line in paf_fh:
                line = line.strip()
                fields = line.split("\t")
                ## see PAF format documentation here:
                ## https://github.com/lh3/miniasm/blob/master/PAF.md
                read_name = fields[0]
                themisto_replicon_ID = int(fields[5].split("|")[0].split("=")[-1]) ## IMPORTANT: turn into an integer
                
                if read_name in multiread_mapping_dict:
                    multiread_mapping_dict[read_name].append(themisto_replicon_ID)
                else:
                    multiread_mapping_dict[read_name] = [themisto_replicon_ID]
    return multiread_mapping_dict


def make_PIRAGenomeDataFrame(
        additional_replicon_reads_dict,
        themisto_ID_to_seq_metadata_dict,
        my_naive_themisto_PCN_df):
    """ Make the DataFrame with the data needed for PIRA on a given genome.
        We have to update the results of the Naive PCN estimates from Themisto
        (results of stage 16) by adding the additional replicon reads found by re-aligning multireads
        with minimap2.
    """

    ## Turn the additional_replicon_reads_dict into a polars dataframe.
    ## I use themisto_ID_to_seq_metadata_dict to initialize the dataframe with zeros.
    ## themisto_ID_to_seq_metadata_dict contains metadata for all replicons in the genome,
    ## regardless of whether any reads mapped to that replicon.
    ## To handle the cases that:
    ## 1) a replicon only contains multireads, therefore it's naive PCN == 0, or
    ## 2) a replicon does not have any multireads, therefore, it is not present in additional_replicon_reads_dict,
    ## we use themisto_ID_to_seq_metadata_dict to initialize a polars dataframe containing
    ## all replicons in our genome.

    themisto_id_col = sorted(themisto_ID_to_seq_metadata_dict.keys())
    replicon_seq_id_col = list()
    replicon_seq_type_col = list()
    replicon_length_col = list()
    additional_readcount_col = list()
    
    for themisto_id in themisto_id_col:
        ## get the SeqID and SeqType given the Themisto ID.
        seq_id, seq_type, seq_length = themisto_ID_to_seq_metadata_dict[themisto_id]
        replicon_seq_id_col.append(seq_id)
        replicon_seq_type_col.append(seq_type)
        replicon_length_col.append(seq_length)
        ## initialize the additional_readcount_col with zeros.
        additional_readcount_col.append(0)

    ## now update the values in additional_readcount_col using additional_replicon_reads_dict.
    for themisto_id, readcount in additional_replicon_reads_dict.items():
        additional_readcount_col[int(themisto_id)] = readcount
        
    ## now make the DataFrame.
    additional_replicon_reads_df = pl.DataFrame(
        {"ThemistoID" : themisto_id_col,
         "SeqID" : replicon_seq_id_col,
         "SeqType" : replicon_seq_type_col,
         "replicon_length" : replicon_length_col,
         "AdditionalReadCount" : additional_readcount_col})

    ## IMPORTANT: the previous code ensures that additional_replicon_reads_df
    ## contains ALL replicons in the genome, even if no reads pseudoaligned or aligned to that replicon.
    ## this may happen for short contigs that are annotated as plasmids or plasmid fragments in the
    ## reference genome (example: genome GCF_017654585.1_ASM1765458v1).

    ## To fill in missing values in the AnnotationAccession column after the merge with additional_replicon_reads_df,
    ## get the unique value in the AnnotationAccession column in my_naive_themisto_PCN_df, and use this
    ## to fill in the missing values.
    my_AnnotationAccession = my_naive_themisto_PCN_df.select(pl.col('AnnotationAccession').unique())
    assert len(my_AnnotationAccession) == 1 ## make sure this value is unique!

    
    ## merge the DataFrames containing the ReadCounts.
    merged_readcount_df = additional_replicon_reads_df.join(
        ## IMPORTANT: my_naive_themisto_PCN_df will NOT contain rows for replicons that didn't have
        ## any reads pseudoalign to it. therefore, we need to left_join my_native_themisto_PCN_df to
        ## additional_replicon_reads_df, which contains the data for ALL replicons in the genome,
        ## even if the Count is zero.
        my_naive_themisto_PCN_df, on="SeqID", how="left").with_columns(
            ##fill in missing values in the AnnotationAccession column after the merge.
            pl.col("AnnotationAccession").fill_null(my_AnnotationAccession)).select(
                ## remove duplicate columns (possibly containing missing values)
                ## that came from my_naive_themisto_PCN_df.
                pl.col("*").exclude("replicon_length_right", "SeqType_right")).with_columns(
                    ## set missing values in the InitialReadCount column to 0.
                    pl.col("InitialReadCount").fill_null(strategy="zero")).with_columns(
                        ## sum those ReadCounts,
                        (pl.col("InitialReadCount") + pl.col("AdditionalReadCount")).alias("ReadCount")).with_columns(
                            ## and re-calculate SequencingCoverage,
                            (pl.col("ReadCount") / pl.col("replicon_length")).alias("SequencingCoverage"))
    
    ## The following is a hack, following code in the function native_themisto_PCN_estimation(),
    ## to recalculate LongestRepliconCoverage and CopyNumber with the additional reads.

    ## find the length of the longest replicon
    max_replicon_length = merged_readcount_df.select(pl.col("replicon_length").max()).item()

    ## filter the DataFrame to get the row for the longest replicon
    longest_replicon_row_df = merged_readcount_df.filter(
        pl.col("replicon_length") == max_replicon_length).with_columns(
            ## and define the LongestRepliconCoverage column for merging back.
            pl.col("SequencingCoverage").alias("LongestRepliconCoverage")).select(
            pl.col("AnnotationAccession", "LongestRepliconCoverage"))

    ## now normalize SequencingCoverage by LongestRepliconCoverage for each genome to calculate PCN.
    PIRAGenomeDataFrame = merged_readcount_df.join(
        ## WEIRD BEHAVIOR in polars: the AnnotationAccession key for both dataframes is preserved,
        ## I don't know why. So remove the newly made AnnotationAccession_right, SeqID_right, SeqType_right columns.
        longest_replicon_row_df, on = "AnnotationAccession").select(
        pl.col("*").exclude("AnnotationAccession_right", "SeqID_right", "SeqType_right")).with_columns(
        (pl.col("SequencingCoverage") / pl.col("LongestRepliconCoverage")).alias("InitialCopyNumberEstimate")).sort(
            ## and sort by the ThemistoID column.
            ## IMPORTANT: 
            "ThemistoID"
    )
    """
    Stage 12 calls generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir), which
    ensures that Themisto Replicon IDs are sorted by replicon length in descending order
    (i.e., Replicon ID == 0 corresponds to the longest chromosome).
    """
    return PIRAGenomeDataFrame


def initializePIRA(multiread_mapping_dict, themisto_ID_to_seq_metadata_dict, my_naive_themisto_PCN_df):
    ## Iterate over the multiread_mapping_dict, and split into two data structures:
    ## 1) reads that map to a single replicon are counted up in a dictionary.
    additional_replicon_reads_dict = dict()
    ## 2) reads that map to multiple replicons are stored in a list of lists,
    ## that will then be turned into the numpy array to store the match matrix M.
    match_matrix_list_of_rows = list()

    for read, replicon_list in multiread_mapping_dict.items():
        replicon_set = set(replicon_list)
        if len(replicon_set) == 0: ## the read does not map to any replicons -- should not occur.
            raise AssertionError(f"READ {read} DID NOT ALIGN TO ANY REPLICON")
        elif len(replicon_set) == 1: ## the read maps to a single replicon.
            my_replicon = replicon_set.pop()
            if my_replicon in additional_replicon_reads_dict:
                additional_replicon_reads_dict[my_replicon] += 1
            else:
                additional_replicon_reads_dict[my_replicon] = 1
        else: ## the read maps to multiple replicons
            ## initialize a row of zeros, based on the number of replicons in this genome.
            match_matrix_rowlist = [0 for k in themisto_ID_to_seq_metadata_dict.keys()]
            for replicon_index in replicon_list:
                match_matrix_rowlist[replicon_index] += 1
            match_matrix_list_of_rows.append(match_matrix_rowlist)

    ## now set up the data structures for PIRA on this genome.
    MatchMatrix = np.array(match_matrix_list_of_rows)
    
    """ Generate the DataFrame containing the ReadCounts,  replicon lengths, and initial PCN estimates.
    We update the results of the Naive PCN estimates from Themisto (results of stage 16)
    by adding the additional replicon reads found by re-aligning multireads
    with minimap2.
    """
    PIRAGenomeDataFrame = make_PIRAGenomeDataFrame(
        additional_replicon_reads_dict, themisto_ID_to_seq_metadata_dict, my_naive_themisto_PCN_df)
    ## return inputs requires for PIRA for this genome.
    return (MatchMatrix, PIRAGenomeDataFrame)


def run_PIRA(M, PIRAGenomeDataFrame, epsilon = 0.00001):
    ## Run PIRA for a genome, assuming that the zero-th index of the match matrix M is the chromosome for normalization.
    print("RUNNING PIRA.")
    print(M)
    print(M.shape)
    print(M.shape[0])
    print()

    print(PIRAGenomeDataFrame)
    print()
    print(PIRAGenomeDataFrame.glimpse())
    print()

    """
    Stage 12 calls generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir), which
    ensures that Themisto Replicon IDs are sorted by replicon length in descending order
    (i.e., Replicon ID == 0 corresponds to the longest chromosome).

    Therefore, we can  extract the InitialCopyNumber Column as a numpy vector as the initial PCN estimate guess for PIRA.
    """
    
    v = PIRAGenomeDataFrame["InitialCopyNumberEstimate"].to_numpy()
    readcount_vector = PIRAGenomeDataFrame["ReadCount"].to_numpy()
    replicon_length_vector = PIRAGenomeDataFrame["replicon_length"].to_numpy()
    
    """
    M may be empty, in the case that all multireads uniquely mapped to chromosome or plasmid.
    In this case, the PIRAGenomeDataFrame has already incorporated all available information,
    so we should just return v, the vector containing the InitialCopyNumberEstimate.

    I handle this logic by only running the PIRA loop if M is non-empty.
    """ 
    
    if M.shape[0] > 0: ## only run PIRA if the Match Matrix M has rows.
        convergence_error = 10000.0
        ## Iterate PIRA until error converges to zero.
        while convergence_error > epsilon:
            print(f"current convergence error: {convergence_error}")
            print(f"current PCN estimate vector: {v}")
            ## Weight M by PCN guess v-- need to turn v into a diagonal matrix first.
            diagonal_v = np.diag(v)
            weighted_M = np.matmul(M, diagonal_v)
            print(weighted_M)

            ## Normalize rows of weighted_M to sum to 1: this the probabilistic read assignment.
            ## Compute the sum of each row
            weighted_M_row_sums = weighted_M.sum(axis=1)
            ## Normalize each row by its sum
            normalized_weighted_M = weighted_M / weighted_M_row_sums[:, np.newaxis]

            print(normalized_weighted_M)
            ## sum over the rows of the normalized and weighted M matrix to generate the
            ## multiread vector.
            multiread_vector = normalized_weighted_M.sum(axis=0)
            print(multiread_vector)

            ## update the PCN estimate vector v using the multiread vector
            unnormalized_v = (multiread_vector + readcount_vector) / replicon_length_vector
            normalization_factor = unnormalized_v[0] ## coverage for the longest chromosome.
            updated_v = unnormalized_v / normalization_factor

            ## update the error estimate
            convergence_error = sum(updated_v - v)
            ## and update the PCN estimate vector v
            v = updated_v
        print(f"final convergence error: {convergence_error}")

    print(f"final PCN estimate vector: {v}")
    return v


def run_PIRA_test_suite():
    ## Test 1: a simple test to check for convergence from very bad initial PCN estimates.
    Test1_DataFrame = pl.DataFrame(
        {
            "replicon_length" : [1000000, 100000, 1000], ## log 10: 6, 5, 3
            "ReadCount" : [1000000, 1000000, 1000000], ## log 10: 6, 6, 6
            "InitialCopyNumberEstimate" : [1, 1, 1] ## should converge to: 1, 10, 1000
        }
    )

    ## Test 1: the Match Matrix is negligible.
    multi_read_row1 = [1, 1, 0]
    test1_match_matrix_list_of_rows = [multi_read_row1]
    Test1_M = np.array(test1_match_matrix_list_of_rows)
    print("*"*80)
    print("PIRA TEST 1: check convergence from bad initial PCN estimates given negligible match matrix")
    run_PIRA(Test1_M, Test1_DataFrame)
    print()


    ## Test 2-- check what happens if the initial PCN vector contains zeros-- can estimates be updated
    ## properly?
    Test2_DataFrame = pl.DataFrame(
        {
            "replicon_length" : [1000000, 100000, 1000], ## log 10: 6, 5, 3
            "ReadCount" : [1000000, 1000000, 1000000], ## log 10: 6, 6, 6
            "InitialCopyNumberEstimate" : [1, 0, 0] ## should converge to: 1, 10, 1000
        }
    )

    ## Test 2: the Match Matrix is negligible, and the initial PCN estimate contains zeros
    print("*"*80)
    print("PIRA TEST 2: check convergence from bad initial PCN estimates with zeros, given negligible match matrix")
    run_PIRA(Test1_M, Test2_DataFrame)
    print()
    
    ## Test 3: Case that Match Matrix is super important for correct PCN estimation.
    Test3_DataFrame = pl.DataFrame(
        {
            "replicon_length" : [1000000, 100000, 1000], ## log 10: 6, 5, 3
            "ReadCount" : [1000000, 0, 0], ## log 10: 6, 6, 6
            "InitialCopyNumberEstimate" : [1, 1, 1] ## should converge to: 1, 10, 1000
        }
    )

    ## Test 3: the Match Matrix is super important for accurate PCN estimation.
    multi_read_row3 = [0, 1, 1]
    test3_match_matrix_list_of_rows = [multi_read_row3] * 1000000
    Test3_M = np.array(test3_match_matrix_list_of_rows)
    print("*"*80)
    print("PIRA TEST 3: check convergence when match matrix needed for accurate PCN estimation")
    run_PIRA(Test3_M, Test3_DataFrame)
    print()
    
    return


def run_PIRA_on_all_genomes(multiread_alignment_dir, themisto_replicon_ref_dir, naive_themisto_PCN_csv_file, PIRA_PCN_csv_file):
    
    ## only run PIRA on genomes with multireads.
    genomes_with_multireads = [x for x in os.listdir(multiread_alignment_dir) if x.startswith("GCF")]

    ## trim the "_genomic" suffix from the genome directories to get the actual IDs needed
    ## for filtering the naive PCN estimates for the genomes with multireads.
    genome_IDs_with_multireads = [x.replace("_genomic", "") for x in genomes_with_multireads]
    
    ## import the results of the Naive PCN estimates from Themisto,
    ## and filter for rows corresponding to genomes with multireads.
    naive_themisto_PCN_df = pl.read_csv(naive_themisto_PCN_csv_file).filter(
        pl.col("AnnotationAccession").is_in(genome_IDs_with_multireads))

    ## Make an empty Polars DataFrame to contain all the results.
    all_PIRA_estimates_DataFrame = pl.DataFrame()
    
    ## now populate the all_PIRA_estimates_DataFrame.
    for genome in genomes_with_multireads:
        ## get the Naive PCN estimates for this particular genome.
        genome_ID = genome.replace("_genomic", "")
        ## trim the "_genomic" suffix from the genome directory to get the actual ID needed
        ## for filtering the naive PCN estimates for this genome with multireads

        
        ## FOR DEBUGGING:
        if genome_ID != "GCF_017654585.1_ASM1765458v1":
            continue
                
        my_naive_themisto_PCN_df = naive_themisto_PCN_df.filter(
            pl.col("AnnotationAccession") == genome_ID).select(
                ## select only the columns we need.
                ['AnnotationAccession', 'SeqID', 'SeqType', 'ReadCount', 'replicon_length']
            ).rename({"ReadCount": "InitialReadCount"}) ## rename ReadCount to InitialReadCount
        
        ## map the themisto replicon ID numbers to a (SeqID, SeqType) tuple.
        themisto_ID_to_seq_metadata_dict = map_themisto_IDs_to_replicon_metadata(themisto_replicon_ref_dir, genome)
        
        ## make a dictionary mapping reads to Themisto replicon IDs.
        genome_dir = os.path.join(multiread_alignment_dir, genome)
        multiread_mapping_dict = parse_multiread_alignments(genome_dir)
        
        ## initialize the data structures for PIRA.
        MatchMatrix, PIRAGenomeDataFrame = initializePIRA(
            multiread_mapping_dict, themisto_ID_to_seq_metadata_dict, my_naive_themisto_PCN_df)
        
        ## now run PIRA for this genome.
        PIRA_PCN_estimate_vector = run_PIRA(MatchMatrix, PIRAGenomeDataFrame)
        print(f"PIRA PCN estimate vector for genome {genome} is: {PIRA_PCN_estimate_vector}")
        print("*****************************************************************************************")

        ## now add the PIRA estimates as a column to the PIRAGenomeDataFrame.
        ## First convert the NumPy array to a Polars Series
        PIRA_PCN_estimate_series = pl.Series("PIRA_CopyNumberEstimate", PIRA_PCN_estimate_vector)
        ## Then add the Polars Series of PIRA estimates  to the DataFrame with initial data
        my_PIRA_PCN_estimate_DataFrame = PIRAGenomeDataFrame.with_columns([PIRA_PCN_estimate_series])

        ## now concatenate the DataFrame for this genome to the big DataFrame for all genomes.
        all_PIRA_estimates_DataFrame = pl.concat(
            [all_PIRA_estimates_DataFrame,
             my_PIRA_PCN_estimate_DataFrame])

    ## now save all_PIRA_estimates_DataFrame to disk.
    all_PIRA_estimates_DataFrame.write_csv(PIRA_PCN_csv_file)
    return





        

################################################################################

def pipeline_main():

    run_log_file = "../results/PCN-pipeline-log.txt"
    ## Configure logging
    logging.basicConfig(filename=run_log_file, level=logging.INFO)
    
    prokaryotes_with_plasmids_file = "../results/prokaryotes-with-chromosomes-and-plasmids.txt"
    RunID_table_csv = "../results/RunID_table.csv"
    reference_genome_dir = "../data/NCBI-reference-genomes/"
    SRA_data_dir = "../data/SRA/"
    ## directories for gene-level copy number estimation with kallisto.
    kallisto_gene_ref_dir = "../results/kallisto_gene_references/"
    kallisto_gene_index_dir = "../results/kallisto_gene_indices/"
    kallisto_gene_quant_results_dir = "../results/kallisto_gene_quant/"
    ## directories for replicon-level copy number estimation with kallisto.
    kallisto_replicon_ref_dir = "../results/kallisto_replicon_references/"
    kallisto_replicon_index_dir = "../results/kallisto_replicon_indices/"
    kallisto_replicon_quant_results_dir = "../results/kallisto_replicon_quant/"

    gene_copy_number_csv_file = "../results/NCBI-gene_copy_numbers.csv"
    ARG_copy_number_csv_file = "../results/NCBI-ARG_copy_numbers.csv"
    replicon_copy_number_csv_file = "../results/NCBI-replicon_copy_numbers.csv"
    calculated_copy_number_csv_file = "../results/NCBI-replicon_copy_numbers_from_genes.csv"
    replicon_length_csv_file = "../results/NCBI-replicon_lengths.csv"

    ## directories for themisto inputs and outputs.
    themisto_replicon_ref_dir = "../results/themisto_replicon_references/"
    themisto_replicon_index_dir = "../results/themisto_replicon_indices/"
    themisto_pseudoalignment_dir = "../results/themisto_replicon_pseudoalignments/"
    themisto_results_csvfile_path = "../results/themisto-replicon-read-counts.csv"

    ## this file contains estimates that throw out multireplicon reads.
    naive_themisto_PCN_csv_file = "../results/naive-themisto-PCN-estimates.csv"
    ## this file contains estimates that equally apportion multireplicon reads
    ## to the relevant plasmids and chromosomes.
    simple_themisto_PCN_csv_file = "../results/simple-themisto-PCN-estimates.csv"

    gbk_annotation_file = "../results/gbk-annotation-table.csv"

    ## directory for filtered multireads.
    multiread_data_dir = "../data/filtered_multireads/"

    ## directory for multiread alignments constructed with minimap2.
    multiread_alignment_dir = "../results/multiread_alignments/"

    ## this file contains PIRA estimates for the genomes that have multireads called by themisto.
    PIRA_PCN_csv_file = "../results/PIRA-PCN-estimates.csv"
    
    #####################################################################################
    ## Stage 1: get SRA IDs and Run IDs for all RefSeq bacterial genomes with chromosomes and plasmids.
    if exists(RunID_table_csv):
        Stage1DoneMessage = f"{RunID_table_csv} exists on disk-- skipping stage 1."
        print(Stage1DoneMessage)
        logging.info(Stage1DoneMessage)
    else: ## This takes 34513 seconds (9.5h) to get RunIDs for 4921 genomes.
        RunID_table_start_time = time.time()  # Record the start time
        create_RefSeq_SRA_RunID_table(prokaryotes_with_plasmids_file, RunID_table_csv)
        RunID_table_end_time = time.time()  # Record the end time
        RunID_table_execution_time = RunID_table_end_time - RunID_table_start_time
        Stage1TimeMessage = f"Stage 1 execution time: {RunID_table_execution_time} seconds"
        print(Stage1TimeMessage)
        logging.info(Stage1TimeMessage)
        quit()

    
    #####################################################################################
    ## Stage 2: download reference genomes for each of the bacterial genomes containing plasmids,
    ## for which we can download Illumina reads from the NCBI Short Read Archive.
    ## first, make a dictionary from RefSeq accessions to ftp paths using the
    ## prokaryotes-with-plasmids.txt file.
    ## NOTE: as of May 27 2024, 4540 reference genomes should be downloaded.
    stage_2_complete_file = "../results/stage2.done"
    reference_genome_log_file = "../results/reference_genome_fetching_log.txt"
    if exists(stage_2_complete_file):
        print(f"{stage_2_complete_file} exists on disk-- skipping stage 2.")
    else:
        refgenome_download_start_time = time.time()  ## Record the start time
        refseq_accession_to_ftp_path_dict = create_refseq_accession_to_ftp_path_dict(prokaryotes_with_plasmids_file)
        ## now download the reference genomes.
        fetch_reference_genomes(RunID_table_csv, refseq_accession_to_ftp_path_dict, reference_genome_dir, reference_genome_log_file)
        refgenome_download_end_time = time.time()  ## Record the end time
        refgenome_download_execution_time = refgenome_download_end_time - refgenome_download_start_time
        Stage2TimeMessage = f"Stage 2 (reference genome download) execution time: {refgenome_download_execution_time} seconds\n"
        with open(stage_2_complete_file, "w") as stage_2_complete_log:
            stage_2_complete_log.write(Stage2TimeMessage)
            stage_2_complete_log.write("reference genomes downloaded successfully.\n")
        quit()
    
    #####################################################################################
    ## Stage 3: download Illumina reads for the genomes from the NCBI Short Read Archive (SRA).
    stage_3_complete_file = "../results/stage3.done"
    if exists(stage_3_complete_file):
        print(f"{stage_3_complete_file} exists on disk-- skipping stage 3.")
    else:
        SRA_download_start_time = time.time()  ## Record the start time
        Run_IDs = get_Run_IDs_from_RunID_table(RunID_table_csv)
        download_fastq_reads(SRA_data_dir, Run_IDs)
        SRA_download_end_time = time.time()  ## Record the end time
        SRA_download_execution_time = SRA_download_end_time - SRA_download_start_time
        Stage3TimeMessage = f"Stage 3 (SRA download) execution time: {SRA_download_execution_time} seconds\n"
        print(Stage3TimeMessage)
        logging.info(Stage3TimeMessage)

        ## check to see if all the expected files exist on disk (does not check for corrupted data).
        if all_fastq_data_exist(Run_IDs, SRA_data_dir):
            with open(stage_3_complete_file, "w") as stage_3_complete_log:
                stage_3_complete_log.write(Stage3TimeMessage)
                stage_3_complete_log.write("SRA read data downloaded successfully.\n")
        quit()
    
    #####################################################################################   
    ## Stage 4: Make gene-level FASTA reference files for copy number estimation for genes in each genome using kallisto.
    stage_4_complete_file = "../results/stage4.done"
    if exists(stage_4_complete_file):
        print(f"{stage_4_complete_file} exists on disk-- skipping stage 4.")
    else:
        make_gene_fasta_ref_start_time = time.time()  ## Record the start time
        make_NCBI_gene_fasta_refs_for_kallisto(reference_genome_dir, kallisto_gene_ref_dir)
        make_gene_fasta_ref_end_time = time.time()  ## Record the end time
        make_gene_fasta_ref_execution_time = make_gene_fasta_ref_end_time - make_gene_fasta_ref_start_time
        Stage4TimeMessage = f"Stage 4 (making gene-level FASTA references for kallisto) execution time: {make_gene_fasta_ref_execution_time} seconds\n"
        print(Stage4TimeMessage)
        logging.info(Stage4TimeMessage)
        with open(stage_4_complete_file, "w") as stage_4_complete_log:
            stage_4_complete_log.write(Stage4TimeMessage)
            stage_4_complete_log.write("Gene-level FASTA reference sequences for kallisto finished successfully.\n")
        quit()

    ## Stage 5: Make replicon-level FASTA reference files for copy number estimation using kallisto.
    stage_5_complete_file = "../results/stage5.done"
    if exists(stage_5_complete_file):
        print(f"{stage_5_complete_file} exists on disk-- skipping stage 5.")
    else:
        make_replicon_fasta_ref_start_time = time.time()  ## Record the start time
        make_NCBI_replicon_fasta_refs_for_kallisto(reference_genome_dir, kallisto_replicon_ref_dir)
        make_replicon_fasta_ref_end_time = time.time()  ## Record the end time
        make_replicon_fasta_ref_execution_time = make_replicon_fasta_ref_end_time - make_replicon_fasta_ref_start_time
        Stage5TimeMessage = f"Stage 5 (making replicon-level FASTA references for kallisto) execution time: {make_replicon_fasta_ref_execution_time} seconds\n"

        print(Stage5TimeMessage)
        logging.info(Stage5TimeMessage)
        with open(stage_5_complete_file, "w") as stage_5_complete_log:
            stage_5_complete_log.write(Stage5TimeMessage)
            stage_5_complete_log.write("Replicon-level FASTA reference sequences for kallisto finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 6: Make gene-level kallisto index files for each genome.
    stage_6_complete_file = "../results/stage6.done"
    if exists(stage_6_complete_file):
        print(f"{stage_6_complete_file} exists on disk-- skipping stage 6.")
    else:
        make_kallisto_gene_index_start_time = time.time()  ## Record the start time
        make_NCBI_kallisto_indices(kallisto_gene_ref_dir, kallisto_gene_index_dir)
        make_kallisto_gene_index_end_time = time.time()  ## Record the end time
        make_kallisto_gene_index_execution_time = make_kallisto_gene_index_end_time - make_kallisto_gene_index_start_time
        Stage6TimeMessage = f"Stage 6 (making gene-indices for kallisto) execution time: {make_kallisto_gene_index_execution_time} seconds\n"
        print(Stage6TimeMessage)
        logging.info(Stage6TimeMessage)
        with open(stage_6_complete_file, "w") as stage_6_complete_log:
            stage_6_complete_log.write(Stage6TimeMessage)
            stage_6_complete_log.write("kallisto gene-index file construction finished successfully.\n")
        quit()

    ## Stage 7: Make replicon-level kallisto index files for each genome.
    stage_7_complete_file = "../results/stage7.done"
    if exists(stage_7_complete_file):
        print(f"{stage_7_complete_file} exists on disk-- skipping stage 7.")
    else:
        make_kallisto_replicon_index_start_time = time.time()  ## Record the start time
        make_NCBI_kallisto_indices(kallisto_replicon_ref_dir, kallisto_replicon_index_dir)
        make_kallisto_replicon_index_end_time = time.time()  ## Record the end time
        make_kallisto_replicon_index_execution_time = make_kallisto_replicon_index_end_time - make_kallisto_replicon_index_start_time
        Stage7TimeMessage = f"Stage 7 (making replicon-indices for kallisto) execution time: {make_kallisto_replicon_index_execution_time} seconds\n"
        print(Stage7TimeMessage)
        logging.info(Stage7TimeMessage)
        with open(stage_7_complete_file, "w") as stage_7_complete_log:
            stage_7_complete_log.write(Stage7TimeMessage)
            stage_7_complete_log.write("kallisto replicon-index file construction finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 8: run kallisto quant on all genome data, on both gene-level and replicon-level indices.
    ## NOTE: right now, this only processes paired-end fastq data-- single-end fastq data is ignored.
    stage_8_complete_file = "../results/stage8.done"
    if exists(stage_8_complete_file):
        print(f"{stage_8_complete_file} exists on disk-- skipping stage 8.")
    else:
        kallisto_quant_start_time = time.time()  ## Record the start time
        RefSeq_to_SRA_RunList_dict = make_RefSeq_to_SRA_RunList_dict(RunID_table_csv)
        
        run_kallisto_quant(RefSeq_to_SRA_RunList_dict, kallisto_gene_index_dir, SRA_data_dir, kallisto_gene_quant_results_dir)
        run_kallisto_quant(RefSeq_to_SRA_RunList_dict, kallisto_replicon_index_dir, SRA_data_dir, kallisto_replicon_quant_results_dir)

        kallisto_quant_end_time = time.time()  ## Record the end time
        kallisto_quant_execution_time = kallisto_quant_end_time - kallisto_quant_start_time
        Stage8TimeMessage = f"Stage 8 (kallisto quant, gene-level and replicon-level) execution time: {kallisto_quant_execution_time} seconds\n"
        print(Stage8TimeMessage)
        logging.info(Stage8TimeMessage)
        with open(stage_8_complete_file, "w") as stage_8_complete_log:
            stage_8_complete_log.write(Stage8TimeMessage)
            stage_8_complete_log.write("kallisto quant, gene-level and replicon-level finished successfully.\n")
        quit()
            
    #####################################################################################
    ## Stage 9: make a table of the estimated copy number and position for all genes in all chromosomes
    ## and plasmids in these genomes. My reasoning is that this may be useful for doing some analyses
    ## like in David Zeevi's science paper about growth rates from chromosomal copy numbers.
    stage_9_complete_file = "../results/stage9.done"
    if exists(stage_9_complete_file):
        print(f"{stage_9_complete_file} exists on disk-- skipping stage 9.")
    else:
        stage9_start_time = time.time()  ## Record the start time
        ## first make a file containing the copy number estimates for each individual gene from kallisto
        measure_NCBI_gene_copy_numbers(kallisto_gene_quant_results_dir, gene_copy_number_csv_file)
        ## then filter that output file for ARGs (faster to do in python than downstream in R).
        filter_gene_copy_number_file_for_ARGs(gene_copy_number_csv_file, ARG_copy_number_csv_file)

        stage9_end_time = time.time()  ## Record the end time
        stage9_execution_time = stage9_end_time - stage9_start_time
        Stage9TimeMessage = f"Stage 9 (tabulate all gene copy numbers) execution time: {stage9_execution_time} seconds\n"
        print(Stage9TimeMessage)
        logging.info(Stage9TimeMessage)
        with open(stage_9_complete_file, "w") as stage_9_complete_log:
            stage_9_complete_log.write(Stage9TimeMessage)
            stage_9_complete_log.write("stage 9 (tabulating all gene copy numbers) finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 10: make a table of the estimated copy number for all chromosomes and plasmids.
    stage_10_complete_file = "../results/stage10.done"
    if exists(stage_10_complete_file):
        print(f"{stage_10_complete_file} exists on disk-- skipping stage 10.")
    else:
        stage10_start_time = time.time()  ## Record the start time
        calculate_NCBI_replicon_copy_numbers_from_genes(kallisto_gene_quant_results_dir, calculated_copy_number_csv_file)
        measure_NCBI_replicon_copy_numbers(kallisto_replicon_quant_results_dir, replicon_copy_number_csv_file)

        stage10_end_time = time.time()  ## Record the end time
        stage10_execution_time = stage10_end_time - stage10_start_time
        Stage10TimeMessage = f"Stage 10 (tabulate all replicon copy numbers) execution time: {stage10_execution_time} seconds\n"
        print(Stage10TimeMessage)
        logging.info(Stage10TimeMessage)
        with open(stage_10_complete_file, "w") as stage_10_complete_log:
            stage_10_complete_log.write(Stage10TimeMessage)
            stage_10_complete_log.write("stage 10 (tabulating all replicon copy numbers) finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 11: tabulate the length of all chromosomes and plasmids.
    stage_11_complete_file = "../results/stage11.done"
    if exists(stage_11_complete_file):
        print(f"{stage_11_complete_file} exists on disk-- skipping stage 11.")
    else:
        stage11_start_time = time.time()  ## Record the start time
        tabulate_NCBI_replicon_lengths(reference_genome_dir, replicon_length_csv_file)
        stage11_end_time = time.time()  ## Record the end time
        stage11_execution_time = stage11_end_time - stage11_start_time
        Stage11TimeMessage = f"Stage 11 (tabulate all replicon lengths) execution time: {stage11_execution_time} seconds\n"
        print(Stage11TimeMessage)
        logging.info(Stage11TimeMessage)
        with open(stage_11_complete_file, "w") as stage_11_complete_log:
            stage_11_complete_log.write(Stage11TimeMessage)
            stage_11_complete_log.write("stage 11 (tabulating all replicon lengths) finished successfully.\n")
        quit()
            
    #####################################################################################
    ## Stage 12: Make FASTA input files for Themisto.
    ## Write out separate fasta files for each replicon in each genome, in a directory for each genome.
    ## Then, write out a text file that contains the paths to the FASTA files of the genomes, one file per line.
    ## See documentation here: https://github.com/algbio/themisto.
    stage_12_complete_file = "../results/stage12.done"
    if exists(stage_12_complete_file):
        print(f"{stage_12_complete_file} exists on disk-- skipping stage 12.")
    else:
        stage12_start_time = time.time()  ## Record the start time
        make_NCBI_replicon_fasta_refs_for_themisto(reference_genome_dir, themisto_replicon_ref_dir)
        stage12_end_time = time.time()  ## Record the end time
        stage12_execution_time = stage12_end_time - stage12_start_time
        Stage12TimeMessage = f"Stage 12 (making fasta references for themisto) execution time: {stage12_execution_time} seconds\n"
        print(Stage12TimeMessage)
        logging.info(Stage12TimeMessage)
        with open(stage_12_complete_file, "w") as stage_12_complete_log:
            stage_12_complete_log.write(Stage12TimeMessage)
            stage_12_complete_log.write("stage 12 (making fasta references for themisto) finished successfully.\n")
        quit()
        
    #####################################################################################
    ## Stage 13: Build separate Themisto indices for each genome.
    stage_13_complete_file = "../results/stage13.done"
    if exists(stage_13_complete_file):
        print(f"{stage_13_complete_file} exists on disk-- skipping stage 13.")
    else:
        stage13_start_time = time.time()  ## Record the start time
        make_NCBI_themisto_indices(themisto_replicon_ref_dir, themisto_replicon_index_dir)
        stage13_end_time = time.time()  ## Record the end time
        stage13_execution_time = stage13_end_time - stage13_start_time
        Stage13TimeMessage = f"Stage 13 (making indices for themisto) execution time: {stage13_execution_time} seconds\n"
        print(Stage13TimeMessage)
        logging.info(Stage13TimeMessage)
        with open(stage_13_complete_file, "w") as stage_13_complete_log:
            stage_13_complete_log.write(Stage13TimeMessage)
            stage_13_complete_log.write("stage 13 (making indices for themisto) finished successfully.\n")
        quit()
        
    #####################################################################################
    ## Stage 14: Pseudoalign reads for each genome against each Themisto index.

    stage_14_complete_file = "../results/stage14.done"
    if exists(stage_14_complete_file):
        print(f"{stage_14_complete_file} exists on disk-- skipping stage 14.")
    else:
        stage14_start_time = time.time()  ## Record the start time
        RefSeq_to_SRA_RunList_dict = make_RefSeq_to_SRA_RunList_dict(RunID_table_csv)        
        run_themisto_pseudoalign(RefSeq_to_SRA_RunList_dict, themisto_replicon_index_dir, SRA_data_dir, themisto_pseudoalignment_dir)
        stage14_end_time = time.time()  ## Record the end time
        stage14_execution_time = stage14_end_time - stage14_start_time
        Stage14TimeMessage = f"Stage 14 (themisto pseudoalignment) execution time: {stage14_execution_time} seconds\n"
        print(Stage14TimeMessage)
        logging.info(Stage14TimeMessage)
        with open(stage_14_complete_file, "w") as stage_14_complete_log:
            stage_14_complete_log.write(Stage14TimeMessage)
            stage_14_complete_log.write("stage 14 (themisto pseudoalignment) finished successfully.\n")
        quit()
        
    #####################################################################################
    ## Stage 15: generate a large CSV file summarizing the themisto pseudoalignment read counts.
    stage_15_complete_file = "../results/stage15.done"
    if exists(stage_15_complete_file):
        print(f"{stage_15_complete_file} exists on disk-- skipping stage 15.")
    else:
        stage15_start_time = time.time()  ## Record the start time
        summarize_themisto_pseudoalignment_results(themisto_replicon_ref_dir, themisto_pseudoalignment_dir, themisto_results_csvfile_path)
        stage15_end_time = time.time()  ## Record the end time
        stage15_execution_time = stage15_end_time - stage15_start_time
        Stage15TimeMessage = f"Stage 15 (themisto pseudoalignment summarization) execution time: {stage15_execution_time} seconds\n"
        print(Stage15TimeMessage)
        logging.info(Stage15TimeMessage)
        with open(stage_15_complete_file, "w") as stage_15_complete_log:
            stage_15_complete_log.write(Stage15TimeMessage)
            stage_15_complete_log.write("stage 15 (themisto pseudoalignment summarization) finished successfully.\n")
        quit()
        
    #####################################################################################
    ## Stage 16: estimate plasmid copy numbers using the themisto read counts.
    stage_16_complete_file = "../results/stage16.done"
    if exists(stage_16_complete_file):
        print(f"{stage_16_complete_file} exists on disk-- skipping stage 16.")
    else:
        stage16_start_time = time.time()  ## Record the start time
        ## Naive PCN calculation, ignoring multireplicon reads.
        naive_themisto_PCN_estimation(themisto_results_csvfile_path, replicon_length_csv_file, naive_themisto_PCN_csv_file)
        ## Simple PCN calculation, evenly distributing multireplicon reads over chromosomes and plasmids.
        simple_themisto_PCN_estimation(themisto_results_csvfile_path, replicon_length_csv_file, simple_themisto_PCN_csv_file)
        stage16_end_time = time.time()  ## Record the end time
        stage16_execution_time = stage16_end_time - stage16_start_time
        Stage16TimeMessage = f"Stage 16 (themisto PCN estimates) execution time: {stage16_execution_time} seconds\n"
        print(Stage16TimeMessage)
        logging.info(Stage16TimeMessage)
        with open(stage_16_complete_file, "w") as stage_16_complete_log:
            stage_16_complete_log.write(Stage16TimeMessage)
            stage_16_complete_log.write("stage 16 (themisto PCN estimates) finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 17: make gbk ecological annotation file.
    stage_17_complete_file = "../results/stage17.done"
    if exists(stage_17_complete_file):
        print(f"{stage_17_complete_file} exists on disk-- skipping stage 17.")
    else:
        stage17_start_time = time.time()  ## Record the start time
        make_gbk_annotation_table(reference_genome_dir, gbk_annotation_file)
        stage17_end_time = time.time()  ## Record the end time
        stage17_execution_time = stage17_end_time - stage17_start_time
        Stage17TimeMessage = f"Stage 17 (gbk ecological annotation) execution time: {stage17_execution_time} seconds\n"
        print(Stage17TimeMessage)
        logging.info(Stage17TimeMessage)
        with open(stage_17_complete_file, "w") as stage_17_complete_log:
            stage_17_complete_log.write(Stage17TimeMessage)
            stage_17_complete_log.write("stage 17 (gbk ecological annotation) finished successfully.\n")
        quit()
    
    #####################################################################################
    ## Use Probabilistic Iterative Read Assignment (PIRA) to improve PCN estimates.
    #####################################################################################
    ## The Naive PCN calculation in Stage 16 generates the initial PCN vectors and saves them on disk.
    
    #####################################################################################
    ## Stage 18: filter fastq reads for multireads.
    stage_18_complete_file = "../results/stage18.done"
    if exists(stage_18_complete_file):
        print(f"{stage_18_complete_file} exists on disk-- skipping stage 18.")
    else:
        stage18_start_time = time.time()  ## Record the start time
        filter_fastq_files_for_multireads(multiread_data_dir, themisto_pseudoalignment_dir, SRA_data_dir)
        stage18_end_time = time.time()  ## Record the end time
        stage18_execution_time = stage18_end_time - stage18_start_time
        Stage18TimeMessage = f"Stage 18 (fastq read filtering) execution time: {stage18_execution_time} seconds\n"
        print(Stage18TimeMessage)
        logging.info(Stage18TimeMessage)
        with open(stage_18_complete_file, "w") as stage_18_complete_log:
            stage_18_complete_log.write(Stage18TimeMessage)
            stage_18_complete_log.write("stage 18 (fastq read filtering) finished successfully.\n")
        quit()
    
    #####################################################################################
    ## Stage 19: make FASTA reference genomes with Themisto Replicon IDs for multiread mapping with minimap2.
    stage_19_complete_file = "../results/stage19.done"
    if exists(stage_19_complete_file):
        print(f"{stage_19_complete_file} exists on disk-- skipping stage 19.")
    else:
        stage19_start_time = time.time()  ## Record the start time
        make_fasta_reference_genomes_for_minimap2(themisto_replicon_ref_dir)
        stage19_end_time = time.time()  ## Record the end time
        stage19_execution_time = stage19_end_time - stage19_start_time
        Stage19TimeMessage = f"Stage 19 (making FASTA reference genomes for multiread alignment) execution time: {stage19_execution_time} seconds\n"
        print(Stage19TimeMessage)
        logging.info(Stage19TimeMessage)
        with open(stage_19_complete_file, "w") as stage_19_complete_log:
            stage_19_complete_log.write(Stage19TimeMessage)
            stage_19_complete_log.write("stage 19 (making FASTA reference genomes for multiread alignment) finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 20: for each genome, align multireads to the replicons with minimap2.
    stage_20_complete_file = "../results/stage20.done"
    if exists(stage_20_complete_file):
        print(f"{stage_20_complete_file} exists on disk-- skipping stage 20.")
    else:
        stage20_start_time = time.time()  ## Record the start time
        align_multireads_with_minimap2(themisto_replicon_ref_dir, multiread_data_dir, multiread_alignment_dir)
        stage20_end_time = time.time()  ## Record the end time
        stage20_execution_time = stage20_end_time - stage20_start_time
        Stage20TimeMessage = f"Stage 20 (aligning multireads with minimap2) execution time: {stage20_execution_time} seconds\n"
        print(Stage20TimeMessage)
        logging.info(Stage20TimeMessage)
        with open(stage_20_complete_file, "w") as stage_20_complete_log:
            stage_20_complete_log.write(Stage20TimeMessage)
            stage_20_complete_log.write("stage 20 (aligning multireads with minimap2) finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 21: Run PIRA.
    ## For each genome, parse minimap2 results to form the match matrix and refine the initial PCN guesses.

    stage_21_complete_file = "../results/stage21.done"
    if exists(stage_21_complete_file):
        print(f"{stage_21_complete_file} exists on disk-- skipping stage 21.")
    else:
        stage21_start_time = time.time()  ## Record the start time
        run_PIRA_on_all_genomes(multiread_alignment_dir, themisto_replicon_ref_dir, naive_themisto_PCN_csv_file, PIRA_PCN_csv_file)
        quit() ## FOR DEBUGGING
        
        stage21_end_time = time.time()  ## Record the end time
        stage21_execution_time = stage21_end_time - stage21_start_time
        Stage21TimeMessage = f"Stage 21 (parsing multiread alignments) execution time: {stage21_execution_time} seconds\n"
        print(Stage21TimeMessage)
        logging.info(Stage21TimeMessage)
        with open(stage_21_complete_file, "w") as stage_21_complete_log:
            stage_21_complete_log.write(Stage21TimeMessage)
            stage_21_complete_log.write("stage 21 (parsing multiread alignments) finished successfully.\n")
        quit()
    
    return


pipeline_main()

