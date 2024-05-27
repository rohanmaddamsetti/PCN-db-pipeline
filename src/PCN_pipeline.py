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
TODO: analyze long-read data as well, using Themisto (published 2023 in Bioinformatics).

CRITICAL TODO: utilize reads shared among replicons (multireads) to more accurately infer PCN.
In principle we can leverage reads mapping to multiple replicons (multireads) in two ways.

TODO 1) estimate PCN more accurately, using the PIRA algorithm that Lingchong and I devised.

TODO 2) we can use multireads to infer the size of duplicated/repeat regions
shared among replicons. This in itself is valuable data that can help us understand plasmid biology,
for instance, determining the extent of HGT within genomes between replicons.

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


def get_Run_IDs(sra_id):
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
        ## the Run_ID is the 3rd field from the end.
        run_accessions = [row.split("\t")[-3] for row in rows if ("Illumina") in row and ("WGS") in row]
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
            Run_IDs = get_Run_IDs(my_SRA_ID)
            for my_Run_ID in Run_IDs:
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
    refseq_ids = [line.split(",")[0] for line in RunID_table_data]
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
                    print(f"{gbff_gz_file} SUCCEEDED.", file=log_fh)
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
                        print(f"{gbff_gz_file} FAILED.", file=log_fh)
                    ## delete the corrupted files if they exist.
                    if exists(gbff_gz_file):
                        os.remove(gbff_gz_file)
                    if exists(md5_file):
                        os.remove(md5_file)
                ## if we are here, then assume the try block worked.
                if exists(gbff_gz_file) and exists(md5_file): ## then check whether the reference genome is OK.
                    if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
                        print(f"{gbff_gz_file} SUCCEEDED.", file=log_fh)
                        gbff_fetched = True  ## assume success if the checksum matches,
                        gbff_fetch_attempts = 0  ## and don't try again.
                    else:
                        os.remove(gbff_gz_file)
                        os.remove(md5_file)
    return
 

def download_fastq_reads(SRA_data_dir, RunID_table_file):
        """
        the Run_ID has to be the last part of the directory.
        see documentation here:
        https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
        """
        Run_IDs = []
        with open(RunID_table_file, "r") as RunID_table_file_obj:
            table_csv = csv.DictReader(RunID_table_file_obj)
            Run_IDs = [row["Run_ID"] for row in table_csv]
        for Run_ID in Run_IDs:
            prefetch_dir_path = os.path.join(SRA_data_dir, Run_ID)
            if os.path.exists(prefetch_dir_path): ## skip if we have already prefetched the read data.
                continue
            ## prefetch will create the prefetch_dir_path automatically-- give it the SRA_data_dir.
            prefetch_args = ["prefetch", Run_ID, "-O", SRA_data_dir]
            print (" ".join(prefetch_args))
            subprocess.run(prefetch_args)
        print("prefetch step completed.")
        my_cwd = os.getcwd()
        os.chdir(SRA_data_dir)
        for Run_ID in Run_IDs:
            sra_fastq_file_1 = Run_ID + "_1.fastq"
            sra_fastq_file_2 = Run_ID + "_2.fastq"
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


def generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir):
    genome_id = os.path.basename(fasta_outdir)
    replicon_fasta_filelist = [x for x in os.listdir(fasta_outdir) if x.endswith(".fna")]
    replicon_listfile = os.path.join(fasta_outdir, genome_id + ".txt")
    with open(replicon_listfile, "w") as fastatxtfile_fh:
        for fastafile in replicon_fasta_filelist:
            my_replicon_fasta_path = os.path.join(fasta_outdir, fastafile)
            fastatxtfile_fh.write(my_replicon_fasta_path + "\n")
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
        ## make sure that this path is real, and not an artifact of some weird file in this directory
        if not os.path.isdir(ref_fasta_dir):
            continue
        index_input_filelist = os.path.join(ref_fasta_dir, genome_id + ".txt")
        index_prefix = os.path.join(themisto_index_dir, genome_id)
        tempdir = os.path.join(themisto_index_dir, "temp")
        ## make the temp directory if it doesn't exist.
        if not exists(tempdir):
            os.mkdir(tempdir)
        themisto_build_args = ["themisto", "build", "-k","31", "-i", index_input_filelist, "--index-prefix", index_prefix, "--temp-dir", tempdir, "--mem-gigas", "2", "--n-threads", "4", "--file-colors"]
        themisto_build_string = " ".join(themisto_build_args)
        slurm_string = "sbatch -p scavenger --mem=2G --cpus-per-task=4 --wrap=" + "\"" + themisto_build_string + "\""
        print(slurm_string)
        subprocess.run(slurm_string, shell=True)
    return


def run_themisto_pseudoalign(RefSeq_to_SRA_RunList_dict, themisto_index_dir, SRA_data_dir, themisto_pseudoalignment_dir):
    ## make the output directory if it does not exist.
    if not exists(themisto_pseudoalignment_dir):
        os.mkdir(themisto_pseudoalignment_dir)

    ## each directory in themisto_index_dir is named after the genome_id of the given genome.
    for genome_id in os.listdir(themisto_index_dir):
        my_index_dir = os.path.join(themisto_index_dir, genome_id)
        ## make sure that this path is real, and not an artifact of some weird file in this directory
        if not os.path.isdir(my_index_dir):
            continue
        ## make the arguments relating to the paths to the index files for this genome.
        my_index_prefix = os.path.join(themisto_index_dir, genome_id)

        ## make the file containing the paths to the sequencing read data for this genome.
        refseq_id = "_".join(genome_id.split("_")[:2])
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

        tempdir = os.path.join(themisto_pseudoalignment_dir, "temp")
        ## make the temp directory if it doesn't exist.
        if not exists(tempdir):
            os.mkdir(tempdir)

        ## make the directory for the pseudoalignments.
        my_pseudoalignment_output_dir = os.path.join(themisto_pseudoalignment_dir, genome_id)
        if not exists(my_pseudoalignment_output_dir):
            os.mkdir(my_pseudoalignment_output_dir)
            
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
        slurm_string = "sbatch -p scavenger --mem=2G --cpus-per-task=4 --wrap=" + "\"" + themisto_pseudoalign_string + "\""
        print(slurm_string)
        subprocess.run(slurm_string, shell=True)
    return


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
                        replicon_set_string = " ".join(sorted(line.strip().split()[1:]))
                        if replicon_set_string in pseudoalignment_read_count_dict:
                            pseudoalignment_read_count_dict[replicon_set_string] += 1
                        else:
                            pseudoalignment_read_count_dict[replicon_set_string] = 1
            ## now, let's map the themisto replicon ID numbers to a (SeqID, SeqType) tuple.
            themisto_ID_to_seq_metadata_dict = dict()
            my_themisto_replicon_ID_mapping_file = os.path.join(themisto_replicon_ref_dir, my_genome_dirname, my_genome_dirname + ".txt")
            with open(my_themisto_replicon_ID_mapping_file, "r") as replicon_ID_mapping_fh:
                for i,fasta_file_path in enumerate(replicon_ID_mapping_fh):
                    fasta_file_path = fasta_file_path.strip()
                    ## get the header from this fasta file.
                    with open(fasta_file_path, "r") as fasta_fh:
                        my_header = fasta_fh.readline().strip()
                    fields = my_header.strip(">").split("|")
                    SeqID = fields[0].split("=")[-1]
                    SeqType = fields[1].split("=")[-1]
                    themisto_ID_to_seq_metadata_dict[str(i)] = (SeqID, SeqType)
            ## now write the pseudoalignment counts to file.
            for replicon_set_string in sorted(pseudoalignment_read_count_dict.keys()):
                read_count = pseudoalignment_read_count_dict[replicon_set_string]
                replicon_ID_list = replicon_set_string.split()
                if len(replicon_ID_list) == 0:
                    SeqID = "NA"
                    SeqType = "NA"
                elif len(replicon_ID_list) == 1:
                    my_replicon_ID = replicon_ID_list[0]
                    SeqID, SeqType = themisto_ID_to_seq_metadata_dict[my_replicon_ID]
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


def filter_fastq_files_for_multireads():
    pass


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
    themisto_replicon_index_dir = "../results/themisto_replicon_indices"
    themisto_pseudoalignment_dir = "../results/themisto_replicon_pseudoalignments"
    themisto_results_csvfile_path = "../results/themisto-replicon-read-counts.csv"

    ## this file contains estimates that throw out multireplicon reads.
    naive_themisto_PCN_csv_file = "../results/naive-themisto-PCN-estimates.csv"
    ## this file contains estimates that equally apportion multireplicon reads
    ## to the relevant plasmids and chromosomes.
    simple_themisto_PCN_csv_file = "../results/simple-themisto-PCN-estimates.csv"

    gbk_annotation_file = "../results/gbk-annotation-table.csv"
    

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
        refseq_accession_to_ftp_path_dict = create_refseq_accession_to_ftp_path_dict(prokaryotes_with_plasmids_file)
        ## now download the reference genomes.
        fetch_reference_genomes(RunID_table_csv, refseq_accession_to_ftp_path_dict, reference_genome_dir, reference_genome_log_file)
        with open(stage_2_complete_file, "w") as stage_2_complete_log:
            stage_2_complete_log.write("reference genomes downloaded successfully.\n")
        quit()
    
    #####################################################################################
    ## Stage 3: download Illumina reads for the genomes from the NCBI Short Read Archive (SRA).
    stage_3_complete_file = "../results/stage3.done"
    if exists(stage_3_complete_file):
        print(f"{stage_3_complete_file} exists on disk-- skipping stage 3.")
    else:
        SRA_download_start_time = time.time()  # Record the start time
        download_fastq_reads(SRA_data_dir, RunID_table_csv)
        SRA_download_end_time = time.time()  # Record the end time
        SRA_download_execution_time = SRA_download_end_time - SRA_download_start_time
        Stage3TimeMessage = f"Stage 3 (SRA download) execution time: {SRA_download_execution_time} seconds"
        print(Stage3TimeMessage)
        logging.info(Stage3TimeMessage)
        with open(stage_3_complete_file, "w") as stage_3_complete_log:
            stage_3_complete_log.write("SRA read data downloaded successfully.\n")
        quit()
    
    #####################################################################################   
    ## Stage 4: Make gene-level FASTA reference files for copy number estimation for genes in each genome using kallisto.
    stage_4_complete_file = "../results/stage4.done"
    if exists(stage_4_complete_file):
        print(f"{stage_4_complete_file} exists on disk-- skipping stage 4.")
    else:
        make_gene_fasta_ref_start_time = time.time()  # Record the start time
        make_NCBI_gene_fasta_refs_for_kallisto(reference_genome_dir, kallisto_gene_ref_dir)
        make_gene_fasta_ref_end_time = time.time()  # Record the end time
        make_gene_fasta_ref_execution_time = make_gene_fasta_ref_end_time - make_gene_fasta_ref_start_time
        Stage4TimeMessage = f"Stage 4 (making gene-level FASTA references for kallisto) execution time: {make_gene_fasta_ref_execution_time} seconds"
        print(Stage4TimeMessage)
        logging.info(Stage4TimeMessage)
        with open(stage_4_complete_file, "w") as stage_4_complete_log:
            stage_4_complete_log.write("Gene-level FASTA reference sequences for kallisto finished successfully.\n")
        quit()

    ## Stage 5: Make replicon-level FASTA reference files for copy number estimation using kallisto.
    stage_5_complete_file = "../results/stage5.done"
    if exists(stage_5_complete_file):
        print(f"{stage_5_complete_file} exists on disk-- skipping stage 5.")
    else:
        make_replicon_fasta_ref_start_time = time.time()  # Record the start time
        make_NCBI_replicon_fasta_refs_for_kallisto(reference_genome_dir, kallisto_replicon_ref_dir)
        make_replicon_fasta_ref_end_time = time.time()  # Record the end time
        make_replicon_fasta_ref_execution_time = make_replicon_fasta_ref_end_time - make_replicon_fasta_ref_start_time
        Stage5TimeMessage = f"Stage 5 (making replicon-level FASTA references for kallisto) execution time: {make_replicon_fasta_ref_execution_time} seconds"

        print(Stage5TimeMessage)
        logging.info(Stage5TimeMessage)
        with open(stage_5_complete_file, "w") as stage_5_complete_log:
            stage_5_complete_log.write("Replicon-level FASTA reference sequences for kallisto finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 6: Make gene-level kallisto index files for each genome.
    stage_6_complete_file = "../results/stage6.done"
    if exists(stage_6_complete_file):
        print(f"{stage_6_complete_file} exists on disk-- skipping stage 6.")
    else:
        make_kallisto_gene_index_start_time = time.time()  # Record the start time
        make_NCBI_kallisto_indices(kallisto_gene_ref_dir, kallisto_gene_index_dir)
        make_kallisto_gene_index_end_time = time.time()  # Record the end time
        make_kallisto_gene_index_execution_time = make_kallisto_gene_index_end_time - make_kallisto_gene_index_start_time
        Stage6TimeMessage = f"Stage 6 (making gene-indices for kallisto) execution time: {make_kallisto_gene_index_execution_time} seconds"
        print(Stage6TimeMessage)
        logging.info(Stage6TimeMessage)
        with open(stage_6_complete_file, "w") as stage_6_complete_log:
            stage_6_complete_log.write("kallisto gene-index file construction finished successfully.\n")
        quit()

    ## Stage 7: Make replicon-level kallisto index files for each genome.
    stage_7_complete_file = "../results/stage7.done"
    if exists(stage_7_complete_file):
        print(f"{stage_7_complete_file} exists on disk-- skipping stage 7.")
    else:
        make_kallisto_replicon_index_start_time = time.time()  # Record the start time
        make_NCBI_kallisto_indices(kallisto_replicon_ref_dir, kallisto_replicon_index_dir)
        make_kallisto_replicon_index_end_time = time.time()  # Record the end time
        make_kallisto_replicon_index_execution_time = make_kallisto_replicon_index_end_time - make_kallisto_replicon_index_start_time
        Stage7TimeMessage = f"Stage 7 (making replicon-indices for kallisto) execution time: {make_kallisto_replicon_index_execution_time} seconds"
        print(Stage7TimeMessage)
        logging.info(Stage7TimeMessage)
        with open(stage_7_complete_file, "w") as stage_7_complete_log:
            stage_7_complete_log.write("kallisto replicon-index file construction finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 8: run kallisto quant on all genome data, on both gene-level and replicon-level indices.
    ## NOTE: right now, this only processes paired-end fastq data-- single-end fastq data is ignored.
    stage_8_complete_file = "../results/stage8.done"
    if exists(stage_8_complete_file):
        print(f"{stage_8_complete_file} exists on disk-- skipping stage 8.")
    else:
        kallisto_quant_start_time = time.time()  # Record the start time
        RefSeq_to_SRA_RunList_dict = make_RefSeq_to_SRA_RunList_dict(RunID_table_csv)
        
        run_kallisto_quant(RefSeq_to_SRA_RunList_dict, kallisto_gene_index_dir, SRA_data_dir, kallisto_gene_quant_results_dir)
        run_kallisto_quant(RefSeq_to_SRA_RunList_dict, kallisto_replicon_index_dir, SRA_data_dir, kallisto_replicon_quant_results_dir)

        kallisto_quant_end_time = time.time()  # Record the end time
        kallisto_quant_execution_time = kallisto_quant_end_time - kallisto_quant_start_time
        Stage8TimeMessage = f"Stage 8 (kallisto quant, gene-level and replicon-level) execution time: {kallisto_quant_execution_time} seconds"
        print(Stage8TimeMessage)
        logging.info(Stage8TimeMessage)
        with open(stage_8_complete_file, "w") as stage_8_complete_log:
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
        stage9_start_time = time.time()  # Record the start time
        ## first make a file containing the copy number estimates for each individual gene from kallisto
        measure_NCBI_gene_copy_numbers(kallisto_gene_quant_results_dir, gene_copy_number_csv_file)
        ## then filter that output file for ARGs (faster to do in python than downstream in R).
        filter_gene_copy_number_file_for_ARGs(gene_copy_number_csv_file, ARG_copy_number_csv_file)

        stage9_end_time = time.time()  # Record the end time
        stage9_execution_time = stage9_end_time - stage9_start_time
        Stage9TimeMessage = f"Stage 9 (tabulate all gene copy numbers) execution time: {stage9_execution_time} seconds"
        print(Stage9TimeMessage)
        logging.info(Stage9TimeMessage)
        with open(stage_9_complete_file, "w") as stage_9_complete_log:
            stage_9_complete_log.write("stage 9 (tabulating all gene copy numbers) finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 10: make a table of the estimated copy number for all chromosomes and plasmids.
    stage_10_complete_file = "../results/stage10.done"
    if exists(stage_10_complete_file):
        print(f"{stage_10_complete_file} exists on disk-- skipping stage 10.")
    else:
        stage10_start_time = time.time()  # Record the start time
        calculate_NCBI_replicon_copy_numbers_from_genes(kallisto_gene_quant_results_dir, calculated_copy_number_csv_file)
        measure_NCBI_replicon_copy_numbers(kallisto_replicon_quant_results_dir, replicon_copy_number_csv_file)

        stage10_end_time = time.time()  # Record the end time
        stage10_execution_time = stage10_end_time - stage10_start_time
        Stage10TimeMessage = f"Stage 10 (tabulate all replicon copy numbers) execution time: {stage10_execution_time} seconds"
        print(Stage10TimeMessage)
        logging.info(Stage10TimeMessage)
        with open(stage_10_complete_file, "w") as stage_10_complete_log:
            stage_10_complete_log.write("stage 10 (tabulating all replicon copy numbers) finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 11: tabulate the length of all chromosomes and plasmids.
    stage_11_complete_file = "../results/stage11.done"
    if exists(stage_11_complete_file):
        print(f"{stage_11_complete_file} exists on disk-- skipping stage 11.")
    else:
        stage11_start_time = time.time()  # Record the start time
        tabulate_NCBI_replicon_lengths(reference_genome_dir, replicon_length_csv_file)
        stage11_end_time = time.time()  # Record the end time
        stage11_execution_time = stage11_end_time - stage11_start_time
        Stage11TimeMessage = f"Stage 11 (tabulate all replicon lengths) execution time: {stage11_execution_time} seconds"
        print(Stage11TimeMessage)
        logging.info(Stage11TimeMessage)
        with open(stage_11_complete_file, "w") as stage_11_complete_log:
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
        stage12_start_time = time.time()  # Record the start time
        make_NCBI_replicon_fasta_refs_for_themisto(reference_genome_dir, themisto_replicon_ref_dir)
        stage12_end_time = time.time()  # Record the end time
        stage12_execution_time = stage12_end_time - stage12_start_time
        Stage12TimeMessage = f"Stage 12 (making fasta references for themisto) execution time: {stage12_execution_time} seconds"
        print(Stage12TimeMessage)
        logging.info(Stage12TimeMessage)
        with open(stage_12_complete_file, "w") as stage_12_complete_log:
            stage_12_complete_log.write("stage 12 (making fasta references for themisto) finished successfully.\n")
        quit()
        
    #####################################################################################
    ## Stage 13: Build separate Themisto indices for each genome.
    stage_13_complete_file = "../results/stage13.done"
    if exists(stage_13_complete_file):
        print(f"{stage_13_complete_file} exists on disk-- skipping stage 13.")
    else:
        stage13_start_time = time.time()  # Record the start time
        make_NCBI_themisto_indices(themisto_replicon_ref_dir, themisto_replicon_index_dir)
        stage13_end_time = time.time()  # Record the end time
        stage13_execution_time = stage13_end_time - stage13_start_time
        Stage13TimeMessage = f"Stage 13 (making indices for themisto) execution time: {stage13_execution_time} seconds"
        print(Stage13TimeMessage)
        logging.info(Stage13TimeMessage)
        with open(stage_13_complete_file, "w") as stage_13_complete_log:
            stage_13_complete_log.write("stage 13 (making indices for themisto) finished successfully.\n")
        quit()
        
    #####################################################################################
    ## Stage 14: Pseudoalign reads for each genome against each Themisto index.

    stage_14_complete_file = "../results/stage14.done"
    if exists(stage_14_complete_file):
        print(f"{stage_14_complete_file} exists on disk-- skipping stage 14.")
    else:
        stage14_start_time = time.time()  # Record the start time
        RefSeq_to_SRA_RunList_dict = make_RefSeq_to_SRA_RunList_dict(RunID_table_csv)
        run_themisto_pseudoalign(RefSeq_to_SRA_RunList_dict, themisto_replicon_index_dir, SRA_data_dir, themisto_pseudoalignment_dir)
        stage14_end_time = time.time()  # Record the end time
        stage14_execution_time = stage14_end_time - stage14_start_time
        Stage14TimeMessage = f"Stage 14 (themisto pseudoalignment) execution time: {stage14_execution_time} seconds"
        print(Stage14TimeMessage)
        logging.info(Stage14TimeMessage)
        with open(stage_14_complete_file, "w") as stage_14_complete_log:
            stage_14_complete_log.write("stage 14 (themisto pseudoalignment) finished successfully.\n")
        quit()
        
    #####################################################################################
    ## Stage 15: generate a large CSV file summarizing the themisto pseudoalignment read counts.
    stage_15_complete_file = "../results/stage15.done"
    if exists(stage_15_complete_file):
        print(f"{stage_15_complete_file} exists on disk-- skipping stage 15.")
    else:
        stage15_start_time = time.time()  # Record the start time
        summarize_themisto_pseudoalignment_results(themisto_replicon_ref_dir, themisto_pseudoalignment_dir, themisto_results_csvfile_path)
        stage15_end_time = time.time()  # Record the end time
        stage15_execution_time = stage15_end_time - stage15_start_time
        Stage15TimeMessage = f"Stage 15 (themisto pseudoalignment summarization) execution time: {stage15_execution_time} seconds"
        print(Stage15TimeMessage)
        logging.info(Stage15TimeMessage)
        with open(stage_15_complete_file, "w") as stage_15_complete_log:
            stage_15_complete_log.write("stage 15 (themisto pseudoalignment summarization) finished successfully.\n")
        quit()
        
    #####################################################################################
    ## Stage 16: estimate plasmid copy numbers using the themisto read counts.
    stage_16_complete_file = "../results/stage16.done"
    if exists(stage_16_complete_file):
        print(f"{stage_16_complete_file} exists on disk-- skipping stage 16.")
    else:
        stage16_start_time = time.time()  # Record the start time
        ## Naive PCN calculation, ignoring multireplicon reads.
        naive_themisto_PCN_estimation(themisto_results_csvfile_path, replicon_length_csv_file, naive_themisto_PCN_csv_file)
        ## Simple PCN calculation, evenly distributing multireplicon reads over chromosomes and plasmids.
        simple_themisto_PCN_estimation(themisto_results_csvfile_path, replicon_length_csv_file, simple_themisto_PCN_csv_file)        
        stage16_end_time = time.time()  # Record the end time
        stage16_execution_time = stage16_end_time - stage16_start_time
        Stage16TimeMessage = f"Stage 16 (themisto PCN estimates) execution time: {stage16_execution_time} seconds"
        print(Stage16TimeMessage)
        logging.info(Stage16TimeMessage)
        with open(stage_16_complete_file, "w") as stage_16_complete_log:
            stage_16_complete_log.write("stage 16 (themisto PCN estimates) finished successfully.\n")
        quit()

    #####################################################################################
    ## Stage 17: make gbk ecological annotation file.
    stage_17_complete_file = "../results/stage17.done"
    if exists(stage_17_complete_file):
        print(f"{stage_17_complete_file} exists on disk-- skipping stage 17.")
    else:
        stage17_start_time = time.time()  # Record the start time
        make_gbk_annotation_table(reference_genome_dir, gbk_annotation_file)
        stage17_end_time = time.time()  # Record the end time
        stage17_execution_time = stage17_end_time - stage17_start_time
        Stage17TimeMessage = f"Stage 17 (gbk ecological annotation) execution time: {stage17_execution_time} seconds"
        print(Stage17TimeMessage)
        logging.info(Stage17TimeMessage)
        with open(stage_17_complete_file, "w") as stage_17_complete_log:
            stage_17_complete_log.write("stage 17 (gbk ecological annotation) finished successfully.\n")
        quit()

    #####################################################################################
    ## Use Probabilistic Iterative Read Assignment (PIRA) to improve PCN estimates.
    #####################################################################################
    ## Stage 18: filter fastq reads for multireads.
    stage_18_complete_file = "../results/stage18.done"
    if exists(stage_18_complete_file):
        print(f"{stage_18_complete_file} exists on disk-- skipping stage 18.")
    else:
        stage18_start_time = time.time()  # Record the start time

        filter_fastq_files_for_multireads()
        quit()

        
        stage18_end_time = time.time()  # Record the end time
        stage18_execution_time = stage18_end_time - stage18_start_time
        Stage18TimeMessage = f"Stage 18 (fastq read filtering) execution time: {stage18_execution_time} seconds"
        print(Stage18TimeMessage)
        logging.info(Stage18TimeMessage)
        with open(stage_18_complete_file, "w") as stage_18_complete_log:
            stage_18_complete_log.write("stage 18 (fastq read filtering) finished successfully.\n")
        quit()

    
    return


pipeline_main()


