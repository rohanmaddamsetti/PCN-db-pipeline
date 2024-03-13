#!/usr/bin/env python

"""
PCN_pipeline.py by Maggie Wilson and Rohan Maddamsetti.

For this pipeline to work, ncbi datasets, pysradb, and kallisto must be in the $PATH.
On DCC, run the following to get these programs into the path:
conda activate PCNdb-env

Currently, this pipeline only analyzes Illumina short-read data with kallisto.

In the future, one could additionally analyze long-read data as well, by using the
state-of-the-art in pseudoalignment analyses, Themisto (published 2023 in Bioinformatics).

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
    my_md5_checksum = md5_call.stdout.split("=")[-1].strip()
    ## verify that the checksums match.
    return my_md5_checksum == my_target_checksum


def fetch_reference_genomes(RunID_table_file, refseq_accession_to_ftp_path_dict, reference_genome_dir):
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

    for ftp_path in ftp_paths:
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
                ## delete the corrupted files if they exist.
                if exists(gbff_gz_file):
                    os.remove(gbff_gz_file)
                if exists(md5_file):
                    os.remove(md5_file)
            ## if we are here, then assume the try block worked.
            if exists(gbff_gz_file) and exists(md5_file): ## then check whether the reference genome is OK.
                if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
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


def generate_fasta_reference_for_kallisto(gbk_gz_path, outfile):
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


def make_NCBI_fasta_refs_for_kallisto(refgenomes_dir, kallisto_ref_outdir):    
    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff.gz")]
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        genome_id = gzfile.split(".gbff.gz")[0]
        fasta_outfile = os.path.join(kallisto_ref_outdir, genome_id+".fna")
        generate_fasta_reference_for_kallisto(gzpath, fasta_outfile)
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
            slurm_string = "sbatch -p youlab --mem=16G --wrap=" + "\"" + kallisto_quant_string + "\""
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


def parse_metadata_in_header(target_id):
    fields = target_id.split("|")
    SeqID = fields[0].split("=")[-1]
    SeqType = fields[1].split("=")[-1]
    locus_tag = fields[2].split("=")[-1]
    ## convert underscores back into spaces.
    product = fields[3].split("=")[-1].replace("_", " ")
    metadata_tuple = (SeqID, SeqType, locus_tag, product)
    return(metadata_tuple)


def estimate_chr_plasmid_copy_numbers(genecount_tsv_path):
    genome_dict = dict()
    ## keys are SeqIDs.
    ## values are a dict: {SeqType: "chromosome", total_length: 10000, total_est_counts: 100}
    with open(genecount_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, locus_tag, product = parse_metadata_in_header(target_id)
            if SeqID in genome_dict:
                genome_dict[SeqID]["total_length"] += float(length)
                genome_dict[SeqID]["total_est_counts"] += float(est_counts)
            else: ## Initialize the dictionary.
                genome_dict[SeqID] = {"SeqType" : SeqType, "total_length" : float(length), "total_est_counts": float(est_counts)}
    coverage_dict = dict()
    print(genome_dict)
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
    print(copy_number_dict)
    return(copy_number_dict)


def measure_NCBI_replicon_copy_numbers(kallisto_quant_results_dir, copy_number_csv_file):
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
        copy_number_dict = estimate_chr_plasmid_copy_numbers(genome_quantfile_path)
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
            SeqID, SeqType, locus_tag, product = parse_metadata_in_header(target_id)
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


def measure_NCBI_gene_copy_numbers(kallisto_quant_results_dir, gene_copy_number_csv_file):
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
    genomedirectories = [x for x in os.listdir(kallisto_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        refseq_id = "_".join(genomedir.split("_")[:2])
        genome_quantfile_path = os.path.join(kallisto_quant_results_dir, genomedir, "abundance.tsv")
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


def tabulate_NCBI_replicon_lengths(refgenomes_dir, replicon_length_csv_file):
    with open(replicon_length_csv_file, 'w') as outfh:
        header = "AnnotationAccession,SeqID,replicon_length\n"
        outfh.write(header)
        for gbk_gz in os.listdir(refgenomes_dir):
            if not gbk_gz.endswith(".gbff.gz"): continue
            annotation_accession = gbk_gz.split("_genomic.gbff")[0]
            infile = os.path.join(refgenomes_dir, gbk_gz)
            with gzip.open(infile, "rt") as genome_fh:
                for replicon in SeqIO.parse(genome_fh, "gb"):
                    SeqID = replicon.id
                    replicon_length = str(len(replicon))
                    ## now write out the data for the replicon.
                    row = ','.join([annotation_accession, SeqID, replicon_length])
                    outfh.write(row + "\n")
    return


################################################################################

def pipeline_main():

    run_log_file = "../results/PCN-pipeline-log.txt"
    ## Configure logging
    logging.basicConfig(filename=run_log_file, level=logging.INFO)
    
    prokaryotes_with_plasmids_file = "../results/prokaryotes-with-plasmids.txt"
    RunID_table_csv = "../results/RunID_table.csv"
    reference_genome_dir = "../data/NCBI-reference-genomes/"
    SRA_data_dir = "../data/SRA/"
    
    kallisto_ref_dir = "../results/kallisto_references/"
    kallisto_index_dir = "../results/kallisto_indices/"
    kallisto_quant_results_dir = "../results/kallisto_quant/"

    gene_copy_number_csv_file = "../results/NCBI-gene_copy_numbers.csv"
    copy_number_csv_file = "../results/NCBI-chromosome_plasmid_copy_numbers.csv"
    replicon_length_csv_file = "../results/NCBI-replicon_lengths.csv"

    ## Stage 1: get SRA IDs and Run IDs for all complete RefSeq bacterial genomes with plasmids.
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

    ## Stage 2: download reference genomes for each of the bacterial genomes containing plasmids,
    ## for which we can download Illumina reads from the NCBI Short Read Archive.
    ## first, make a dictionary from RefSeq accessions to ftp paths using the
    ## prokaryotes-with-plasmids.txt file.
    stage_2_complete_file = "../results/stage2.done"
    if exists(stage_2_complete_file):
        print(f"{stage_2_complete_file} exists on disk-- skipping stage 2.")
    else:
        refseq_accession_to_ftp_path_dict = create_refseq_accession_to_ftp_path_dict(prokaryotes_with_plasmids_file)
        ## now download the reference genomes.
        fetch_reference_genomes(RunID_table_csv, refseq_accession_to_ftp_path_dict, reference_genome_dir)
        with open(stage_2_complete_file, "w") as stage_2_complete_log:
            stage_2_complete_log.write("reference genomes downloaded successfully.\n")
        
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

    ## Stage 4: Make FASTA reference files for copy number estimation for genes in each genome using kallisto.
    stage_4_complete_file = "../results/stage4.done"
    if exists(stage_4_complete_file):
        print(f"{stage_4_complete_file} exists on disk-- skipping stage 4.")
    else:
        make_fasta_ref_start_time = time.time()  # Record the start time
        make_NCBI_fasta_refs_for_kallisto(reference_genome_dir, kallisto_ref_dir)
        make_fasta_ref_end_time = time.time()  # Record the end time
        make_fasta_ref_execution_time = make_fasta_ref_end_time - make_fasta_ref_start_time
        Stage4TimeMessage = f"Stage 4 (making FASTA references for kallisto) execution time: {make_fasta_ref_execution_time} seconds"
        print(Stage4TimeMessage)
        logging.info(Stage4TimeMessage)
        with open(stage_4_complete_file, "w") as stage_4_complete_log:
            stage_4_complete_log.write("FASTA reference sequences for kallisto finished successfully.\n")

    ## Stage 5: Make kallisto index files for each genome.
    stage_5_complete_file = "../results/stage5.done"
    if exists(stage_5_complete_file):
        print(f"{stage_5_complete_file} exists on disk-- skipping stage 5.")
    else:
        make_kallisto_index_start_time = time.time()  # Record the start time
        make_NCBI_kallisto_indices(kallisto_ref_dir, kallisto_index_dir)
        make_kallisto_index_end_time = time.time()  # Record the end time
        make_kallisto_index_execution_time = make_kallisto_index_end_time - make_kallisto_index_start_time
        Stage5TimeMessage = f"Stage 5 (making indices for kallisto) execution time: {make_kallisto_index_execution_time} seconds"
        print(Stage5TimeMessage)
        logging.info(Stage5TimeMessage)
        with open(stage_5_complete_file, "w") as stage_5_complete_log:
            stage_5_complete_log.write("kallisto index file construction finished successfully.\n")

    ## Stage 6: run kallisto quant on all genome data.
    ## NOTE: right now, this only processes paired-end fastq data-- single-end fastq data is ignored.
    stage_6_complete_file = "../results/stage6.done"
    if exists(stage_6_complete_file):
        print(f"{stage_6_complete_file} exists on disk-- skipping stage 6.")
    else:
        kallisto_quant_start_time = time.time()  # Record the start time
        RefSeq_to_SRA_RunList_dict = make_RefSeq_to_SRA_RunList_dict(RunID_table_csv)
        run_kallisto_quant(RefSeq_to_SRA_RunList_dict, kallisto_index_dir, SRA_data_dir, kallisto_quant_results_dir)
        kallisto_quant_end_time = time.time()  # Record the end time
        kallisto_quant_execution_time = kallisto_quant_end_time - kallisto_quant_start_time
        Stage6TimeMessage = f"Stage 6 (kallisto quant) execution time: {kallisto_quant_execution_time} seconds"
        print(Stage6TimeMessage)
        logging.info(Stage6TimeMessage)
        with open(stage_6_complete_file, "w") as stage_6_complete_log:
            stage_6_complete_log.write("kallisto quant finished successfully.\n")

    ## Stage 7: make a table of the estimated copy number and position for all genes in all chromosomes
    ## and plasmids in these genomes. My reasoning is that this may be useful for doing some analyses
    ## like in David Zeevi's science paper about growth rates from chromosomal copy numbers.
    stage_7_complete_file = "../results/stage7.done"
    if exists(stage_7_complete_file):
        print(f"{stage_7_complete_file} exists on disk-- skipping stage 7.")
    else:
        stage7_start_time = time.time()  # Record the start time

        measure_NCBI_gene_copy_numbers(kallisto_quant_results_dir, gene_copy_number_csv_file)

        stage7_end_time = time.time()  # Record the end time
        stage7_execution_time = stage7_end_time - stage7_start_time
        Stage7TimeMessage = f"Stage 7 (tabulate all gene copy numbers) execution time: {stage7_execution_time} seconds"
        print(Stage7TimeMessage)
        logging.info(Stage7TimeMessage)
        with open(stage_7_complete_file, "w") as stage_7_complete_log:
            stage_7_complete_log.write("stage 7 (tabulating all gene copy numbers) finished successfully.\n")

    ## Stage 8: make a table of the estimated copy number for all chromosomes and plasmids.
    stage_8_complete_file = "../results/stage8.done"
    if exists(stage_8_complete_file):
        print(f"{stage_8_complete_file} exists on disk-- skipping stage 8.")
    else:
        stage8_start_time = time.time()  # Record the start time
        
        measure_NCBI_replicon_copy_numbers(kallisto_quant_results_dir, copy_number_csv_file)

        stage8_end_time = time.time()  # Record the end time
        stage8_execution_time = stage8_end_time - stage8_start_time
        Stage8TimeMessage = f"Stage 8 (tabulate all replicon copy numbers) execution time: {stage8_execution_time} seconds"
        print(Stage8TimeMessage)
        logging.info(Stage8TimeMessage)
        with open(stage_8_complete_file, "w") as stage_8_complete_log:
            stage_8_complete_log.write("stage 8 (tabulating all replicon copy numbers) finished successfully.\n")

    ## Stage 9: tabulate the length of all chromosomes and plasmids.
    stage_9_complete_file = "../results/stage9.done"
    if exists(stage_9_complete_file):
        print(f"{stage_9_complete_file} exists on disk-- skipping stage 9.")
    else:
        stage9_start_time = time.time()  # Record the start time
        
        tabulate_NCBI_replicon_lengths(reference_genome_dir, replicon_length_csv_file)

        stage9_end_time = time.time()  # Record the end time
        stage9_execution_time = stage9_end_time - stage9_start_time
        Stage9TimeMessage = f"Stage 9 (tabulate all replicon lengths) execution time: {stage9_execution_time} seconds"
        print(Stage9TimeMessage)
        logging.info(Stage9TimeMessage)
        with open(stage_9_complete_file, "w") as stage_9_complete_log:
            stage_9_complete_log.write("stage 9 (tabulating all replicon lengths) finished successfully.\n")

    return


pipeline_main()


