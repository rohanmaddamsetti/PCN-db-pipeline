#!/usr/bin/env python

"""
PCN_pipeline.py by Maggie Wilson and Rohan Maddamsetti.

For this pipeline to work, ncbi datasets, pysradb, and kallisto must be in the $PATH.

Currently, this pipeline only analyzes Illumina short-read data with kallisto.

In the future, one could additionally analyze long-read data as well, by using the
state-of-the-art in pseudoalignment analyses, Themisto (published 2023 in Bioinformatics).

"""

import subprocess
import argparse
import json
import os
import sys
import gzip
import re
import csv
from Bio import SeqIO
from os.path import basename, exists
import urllib.request
import time
import logging


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
    with open(outfile, "w") as outfh:
        with open(gbk_gz_path,'rt') as gbk_fh:
            SeqID = None
            SeqType = None
            for i, record in enumerate(SeqIO.parse(gbk_fh, "genbank")):
                SeqID = record.id
                if "complete" in record.description:
                    if "plasmid" in record.description:
                        SeqType = "plasmid"
                    elif "chromosome" in record.description or i == 0:
                        ## IMPORTANT: we assume here that the first record is a chromosome.
                        SeqType = "chromosome"
                    else:
                        continue
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
    print(gzfilelist)
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        genome_id = gzfile.split(".gbff")[0]
        fasta_outfile = os.path.join(kallisto_ref_outdir, genome_id+".fna")
        print("making: ", fasta_outfile)
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


def run_kallisto_quant(NCBI_genomeID_to_SRA_ID_dict, kallisto_index_dir, SRA_data_dir, results_dir):
    index_list = [x for x in os.listdir(kallisto_index_dir) if x.endswith(".idx")]
    for index_file in index_list:
        index_path = os.path.join(kallisto_index_dir, index_file)
        genome_id = index_file.split(".idx")[0]
        SRA_id = NCBI_genomeID_to_SRA_ID_dict[genome_id]
        if SRA_id == "NA":
            continue
        else:
            read_file1 = SRA_id + "_1.fastq"
            read_file2 = SRA_id + "_2.fastq"
            read_path1 = os.path.join(SRA_data_dir, read_file1)
            read_path2 = os.path.join(SRA_data_dir, read_file2)
            output_path = os.path.join(results_dir, genome_id)
            ## run with 10 threads by default.
            kallisto_quant_args = ["kallisto", "quant", "-t", "8", "-i", index_path, "-o", output_path, "-b", "100", read_path1, read_path2]
            ##print(" ".join(kallisto_quant_args))
            subprocess.run(kallisto_quant_args)
    return


def make_genome_to_SRA_dict(NCBI_metadata_csv):
    genome_to_SRA_dict = dict()
    with open(NCBI_metadata_csv, "r") as csv_fh:
        for i, line in enumerate(csv_fh):
            if i == 0: continue ## skip the header.
            line = line.strip() 
            ReferenceGenome, SRA_Data = line.split(',')
            genome_id = ReferenceGenome.split(".gbff.gz")[0]
            genome_to_SRA_dict[genome_id] = SRA_Data
    return genome_to_SRA_dict


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
        print(value_tuple)
        print(seqtype)
        print(coverage)
        print(chromosome_coverage)
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


def isARG(product_annotation):
    chloramphenicol_keywords = "chloramphenicol|Chloramphenicol"
    tetracycline_keywords = "tetracycline efflux|Tetracycline efflux|TetA|Tet(A)|tetA|tetracycline-inactivating"
    MLS_keywords = "macrolide|lincosamide|streptogramin"
    multidrug_keywords = "Multidrug resistance|multidrug resistance|antibiotic resistance"
    beta_lactam_keywords = "lactamase|LACTAMASE|beta-lactam|oxacillinase|carbenicillinase|betalactam\S*"
    glycopeptide_keywords = "glycopeptide resistance|VanZ|vancomycin resistance|VanA|VanY|VanX|VanH|streptothricin N-acetyltransferase"
    polypeptide_keywords = "bacitracin|polymyxin B|phosphoethanolamine transferase|phosphoethanolamine--lipid A transferase"
    diaminopyrimidine_keywords = "trimethoprim|dihydrofolate reductase|dihydropteroate synthase"
    sulfonamide_keywords = "sulfonamide|Sul1|sul1|sulphonamide"
    quinolone_keywords = "quinolone|Quinolone|oxacin|qnr|Qnr"
    aminoglycoside_keywords = "Aminoglycoside|aminoglycoside|streptomycin|Streptomycin|kanamycin|Kanamycin|tobramycin|Tobramycin|gentamicin|Gentamicin|neomycin|Neomycin|16S rRNA (guanine(1405)-N(7))-methyltransferase|23S rRNA (adenine(2058)-N(6))-methyltransferase|spectinomycin 9-O-adenylyltransferase|Spectinomycin 9-O-adenylyltransferase|Rmt"
    macrolide_keywords = "macrolide|ketolide|Azithromycin|azithromycin|Clarithromycin|clarithromycin|Erythromycin|erythromycin|Erm|EmtA"
    antimicrobial_keywords = "QacE|Quaternary ammonium|quaternary ammonium|Quarternary ammonium|quartenary ammonium|fosfomycin|ribosomal protection|rifampin ADP-ribosyl|azole resistance|antimicrob\S*"
    ARG_regex = "|".join([chloramphenicol_keywords, tetracycline_keywords,
                          MLS_keywords, multidrug_keywords, beta_lactam_keywords,
                          glycopeptide_keywords, polypeptide_keywords, diaminopyrimidine_keywords,
                          sulfonamide_keywords, quinolone_keywords, aminoglycoside_keywords,
                          macrolide_keywords, antimicrobial_keywords])
    if re.search(ARG_regex, product_annotation): return True
    return False


def estimate_ARG_copy_numbers(genecount_tsv_path):

    chromosomal_gene_length = 0.0
    chromosomal_gene_est_counts = 0.0

    ARG_coverage_dict = dict()
    ## get the chromosomal gene coverage, and get the coverage for all ARGs.
    with open(genecount_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, locus_tag, product = parse_metadata_in_header(target_id)
            if SeqType == "chromosome":
                chromosomal_gene_length += float(length)
                chromosomal_gene_est_counts += float(est_counts)
            if isARG(product):
                coverage = float(est_counts) / float(length)
                ARG_coverage_dict[locus_tag] = (SeqID, SeqType, product, coverage)
    ## NOTE: GCF_026154285.1_ASM2615428v1 did not have any reads pseudoalign.
    ## Return an empty dict() when nothing aligns to the chromosome.
    if chromosomal_gene_length == 0:
        print("WARNING: no reads pseudoaligned in file: ", genecount_tsv_path)
        print("estimate_ARG_copy_numbers is returning an empty dict.")
        return(dict())
    chromosome_coverage = chromosomal_gene_est_counts/chromosomal_gene_length
    ## now normalize by chromosome coverage to get copy number estimates.
    ARG_copy_number_dict = dict()
    for locus_tag, value_tuple in ARG_coverage_dict.items():
        my_SeqID, my_SeqType, my_product, my_coverage = value_tuple
        ARG_copy_number_dict[locus_tag] = (my_SeqID, my_SeqType, my_product, my_coverage/chromosome_coverage)
    return(ARG_copy_number_dict)


def measure_NCBI_ARG_copy_numbers(kallisto_quant_results_dir, ARG_copy_number_csv_file):
    """
    define lists to encode the following columns of the table.
    AnnotationAccession, SeqID, SeqType, locus_tag, product, CopyNumber
    """
    AnnotationAccessionVec = []
    SeqIDVec = [] ## this is for the replicon.
    SeqTypeVec = []
    LocusTagVec = []
    ProductVec = []
    CopyNumberVec = []
    
    ## skip .DS_Store and any other weird files.
    genomedirectories = [x for x in os.listdir(kallisto_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        ## I probably should have trimmed the '_genomic' suffix in an earlier step.
        annotation_accession = genomedir.split("_genomic")[0]
        genome_quantfile_path = os.path.join(kallisto_quant_results_dir, genomedir, "abundance.tsv")
        ARG_copy_number_dict = estimate_ARG_copy_numbers(genome_quantfile_path)
        for locus_tag, value_tuple in ARG_copy_number_dict.items():
            SeqID, seqtype, product, copy_number = value_tuple
            AnnotationAccessionVec.append(annotation_accession)
            SeqIDVec.append(SeqID)
            SeqTypeVec.append(seqtype)
            LocusTagVec.append(locus_tag)
            ProductVec.append(product)
            CopyNumberVec.append(copy_number)

    assert len(AnnotationAccessionVec) == len(SeqIDVec) == len(SeqTypeVec) == len(LocusTagVec) == len(ProductVec) == len(CopyNumberVec)
    ## now write the ARG copy number data to file.
    with open(ARG_copy_number_csv_file, "w") as outfh:
        header = "AnnotationAccession,SeqID,SeqType,locus_tag,product,CopyNumber"
        outfh.write(header + "\n")
        for i in range(len(AnnotationAccessionVec)):
            outfh.write(AnnotationAccessionVec[i] + "," + SeqIDVec[i] + "," + SeqTypeVec[i] + "," + LocusTagVec[i] + "," + ProductVec[i] + "," + str(CopyNumberVec[i]) + "\n")
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
    RunID_table_file = "../results/RunID_table.csv"
    reference_genome_dir = "../data/NCBI-reference-genomes/"
    SRA_data_dir = "../data/SRA/"
    
    kallisto_ref_dir = "../results/kallisto_references/"
    kallisto_index_dir = "../results/kallisto_indices/"
    kallisto_quant_results_dir = "../results/kallisto_quant/"
    genome_metadatacsv = "../results/genome_metadata.csv"
    copy_number_csv_file = "../results/chromosome_plasmid_copy_numbers.csv"
    ARG_copy_number_csv_file = "../results/ARG_copy_numbers.csv"
    replicon_length_csv_file = "../results/replicon_lengths.csv"

    ## Stage 1: get SRA IDs and Run IDs for all complete RefSeq bacterial genomes with plasmids.
    if exists(RunID_table_file):
        Stage1DoneMessage = f"{RunID_table_file} exists on disk-- skipping stage 1."
        print(Stage1DoneMessage)
        logging.info(Stage1DoneMessage)
    else: ## This takes 34513 seconds (9.5h) to get RunIDs for 4921 genomes.
        RunID_table_start_time = time.time()  # Record the start time
        create_RefSeq_SRA_RunID_table(prokaryotes_with_plasmids_file, RunID_table_file)
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
        fetch_reference_genomes(RunID_table_file, refseq_accession_to_ftp_path_dict, reference_genome_dir)
        with open(stage_2_complete_file, "w") as stage_2_complete_log:
            stage_2_complete_log.write("reference genomes downloaded successfully.\n")
        
    ## Stage 3: download Illumina reads for the genomes from the NCBI Short Read Archive (SRA).
    stage_3_complete_file = "../results/stage3.done"
    if exists(stage_3_complete_file):
        print(f"{stage_3_complete_file} exists on disk-- skipping stage 3.")
    else:
        SRA_download_start_time = time.time()  # Record the start time
        download_fastq_reads(SRA_data_dir, RunID_table_file)
        SRA_download_end_time = time.time()  # Record the end time
        SRA_download_execution_time = SRA_download_end_time - SRA_download_start_time
        Stage3TimeMessage = f"Stage 3 (SRA download) execution time: {SRA_download_execution_time} seconds"
        print(Stage3TimeMessage)
        logging.info(Stage3TimeMessage)
        with open(stage_3_complete_file, "w") as stage_3_complete_log:
            stage_3_complete_log.write("SRA read data downloaded successfully.\n")


    ## Stage 4:
            
    ##make_NCBI_fasta_refs_for_kallisto(reference_genome_dir, kallisto_ref_dir)
    ##make_NCBI_kallisto_indices(kallisto_ref_dir, kallisto_index_dir)
    ##run_kallisto_quant(NCBI_genomeID_to_SRA_ID_dict, kallisto_index_dir, SRA_data_dir, kallisto_quant_results_dir)
    ##measure_NCBI_replicon_copy_numbers(kallisto_quant_results_dir, copy_number_csv_file)
    ##measure_NCBI_ARG_copy_numbers(kallisto_quant_results_dir, ARG_copy_number_csv_file)
    ##tabulate_NCBI_replicon_lengths(reference_genome_dir, replicon_length_csv_file)

    ## TODO: make a table of the estimated copy number and position for all genes in all chromosomes
    ## and plasmids in these genomes. My reasoning is that this may be useful for doing some analyses
    ## like in David Zeevi's science paper about growth rates from chromosomal copy numbers.

    return


pipeline_main()


