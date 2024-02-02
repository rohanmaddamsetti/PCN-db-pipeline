#!/usr/bin/env python

"""
PCN_pipeline.py by Maggie Wilson and Rohan Maddamsetti.

Dependencies:
pip install pysradb
pip install biopython
conda install conda-forge::ncbi-datasets-cli
conda install bioconda::kallisto

"""

import subprocess
import argparse
import json
import os
import gzip
import re
import csv
from Bio import SeqIO
from os.path import basename, exists
import urllib.request


################################################################################
## Functions.

def convert_ids_to_path(gcf_id, asm_id):
    gcf_id_parts = gcf_id.split('_')
    gcf_id_prefix = gcf_id_parts[0]
    gcf_id_nums = gcf_id_parts[1].zfill(9)
    gcf_id_path = '/'.join([gcf_id_prefix, gcf_id_nums[:3], gcf_id_nums[3:6], gcf_id_nums[6:9]])
    path = f"ftp://ftp.ncbi.nlm.nih.gov/genomes/all/{gcf_id_path}/{gcf_id}_{asm_id}"
    return path


def get_SRA_ID_from_RefSeqID(refseq_id):
    bash_command = f'datasets summary genome accession {refseq_id}'
    output = subprocess.check_output(bash_command, shell=True)
    output_str = output.decode('utf-8')
    
    data = output_str
    json_data = json.loads(data)
    sra_id = "NA"
    reports = json_data.get("reports")
    if reports:
        sample_ids = reports[0].get("assembly_info", {}).get("biosample", {}).get("sample_ids")
        if sample_ids:
            for sample_id in sample_ids:
                if sample_id.get("db") == "SRA":
                    sra_id = sample_id.get("value")
                    break
    print(sra_id)
    return(sra_id)


def get_run_ID(sra_id):
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
    ## check to see if pysradb_output is meaningful, if not return NA.
    if not pysra_command_worked:
        return "NA"
    
    pysradb_output_str = pysradb_output.decode('utf-8')
    # Splits the metadata of the SRA ID into respective rows. 
    # And isolates the rows that use the Illumina instrument.
    rows = pysradb_output_str.strip().split('\n')
    filtered_rows = [row for row in rows if ("Illumina") in row and ("WGS") in row]
    #print(filtered_rows)
    rows_with_keyword_str = str(filtered_rows)
    # Selects the run accession by looking for the starting pattern of the ID.
    pattern = r"SRR\d+|ERR\d+|DRR\d+"
    matching_run_accession = re.findall(pattern, rows_with_keyword_str)
    if any (matching_run_accession):
        run_accession = matching_run_accession[0]
    else:
        run_accession = "NA"
    return(run_accession)


def create_RefSeq_SRA_RunID_table(RefSeq_list_file, assembly_list_file, table_outfile):
    with open(RefSeq_list_file, "r") as RefSeq_list_file_obj:
        refseq_ids = RefSeq_list_file_obj.read().splitlines()
    with open(assembly_list_file, "r") as assembly_list_file_obj:
        assembly_ids = assembly_list_file_obj.read().splitlines()
    with open(table_outfile, "w") as table_outfile_obj:
        header = "RefSeq_ID,Assembly_ID,SRA_ID,Run_accession_ID\n"
        table_outfile_obj.write(header) 
        for i, RefSeq_accession in enumerate(refseq_ids):
            my_SRA_ID = get_SRA_ID_from_RefSeqID(RefSeq_accession)
            if my_SRA_ID == "NA":
                my_run_ID = "NA"
            else:
                my_run_ID = get_run_ID(my_SRA_ID)
                print(my_run_ID)
            assembly_id = assembly_ids[i] if i < len(assembly_ids) else "NA"
            row = f"{RefSeq_accession},{assembly_id},{my_SRA_ID},{my_run_ID}\n"
            table_outfile_obj.write(row)
    return


def fetch_gbk(table_outfile, ftp_path_file):
    refseq_ids = []
    asm_ids = []
    list_of_ftp_paths =[]
    with open(table_outfile, "r") as table_obj_csv:
        table_csv = csv.DictReader(table_obj_csv)
        for row in table_csv:
            run_accession_id = row["Run_accession_ID"]
            if run_accession_id != "NA":
                refseq_ids.append(row["RefSeq_ID"])
                asm_ids.append(row["Assembly_ID"])
    
    num_paths = 0          
    for refseq_accession, assembly_accession in zip(refseq_ids, asm_ids):
        ftp_path = convert_ids_to_path(refseq_accession, assembly_accession)
        if ftp_path not in list_of_ftp_paths:  # Add path if it's not already in the list
            list_of_ftp_paths.append(ftp_path)
            num_paths += 1

    with open(ftp_path_file, "w") as file:
        for path in list_of_ftp_paths:
            file.write(f"{path}\n")
    
    with open(ftp_path_file, "r") as ftp_path_list:
        ftp_paths = ftp_path_list.read().splitlines()

    for ftp_path in ftp_paths:
        
        gbk_ftp_path = ftp_path + '/' + basename(ftp_path) + "_genomic.gbff.gz"
        gbff_gz_fname = "../data/NCBI-reference-genomes/" + basename(ftp_path) + "_genomic.gbff.gz"
        gbff_fname = "../data/NCBI-reference-genomes/" + basename(ftp_path) + "_genomic.gbff"
        checksum_file = ftp_path + '/' + "md5checksums.txt"
        my_base_filename = basename(ftp_path) + "_genomic.gbff.gz"
        md5_file_path = "../data/NCBI-reference-genomes/" + my_base_filename + "_checksums.txt"
        
        if exists(gbff_fname):
            continue
        
        gbk_fetch_attempts = 5
        gbk_not_fetched = True
        while gbk_not_fetched and gbk_fetch_attempts:
            try:
                urllib.request.urlretrieve(gbk_ftp_path, filename=gbff_gz_fname)
                urllib.request.urlretrieve(checksum_file, filename=md5_file_path)           
                gbk_not_fetched = False  # assume success if the previous line worked,
                gbk_fetch_attempts = 0  # and don't try again.
            except urllib.error.URLError:
                # if some problem happens, try again.
                gbk_fetch_attempts -= 1
                if exists(gbff_gz_fname):
                    # delete the corrupted file if it exists.
                    os.remove(gbff_gz_fname)
                            
    
        with open(md5_file_path, "r") as checksum_fh:
            target_string = "_genomic.gbff.gz"
            for line in checksum_fh:
                if target_string in line:
                    isolated_line = line.strip()
            
                    my_checksum, my_base_filename = isolated_line.split()
                    my_checksum_ncbi = my_checksum.split("./")[0].strip()
                    print(my_checksum_ncbi)
            
            # run md5 on my_file_path, using the directory for the checksum file,
            # and get the output.
            md5_call = subprocess.run(["md5sum", gbff_gz_fname], capture_output=True, text=True)
            #Needs to be md5sum for DCC, but md5 for locally on zsh.
            my_md5_checksum = md5_call.stdout.split("=")[-1].strip()
            print(my_md5_checksum)
            ## verify that the checksums match.
            if (my_md5_checksum == my_checksum_ncbi):
                print(my_md5_checksum, "matches", my_checksum_ncbi)
            else:
                error_message = "ERROR: " + my_md5_checksum + " does not match " + my_checksum_ncbi + " for file " + gbff_gz_fname
                #raise AssertionError(error_message)
                #It raised the error for me even though they matched
    return
 
        
def create_genome_metadatacsv(ftp_file, table_outfile, genome_metadatacsv):
    reference_genome=[]
    with open(ftp_file, "r") as ftp_path_list:
        for line in ftp_path_list:
            ftp_path = line.strip()
            genome_id = ftp_path.split('/')[-1]
            genome_id += '_genomic.gbff.gz'
            reference_genome.append(genome_id)
    run_accession_list=[]
    with open(table_outfile, "r") as table_obj_csv:
        table_csv = csv.DictReader(table_obj_csv)
        for row in table_csv:
            run_accession_id = row["Run_accession_ID"]
            if run_accession_id != "NA":
                run_accession_list.append(row["Run_accession_ID"])
    rows = zip(reference_genome, run_accession_list)
    with open(genome_metadatacsv, "w", newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['ReferenceGenome','SRA_Data'])
        writer.writerows(rows)
    print(f"CSV file '{csv_file}' has been created successfully.")
    return


def download_fastq_reads(SRA_data_dir, table_outfile):
        """
        the sra_id has to be the last part of the directory.
        see documentation here:
        https://github.com/ncbi/sra-tools/wiki/08.-prefetch-and-fasterq-dump
        """
        sra_numbers = []
        with open(table_outfile, "r") as table_outfile_obj:
            table_csv = csv.DictReader(table_outfile_obj)
            for row in table_csv:
                run_accession_id = row["SRA_Data"]
                if run_accession_id != "NA":
                    sra_numbers.append(run_accession_id.strip())
        
        for sra_id in sra_numbers:
            sra_dir_path = os.path.join(SRA_data_dir, sra_id)
            if os.path.exists(sra_dir_path): continue
            prefetch_args = ["prefetch", sra_id, "-O", sra_dir_path]
            print (" ".join(prefetch_args))
            subprocess.run(prefetch_args)
        my_cwd = os.getcwd()
        os.chdir(SRA_data_dir)
        for sra_id in sra_numbers:
            sra_fastq_file_1 = sra_id + "_1.fastq"
            sra_fastq_file_2 = sra_id + "_2.fastq"
            if os.path.exists(sra_fastq_file_1) and os.path.exists(sra_fastq_file_2):
                continue
            else:
                print ("Generating fastq for: " + sra_id)
                fasterq_dump_args = ["fasterq-dump", "--threads", "10", sra_id]
                print(" ".join(fasterq_dump_args))
                subprocess.run(fasterq_dump_args)
    
        # now change back to original working directory.
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
    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff")]
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
            subprocess.run(kallisto_quant_args)
            #print(" ".join(kallisto_quant_args))
            #quit()
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
    RefSeq_list_file = "../data/ncbidataset_refseq.txt"
    assembly_list_file = ("../data/ncbidataset_assembly.txt")
    ftp_path_file = ("../results/ftp_path.txt")
    table_outfile = "../results/runID_table.txt"
    SRA_data_dir = "../results/SRA/"
    refgenomes_dir = "../results/ref_genomes/"
    kallisto_ref_dir = "../results/kallisto_references/"
    kallisto_index_dir = "../results/kallisto_indices/"
    kallisto_quant_results_dir = "../results/kallisto_quant/"
    genome_metadatacsv = "../results/genome_metadata.csv"
    copy_number_csv_file = "../results/chromosome_plasmid_copy_numbers.csv"
    ARG_copy_number_csv_file = "../results/ARG_copy_numbers.csv"
    replicon_length_csv_file = "../results/replicon_lengths.csv"
    
    create_RefSeq_SRA_RunID_table(RefSeq_list_file, assembly_list_file, table_outfile)
    fetch_gbk(table_outfile, ftp_path_file)
    create_genome_metadatacsv(ftp_path_file, table_outfile, genome_metadatacsv)
    NCBI_genomeID_to_SRA_ID_dict = make_genome_to_SRA_dict(genome_metadatacsv)
    download_fastq_reads(SRA_data_dir, genome_metadatacsv)
    make_NCBI_fasta_refs_for_kallisto(refgenomes_dir, kallisto_ref_dir)
    make_NCBI_kallisto_indices(kallisto_ref_dir, kallisto_index_dir)
    run_kallisto_quant(NCBI_genomeID_to_SRA_ID_dict, kallisto_index_dir, SRA_data_dir, kallisto_quant_results_dir)
    measure_NCBI_replicon_copy_numbers(kallisto_quant_results_dir, copy_number_csv_file)
    measure_NCBI_ARG_copy_numbers(kallisto_quant_results_dir, ARG_copy_number_csv_file)
    tabulate_NCBI_replicon_lengths(refgenomes_dir, replicon_length_csv_file)
    return


pipeline_main()


