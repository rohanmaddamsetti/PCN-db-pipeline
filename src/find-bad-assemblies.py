#!/usr/bin/env python

"""
find-bad-assemblies.py by Rohan Maddamsetti.

This helper script finds the names of bad reference genomes (those without chromosomes)
and their associated SRA data files,

and deletes the SRA data and reference genome files and checksum files corresponding to these bad quality
genomes (contig assemblies without chromosomes).

"""

import os
import subprocess


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

def find_bad_ones():
    good_ones_txt = "../results/prokaryotes-with-chromosomes-and-plasmids.txt"
    good_ones_dict = create_refseq_accession_to_ftp_path_dict(good_ones_txt)
    with open("../results/BAD-RunID_table.csv", "w") as bad_runids_fh, open("../results/OLD-RunID_table.csv", "r") as runids_fh:
        for i, line in enumerate(runids_fh):
            if i == 0:
                bad_runids_fh.write(line)
                print(line)
            else:
                fields = line.split(",")
                cur_ref_seq_id = fields[0]
                if cur_ref_seq_id not in good_ones_dict.keys():
                    bad_runids_fh.write(line)
                    print(line)
    return


def find_good_ones():
    good_ones_txt = "../results/prokaryotes-with-chromosomes-and-plasmids.txt"
    good_ones_dict = create_refseq_accession_to_ftp_path_dict(good_ones_txt)
    with open("../results/GOOD-RunID_table.csv", "w") as good_runids_fh, open("../results/OLD-RunID_table.csv", "r") as runids_fh:
        for i, line in enumerate(runids_fh):
            if i == 0:
                good_runids_fh.write(line)
                print(line)
            else:
                fields = line.split(",")
                cur_ref_seq_id = fields[0]
                if cur_ref_seq_id in good_ones_dict.keys():
                    good_runids_fh.write(line)
                    print(line)
    return


def delete_bad_SRA_files():
    with open("../results/BAD-RunID_table.csv", "r") as badones_fh:
        for i, line in enumerate(badones_fh):
            if i == 0: continue ## skip header
            bad_SRA_ID = line.strip().split(",")[-1]
            bad_SRA_ID_argument = "../data/SRA/"+ bad_SRA_ID + "*"
            rm_args = ["rm ", "-rf", bad_SRA_ID_argument]
            rm_arg_string = " ".join(rm_args)
            print(rm_arg_string)
            subprocess.run(rm_arg_string, shell=True)
    return


def delete_bad_NCBI_reference_genome_files():
    refgenomes_dir = "../data/NCBI-reference-genomes/"

    bad_genome_IDs = []
    with open("../results/BAD-RunID_table.csv", "r") as badones_fh:
        for i, line in enumerate(badones_fh):
            if i == 0: continue ## skip header
            bad_genome_ID = line.strip().split(",")[0]
            bad_genome_IDs.append(bad_genome_ID)
    
    GCF_filelist = [x for x in os.listdir(refgenomes_dir) if x.startswith("GCF_")]
    for GCF_file in GCF_filelist:
        my_GCF_filepath = os.path.join(refgenomes_dir, GCF_file)
        for bad_genome_ID in bad_genome_IDs:
            if GCF_file.startswith(bad_genome_ID):
                rm_args = ["rm ", my_GCF_filepath]
                rm_arg_string = " ".join(rm_args)
                print(rm_arg_string)
                subprocess.run(rm_arg_string, shell=True)
                break
    return


def main():

    ##find_good_ones()
    ## find_bad_ones()
    ##delete_bad_SRA_files()
    ##delete_bad_NCBI_reference_genome_files()

    ## understand the discrepancy between the SRA data downloaded for ~6000 genomes,
    ## but only ~4500 reference genomes downloaded.
    downloaded_reference_genome_IDs = []
    with open("../results/downloaded-genome-ids.txt") as downloaded_ref_genomes_fh:
        for line in downloaded_ref_genomes_fh:
            my_id = line.strip()
            downloaded_reference_genome_IDs.append(my_id)

    downloaded_SRA_data_genome_IDs = []
    with open("../results/RunID_table.csv", "r") as goodones_fh:
        for i, line in enumerate(goodones_fh):
            if i == 0: continue ## skip header
            genome_ID = line.strip().split(",")[0]
            downloaded_SRA_data_genome_IDs.append(genome_ID)

    missing_reference_genome_ids = list(set(downloaded_SRA_data_genome_IDs) - set(downloaded_reference_genome_IDs))
    
    with open("../results/missing-reference-genomes.txt", "w") as missing_genomes_fh:
        for x in missing_reference_genome_ids:
            missing_genomes_fh.write(x + "\n")
    

main()
