#Use environment breseqwork or an environment where you have downloaded breseq

import subprocess
import os
import gzip
import re
from Bio import SeqIO
from bs4 import BeautifulSoup



def run_breseq_on_DCC(metadata_csv, unzipped_refgenomes_dir, SRA_data_dir, breseq_results_dir):
    with open(metadata_csv, "r") as metadata_fh:
        for i, line in enumerate(metadata_fh):
            line = line.strip()
            if i == 0: continue ## skip the header.
            gbk_gz, sra_id = line.split(',')
            unzipped_gbk = gbk_gz.split(".gz")[0]
            annotation_accession = gbk_gz.split("_genomic.gbff.gz")[0]

            breseq_outdir = os.path.join(breseq_results_dir, annotation_accession)
            cur_reference_gbk = os.path.join(unzipped_refgenomes_dir, unzipped_gbk)
            fastq1 = os.path.join(SRA_data_dir, sra_id+"_1.fastq")
            fastq2 = os.path.join(SRA_data_dir, sra_id+"_2.fastq")
    
            breseq_args = ["breseq", "-o", breseq_outdir, "-r", cur_reference_gbk, fastq1, fastq2]
            breseq_string = " ".join(breseq_args)
            #sbatch_string = "sbatch -p scavenger -t 2:30:00 --mem=8G --wrap=\"" + breseq_string + "\""
            #print(sbatch_string)
            subprocess.run(breseq_string, shell=True)
    return

def run_random_breseq_on_DCC(metadata_csv, unzipped_refgenomes_dir, SRA_data_dir, breseq_results_dir):
    with open(metadata_csv, "r") as metadata_fh:
        metadata_list = []
        for i,line in enumerate(metadata_fh):
            line = line.strip()
            if i == 0: continue ## skip the header.
            metadata_list.append(line)
        #Here I will randomize the list.
        randomized_metadata_list = metadata_list #make random
        truncated_randomized_metadata_list = randomized_metadata_list[100:]
            
        for line in truncated_randomized_metadata_list:
            gbk_gz, sra_id = line.split(',')
            unzipped_gbk = gbk_gz.split(".gz")[0]
            annotation_accession = gbk_gz.split("_genomic.gbff.gz")[0]

            breseq_outdir = os.path.join(breseq_results_dir, annotation_accession)
            cur_reference_gbk = os.path.join(unzipped_refgenomes_dir, unzipped_gbk)
            fastq1 = os.path.join(SRA_data_dir, sra_id+"_1.fastq")
            fastq2 = os.path.join(SRA_data_dir, sra_id+"_2.fastq")
    
            breseq_args = ["breseq", "-o", breseq_outdir, "-r", cur_reference_gbk, fastq1, fastq2]
            breseq_string = " ".join(breseq_args)
            #sbatch_string = "sbatch -p scavenger -t 2:30:00 --mem=8G --wrap=\"" + breseq_string + "\""
            #print(sbatch_string)
            subprocess.run(breseq_string, shell=True)
    return


def parse_breseq_results(breseq_outdir, results_csv_path):
    with open(results_csv_path, "w") as csv_fh:
        ## print a header string.
        print("AnnotationAccession,SeqID,replicon_length,mean_coverage,description,SeqType,CopyNumber", file=csv_fh)
        ## filter on IDs starting with "GCF_"
        for genome_accession in [x for x in os.listdir(breseq_outdir) if x.startswith("GCF_")]:
            breseq_summary_path = os.path.join(breseq_outdir, genome_accession, "output", "summary.html")
            ## Read the HTML file if it exists.
            if not os.path.exists(breseq_summary_path): continue
            with open(breseq_summary_path, 'r') as summary_fh:
                html_content = summary_fh.read()

            ## Create a BeautifulSoup object
            soup = BeautifulSoup(html_content, 'html.parser')

            ## Find the table by its section heading
            table_section = soup.find('h2', string='Reference Sequence Information')
            reference_table = table_section.find_next('table')

            ## Extract the table data
            table_data = []
            for row in reference_table.find_all('tr'):
                row_data = [cell.get_text(strip=True) for cell in row.find_all('td')]
                if row_data and row_data[0] == "coverage":
                    table_data.append(row_data)

            ## Print the extracted table data
            for i, row in enumerate(table_data):
                SeqID = row[2]
                replicon_length = row[3].replace(",", "") ## remove commas in the length.
                mean_coverage = row[4]
                ## for CSV compatibility, replace commas with semicolons.
                description = row[-1].replace(",",";")
                if i == 0:
                    SeqType = "chromosome"
                    chromosome_coverage = float(mean_coverage)
                else:
                    SeqType = "plasmid"
                ## handle cases where coverage fit failed.
                if mean_coverage == "NA" or mean_coverage.startswith('['):
                    copy_number = "NA"
                else:
                    copy_number = str(float(mean_coverage)/chromosome_coverage)
                csv_string = ",".join([genome_accession, SeqID, replicon_length, mean_coverage, description, SeqType, copy_number])
                print(csv_string, file=csv_fh)
    return


def main():
    SRA_data_dir = "../results/PCN_pipeline/SRA/"
    metadata_csv = "../results/PCN_pipeline/genome_metadata.csv"
    unzipped_refgenomes_dir = "../results/PCN_pipeline/ref_genomes/"
    breseq_outdir = "../results/PCN_pipeline/breseq/"
    breseq_results_csv = "../results/PCN_pipeline/breseq_copy_num.csv"
    run_breseq_on_DCC(metadata_csv, unzipped_refgenomes_dir, SRA_data_dir, breseq_outdir)
    parse_breseq_results(breseq_outdir, breseq_results_csv)
    return

main()
        
