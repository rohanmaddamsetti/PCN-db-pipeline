# PCN-db-pipeline by Rohan Maddamsetti, Maggie Wilson, and Irida Shyti.
## Requirements: biopython, kallisto, SRA-Toolkit, pysradb, and ncbi-datasets-cli

### This program can either be run locally or on Duke Compute Cluster (DCC). DCC is suggested to use due to the large amount of data that is downloaded from the thousands of samples. 

### Make a top-level directory with three directories inside, named "data", "results", and "src". Now copy all source code files in this repository into "src".

#### The following file needs to be downloaded into the data/ directory.
#### data/prokaryotes.txt should be downloaded from: https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

#### results/prokaryotes-with-chromosomes-and-plasmids.txt should then be generated with the following command (run from src/ directory):
(head -n 1 ../data/prokaryotes.txt && grep "plasmid" ../data/prokaryotes.txt | grep "chromosome") > ../results/prokaryotes-with-chromosomes-and-plasmids.txt  

This command ensures that every genome has both an annotated chromosome and at least one annotated plasmid.

SRA data for ~6000 genomes was downloaded, but only have about ~4500 reference genomes. to understand this discrepancy, I ran:
cd ../data/NCBI-reference-genomes
ls | grep ".gbff.gz" | sed 's/_genomic.gbff.gz$//' > ../../results/downloaded-genome-ids.txt


### Expected Output:
#### The first steps of this pipeline, convert these sample IDs into run accession IDs and create the reference genome path. Two txt files, a table with run accession IDs and a list of the genome paths  will be created. Then, the sample's reference genome will be downloaded if a run accession ID is present, these genomes will be downloaded into the ref-genomes folder in results. There should be 4,540 genomes downloaded, as that is the number of samples that have a run ID.
#### A metadata csv is then created, a master list of all the samples were are downloading data for. Next, fastq files are downloaded for each run ID and put into the SRA folder in results.
#### The next several functions parse through the SRA and reference genomes using kallisto and themisto.
#### Lastly the final tables of results are created as "PIRA-PCN-estimates.csv", and "NCBI-replicon_lengths.csv".

## How to run on DCC or locally.
#### create a new conda environment (here I clone the base environment):
##### conda create --name PCNdb_env --clone base
##### conda activate PCNdb_env
#### install pysradb and biopython in this environment using pip:
##### pip install pysradb
##### pip install biopython
#### install kallisto and ncbi-datasets-cli in this new environment
##### conda install bioconda::kallisto
##### conda install conda-forge::ncbi-datasets-cli
##### conda install bioconda::breseq


#### You will then need to load the SRA-Toolkit module that is available on DCC:
##### module load SRA-Toolkit
#### alternatively, install a pre-built binary of SRA-Toolkit on your machine from here: https://github.com/ncbi/sra-tools

### This will take ~2 weeks to run on DCC.

### The Illumina data download (stage 3) takes 281 hours (12 days) to download 15 TB of raw sequencing data from NCBI Short Read Archive (SRA).  

### sbatch --mem=16G -t 430:00:00 -p youlab --wrap="python PCN_pipeline.py"
#### Use .../work, as ~15 terabytes will be downloaded.
