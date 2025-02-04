#!/bin/bash
#SBATCH --job-name=prokaryotes_analysis
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G

mkdir -p {data,results,src}
wget -O data/prokaryotes.txt https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

grep "plasmid" data/prokaryotes.txt | grep "Complete Genome" | sed -i 's/GCA/GCF/g' > results/complete-prokaryotes-with-plasmids.txt

conda create --name PCNdb_env --clone base
conda activate PCNdb_env
pip install pysradb biopython
conda install -c bioconda kallisto breseq
conda install -c conda-forge ncbi-datasets-cli

module load SRA-Toolkit

python PCN_pipeline.py

