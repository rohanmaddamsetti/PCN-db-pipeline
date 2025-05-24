# PCN-db-pipeline by Rohan Maddamsetti, Maggie Wilson, and Irida Shyti.

A pipeline for analyzing plasmid copy numbers (PCN) in bacterial genomes.


## Table of Contents
1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Setup](#setup)
4. [Running the Pipeline](#running-the-pipeline)
5. [Expected Output](#expected-output)
6. [Performance](#performance)


## Overview

This pipeline analyzes plasmid copy numbers (PCN) in bacterial genomes using:
- NCBI RefSeq genome annotation data
- NCBI SRA Illumina short-read sequencing data
- Themisto for pseudoalignment
- Kallisto for transcript quantification (not critical for PCN estimation, used for control experiments)
- Breseq for mutation analysis (not critical for PCN estimation, used for control experiments)

All classes and functions are in the source code file src/PCN_library.py. Each stage of the pipeline is described in the comments in the source code file src/PCN_pipeline.py.  


## Requirements

### Software
- Python 3.11+
- Biopython 1.85
- SRA-Toolkit 3.2.0
- pysradb 2.2.2
- ncbi-datasets-cli 17.1.0
- Themisto 3.2.2
- Breseq 0.37+
- Kallisto 0.51.1

### Hardware
- Linux High Performance Computing Cluster (HPC) with SLURM job scheduler (https://slurm.schedmd.com/)
- Recommended: Duke Compute Cluster (DCC)
- Storage: ~15TB for raw sequencing data
- Memory: 16GB minimum
- Time: ~2 week for full pipeline: 3 days for Illumina short-read data download, and a few days for PCN estimation on HPC.

This pipeline can be run locally or on Duke Compute Cluster (DCC).  
DCC or your high-performance computing cluster (HPC) is recommended,  
due to the large amount of FASTQ sequencing data downloaded for PCN estimation.  


## Setup

### Notes before setting up the pipeline.

This pipeline should be run on a Linux HPC system with a SLURM job manager.
The pipeline *can* be run on MacOS (we wrote and debugged most of the pipeline on MacOS) using the command “python PCN_pipeline.py”; however, we recommend that users *only* run the pipeline on MacOS for debugging or testing individual stages for the following reasons:
1) Several stages of the pipeline submit thousands of HPC jobs in parallel to speed up computation, and this is not possible on a laptop.
2) The SRA data download is substantial– ~15TB of sequencing reads – and so a full download is not possible.  

1. Create project structure by downloading or cloning this github repository.
Alternatively, create a new top-level project directory
with containing src/, data/, results/ directories.
Then, copy the source code in this github repository into the src/ directory for your project.

   ```bash
   mkdir my-test-PCN-project
   cd my-test-PCN-project
   mkdir -p {data,results,src}
   ```

2. Download required data:
   ```bash
   wget -O data/prokaryotes.txt https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt ## on linux
   curl https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt > data/prokaryotes.txt ## on mac and linux
   ```

3. Filter for complete genomes containing plasmids, and change GenBank IDs (GCA_*) to RefSeq Accessions (GCF_*).
   ```bash
   grep "plasmid" data/prokaryotes.txt | grep "Complete Genome" | sed 's/GCA/GCF/g' > results/complete-prokaryotes-with-plasmids.txt
   ```

4. Set up conda environment and install python dependencies:
   ```bash
   conda create --name PCNdb_env --clone base
   conda activate PCNdb_env
   pip install pysradb biopython HTSeq beautifulsoup4 polars
   conda install -c bioconda kallisto breseq
   conda install -c conda-forge ncbi-datasets-cli
   ```

5. Install SRA-Toolkit:
   - On DCC: `module load SRA-Toolkit`
   - Locally: Download from [SRA-Tools GitHub](https://github.com/ncbi/sra-tools)


## Running the Pipeline

   ```bash
   conda activate PCNdb_env  
   cd src/  
   sbatch --mem=16G -t 430:00:00 -p youlab --wrap="python PCN_pipeline.py"
   ```

Note that the pipeline quits at the end of each stage (progress is saved). Therefore, one has to run the pipeline anew to start the next stage.


## Expected Output

The pipeline generates:
1. Reference genomes in `data/NCBI-reference-genomes/`
2. SRA data in `data/SRA/`
3. Analysis results:
   - `results/PIRA-PCN-estimates.csv`
   - `results/NCBI-replicon_lengths.csv`
   - `results/kallisto-replicon_copy_numbers.csv`
   - `results/themisto-replicon-read-counts.csv`  


## Performance

| Stage | Time Estimate | Data Size |
|-------|---------------|-----------|
| Reference Genome Download | ~1-2 days | ~50GB |
| SRA Data Download | ~3 days | ~15TB |
| Full Pipeline | ~2 weeks | ~15TB |

Note: Run the full pipeline in the `/work` directory on DCC since 15TB+ of disk space is needed.  

