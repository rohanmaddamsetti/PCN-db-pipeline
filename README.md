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
- minimap2 for re-aligning multireads
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
- minimap2 2.29 
- Breseq 0.37+
- Kallisto 0.51.1

### Hardware
- Required: Linux High Performance Computing Cluster (HPC) with SLURM job scheduler (https://slurm.schedmd.com/)
- Storage: ~15TB for raw sequencing data
- Memory: 16GB minimum
- Time: ~2 week for full pipeline: 3 days for Illumina short-read data download, and a few days for PCN estimation on HPC.
- Recommended: Duke Compute Cluster (DCC)

## Notes before setting up the pipeline

The pipeline is currently set to Test Mode, so that a limited number of genomes are downloaded.
In Test Mode, the pipeline *can* be run on MacOS as follows:

```
python PCN_pipeline.py
```

We recommend that users *only* run the pipeline on MacOS for debugging or testing individual stages for the following reasons:

1) Several stages of the pipeline submit thousands of HPC jobs in parallel to speed up computation, and this is not possible on a laptop.
2) The SRA data download is substantial– ~15TB of sequencing reads – and so a full download is not possible.  


To run the full pipeline, open PCN_pipeline.py in your favorite text editor and set `TEST_MODE = False` in line 25.
**The full pipeline should be run on a Linux HPC system with a SLURM job manager.**

## Setup

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

3. Filter for complete genomes containing plasmids, and change GenBank IDs (GCA_\*) to RefSeq Accessions (GCF_\*).
   ```bash
   grep "plasmid" data/prokaryotes.txt | grep "Complete Genome" | sed 's/GCA/GCF/g' > results/complete-prokaryotes-with-plasmids.txt
   ```

4. Set up conda environment and install python dependencies:
   ```bash
   conda create --name PCNdb_env --clone base
   conda activate PCNdb_env
   pip install pysradb biopython HTSeq beautifulsoup4 polars
   conda install -c bioconda kallisto breseq
   conda install -c conda-forge ncbi-datasets-cli tqdm
   ```

5. Install SRA-Toolkit:
   - On DCC: `module load SRA-Toolkit/3.0.0-rhel8` (users on other HPC systems will probably have a different version number).
   - Locally: Download from [SRA-Tools GitHub](https://github.com/ncbi/sra-tools)

6. Install remaining required dependencies and external tools (see requirements above).

- **External Tools:**  
  - [`themisto`](https://github.com/algbio/themisto/releases/) – must be installed and available in your `$PATH`
  - [`minimap2`](https://github.com/lh3/minimap2/releases) – must be installed and available in your `$PATH`  

    For convenience, MacOS and Linux binaries for `themisto` and a Linux binary for `minimap2` are provided in the `bin/` directory. You may have to turn these binaries into executables like so:

    ```bash
    chmod +x ${PWD}/bin/themisto_linux-v3.2.2/themisto  ##make the linux themisto binary into an executable
    ```

    Then, you can add these binaries to the `$PATH` from the command-line as follows, before running the pipeline:

    ```bash
    export PATH="${PWD}/bin/themisto_linux-v3.2.2:$PATH"
    ```

    Alternatively, you can install themisto and minimap2 from github using the links above. Note that we have had difficulty compiling themisto from source. The v3.2.2 release works for us on MacOS and Linux.
    

## Running the Pipeline

### General notes on running the pipeline:

The pipeline quits at the end of each stage (progress is saved). Therefore, one has to run the pipeline anew to start the next stage. (This is helpful for debugging long computations).


### Running the pipeline in Test Mode:

The pipeline is currently set to **Test Mode**, so that a limited number of genomes are downloaded, so that the pipeline can be run locally on a laptop for testing.
One can do so with the following snippet:

```
conda activate PCNdb_env
cd src/
python PCN_pipeline.py
``` 

**To run the full pipeline, open `PCN_pipeline.py` in your favorite text editor and set `TEST_MODE = False` in line 25. Users *must* run the full pipeline on a Linux HPC with the SLURM job manager for the following reasons.**

1) Several stages of the pipeline submit thousands of HPC jobs in parallel to speed up computation, and this is not possible on a laptop.
2) The SRA data download is substantial– ~15TB of sequencing reads – and so a full download is not possible.  


```
conda activate PCNdb_env
cd src/
sbatch --mem=16G -t 430:00:00 --wrap="python PCN_pipeline.py"
```

You can submit to a partition specific to your lab as well, this is what we run on the Duke Compute Cluster in the You lab.

    sbatch --mem=16G -t 430:00:00 -p youlab --wrap="python PCN_pipeline.py"

One our our users uses the MIT compute cluster, and uses the following snippet to load conda, activate their environment, and run the pipeline:

    sbatch --mem=16G -t 00:00:10 -p mit_normal --wrap=“source /home/software/anaconda3/2023.07/etc/profile.d/conda.sh && conda activate PCNdb_env && python PCN_pipeline.py”


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

