FROM continuumio/miniconda3:latest

# Install dependencies that might be needed
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dev \
    wget \
    curl \
    libxml2 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Add conda-forge as a channel and install SRA toolkit
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y sra-tools && \
    conda clean -a -y

# Verify installation
RUN fastq-dump --version

WORKDIR /app

COPY . /app

RUN mkdir data
RUN mkdir results

RUN wget -O data/prokaryotes.txt https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

RUN (head -n 1 data/prokaryotes.txt && grep "plasmid" data/prokaryotes.txt | grep "Complete Genome" | sed 's/GCA/GCF/g') > results/complete-prokaryotes-with-plasmids.txt

# Create conda environment with Python 3.10 instead of cloning base (which has Python 3.12)
RUN conda create --name PCNdb_env python=3.10
# Initialize conda for bash (the default shell in most Docker images)
RUN conda init bash
# Use shell activation to run commands in the conda environment
SHELL ["conda", "run", "-n", "PCNdb_env", "/bin/bash", "-c"]

# Install packages in the PCNdb_env environment
RUN pip install pysradb biopython tqdm beautifulsoup4 polars HTSeq numpy
RUN conda install -c bioconda kallisto breseq
RUN conda install -c conda-forge ncbi-datasets-cli


WORKDIR /app/src

# # Reset shell to default and provide a command that activates the environment when the container starts
SHELL ["/bin/bash", "-c"]
CMD ["conda", "run", "-n", "PCNdb_env", "python3", "PCN_pipeline.py"]