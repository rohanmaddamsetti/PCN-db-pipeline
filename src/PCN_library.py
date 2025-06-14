#!/usr/bin/env python

"""
PCN_library.py by Rohan Maddamsetti, Maggie Wilson, and Irida Shyti.

For this pipeline to work, ncbi datasets, pysradb, kallisto, themisto, minimap2, and breseq 0.39+ must be in the $PATH.
On the Duke Compute Cluster (DCC), run the following to get these programs into the path:
conda activate PCNdb-env

IMPORTANT: this script assumes it is being run on DCC if sys.platform == "linux".
This means that users on a linux machine will need to modify a couple functions if they
are running this code locally, and cannot use slurm to submit many jobs in parallel.

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
import threading ## for keep track of stuck threads.
import logging
import glob
import pprint
import random
import polars as pl
from tqdm import tqdm
import HTSeq ## for filtering fastq multireads.
import numpy as np ## for matrix multiplications for running PIRA.
from bs4 import BeautifulSoup ## for parsing breseq output.

## imports for parallelization.
from tqdm.asyncio import tqdm as async_tqdm
import asyncio
import shutil
from contextlib import asynccontextmanager

"""
TODO list:

1) clean up code to be consistent throughout in the use of
AnnotationAccessions and/or RefSeq_IDs as directory names.

2) fix my polars dataframe code style to use parentheses, and start lines with ".join" and so forth.
Example:

filtered_df = (
    df.groupby("AnnotationAccession")
      .agg([
          (pl.col("SufficientReadCount").all()).alias("AllTrue")
      ])
      .filter(pl.col("AllTrue"))
      .join(df, on="AnnotationAccession", how="inner")
      .drop("AllTrue")
)

Additional notes: the following breseq runs did not finish within 24h with 16GB of memory and 1 core:
GCF_000025625.1_ASM2562v1
GCF_014872735.1_ASM1487273v1
"""


################################################################################
## Classes.

class RateLimiter:
    """Rate limiter to prevent overwhelming NCBI servers."""
    
    def __init__(self, calls_per_minute=10):
        self.calls_per_minute = calls_per_minute
        self.interval = 60.0 / calls_per_minute  # seconds between calls
        self.last_call_time = 0
        self.lock = asyncio.Lock()
    
    @asynccontextmanager
    async def limit(self):
        """Context manager to limit the rate of API calls."""
        async with self.lock:
            # Calculate time since last call
            current_time = time.time()
            time_since_last_call = current_time - self.last_call_time
            
            # If we need to wait to respect rate limit
            if time_since_last_call < self.interval:
                wait_time = self.interval - time_since_last_call
                logging.debug(f"Rate limiting: waiting {wait_time:.2f} seconds")
                await asyncio.sleep(wait_time)
            
            # Update last call time
            self.last_call_time = time.time()
        
        try:
            # Allow the caller to proceed
            yield
        finally:
            pass  # No cleanup needed


################################################################################
## Functions.

def get_SRA_ID_from_RefSeqID(refseq_id):
    """Fetch the SRA ID corresponding to a RefSeq accession ID."""
    logging.info(f"Getting SRA ID for {refseq_id}...")
    ## datasets must be in $PATH.
    bash_command = f"datasets summary genome accession {refseq_id}"

    try:
        cmd_output = subprocess.check_output(bash_command, shell=True, stderr=subprocess.STDOUT)
        json_output = cmd_output.decode("utf-8")
        json_data = json.loads(json_output)
        
        # Log the full response structure
        logging.info(f"Full NCBI response structure for {refseq_id}:")
        logging.info(json.dumps(json_data, indent=2))

        sra_id = "NA"
        reports = json_data.get("reports")
        if reports:
            sample_ids = reports[0].get("assembly_info", {}).get("biosample", {}).get("sample_ids")
            if sample_ids:
                for sample_id in sample_ids:
                    if sample_id.get("db") == "SRA":
                        sra_id = sample_id.get("value")
                        break
        return sra_id

    except subprocess.CalledProcessError as e:
        logging.error(f"Error fetching SRA ID for {refseq_id}: {e.output.decode()}")
        return "NA"


def fetch_Run_IDs_with_pysradb(sra_id):
    """Fetch Run IDs for a given SRA ID, filtering strictly for Illumina WGS data with retry logic."""
    logging.info(f"Fetching Run IDs for {sra_id}...")

    ## pysradb must be in $PATH.
    pysradb_command = ["pysradb", "metadata", "--detailed", sra_id]
    max_retries = 5
    attempt = 0
    pysra_command_worked = False

    while attempt < max_retries:
        try:
            result = subprocess.run(pysradb_command, capture_output=True, text=True, check=True)
            output = result.stdout.strip()
            if output:
                pysra_command_worked = True
                break
        except subprocess.CalledProcessError as e:
            logging.warning(f"Attempt {attempt+1}/{max_retries} failed for {sra_id}: {str(e)}")
        attempt += 1
        time.sleep(0.1)  # Small delay before retrying
        
    if not pysra_command_worked:
        logging.error(f"Failed to fetch metadata for {sra_id} after {max_retries} attempts.")
        return []

    # Parse the output to extract relevant metadata
    lines = output.split("\n")
    if len(lines) <= 1:  # No data, only header
        logging.warning(f"No run data found for SRA ID: {sra_id}")
        return []

    run_ids = list()
    header = lines[0].split("\t")
    for my_SRA_metadata in lines[1:]:
        fields = my_SRA_metadata.split("\t")
        try: ## if the sequencing run in this line of metadata has bad data, just try the next line.
            run_idx = header.index("run_accession")
            total_size_idx = header.index("total_size")
            run_id = fields[run_idx]
            total_size = fields[total_size_idx]
            int_total_size = int(total_size) ## this can fail if total_size == '<NA>'
        except ValueError as e:
            logging.error(f"Error parsing pysradb output: {str(e)}")
            continue

        ## If there is data associated with this accession (total_size > 0), the run_accession is valid,
        ## and this is Illumina WGS data, then add to the list of run_accessions.
        ## to be really safe, only use data with both ILLUMINA and Illumina keywords in its metadata.
        if int_total_size > 0 and run_id != "nan" and "WGS" in my_SRA_metadata and "GENOMIC" in my_SRA_metadata and "ILLUMINA" in my_SRA_metadata and "Illumina" in my_SRA_metadata:
            run_ids.append(run_id)

    if run_ids:
        logging.info(f"Found {len(run_ids)} WGS Illumina runs for SRA ID: {sra_id}: {', '.join(run_ids)}")
    else:
        logging.warning(f"No Illumina WGS data found for SRA ID: {sra_id}")
        
    return run_ids


def create_RefSeq_SRA_RunID_table(prokaryotes_with_plasmids_file, RunID_table_csv, TEST_MODE, TEST_DOWNLOAD_LIMIT):
    """Create a table mapping RefSeq IDs to SRA IDs and Run IDs."""
    logging.info("Creating RefSeq SRA RunID table...")
    
    # Create the output directory if it doesn't exist
    os.makedirs(os.path.dirname(RunID_table_csv), exist_ok=True)
    
    with open(RunID_table_csv, 'w') as f:
        writer = csv.writer(f)
        logging.info("Writing header to RunID table...")
        writer.writerow(["RefSeq_ID", "SRA_ID", "Run_ID"])
        
        # Count successful entries for reporting
        successful_entries = 0

        ## get all RefSeq IDs in the complete-prokaryotes-with-plasmids.txt file.
        with open(prokaryotes_with_plasmids_file, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            next(reader)  # Skip header
            
            for i, row in enumerate(reader):
                
                if TEST_MODE and successful_entries >= TEST_DOWNLOAD_LIMIT:
                    logging.info(f"Test mode: stopping calls once SRA Run IDs for {TEST_DOWNLOAD_LIMIT} genomes collected")
                    print(f"Test mode: stopping RefSeq_SRA_RunID table after {TEST_DOWNLOAD_LIMIT} genome entries")
                    break

                ## should be 23 tab-separated fields in complete-prokaryotes-with-plasmids.txt.
                assert len(row) == 23
                RefSeq_ID = row[18]
                ## filter for valid IDs (some rows have a '-' as a blank placeholder).
                if RefSeq_ID == "-":
                    continue
                assert RefSeq_ID.startswith("GCF") ## make sure we have the right ID
                logging.info(f"Processing genome {i+1}: {RefSeq_ID}")

                SRA_ID = get_SRA_ID_from_RefSeqID(RefSeq_ID)
                if SRA_ID != "NA":
                    logging.info(f"Found SRA ID: {SRA_ID} for {RefSeq_ID}")
                    Run_IDs = fetch_Run_IDs_with_pysradb(SRA_ID)
                    
                    if Run_IDs:
                        for Run_ID in Run_IDs:
                            if Run_ID == "0" or Run_ID == "nan":
                                continue ## skip a bad Run_ID.
                            writer.writerow([RefSeq_ID, SRA_ID, Run_ID])
                            successful_entries += 1
                            logging.info(f"Added entry: {RefSeq_ID}, {SRA_ID}, {Run_ID}")
                    else:
                        logging.warning(f"No Run IDs found for {RefSeq_ID} (SRA ID: {SRA_ID})")
                else:
                    logging.warning(f"No SRA ID found for {RefSeq_ID}")
    
    ## Log summary
    logging.info(f"Finished creating RunID table with {successful_entries} entries")
    
    ## Debug: output the contents of the file
    if exists(RunID_table_csv):
        with open(RunID_table_csv, 'r') as f:
            content = f.read()
            logging.info(f"RunID table contents:\n{content}")
    else:
        logging.error(f"RunID table file {RunID_table_csv} was not created")
    return


def create_refseq_accession_to_ftp_path_dict(prokaryotes_with_plasmids_file):
    refseq_accession_to_ftp_path_dict = dict()
    
    with open(prokaryotes_with_plasmids_file, "r") as prok_with_plasmids_file_obj:
        for i, line in enumerate(prok_with_plasmids_file_obj):
            if i == 0: continue  # skip the header
            
            fields = line.strip().split("\t")                
            # Get the accession field (5th from end) and turn GCA Genbank IDs into GCF RefSeq IDs
            refseq_id = fields[-5].replace("GCA", "GCF")
            
            # Get the ftp_url field (3rd from end) and make sure we turn the GCA Genbank URL
            # into the GCF RefSeq FTP URL
            ftp_url = fields[-3].replace("GCA", "GCF")
            
            # Check for valid IDs and URLs (suppressed RefSeq IDs have a '-' as a blank placeholder)
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

    if sys.platform == "darwin": ## then the checksum is the LAST 'word' in the output.
        my_md5_checksum = md5_call.stdout.split()[-1].strip()
    elif sys.platform == "linux": ## then the checksum is the FIRST 'word' in the output.
        my_md5_checksum = md5_call.stdout.split()[0].strip()
    else:
        raise AssertionError("UNKNOWN PLATFORM")

    ## verify that the checksums match.
    return my_md5_checksum == my_target_checksum


async def download_single_genome(ftp_path, reference_genome_dir):
    """Download a single genome and its MD5 file"""
    my_full_accession = basename(ftp_path)
    my_base_filename = my_full_accession + "_genomic.gbff.gz"
    
    # Files on the NCBI FTP site to download
    gbff_ftp_path = os.path.join(ftp_path, my_base_filename)
    md5_ftp_path = os.path.join(ftp_path, "md5checksums.txt")
    
    # Local paths
    gbff_gz_file = os.path.join(reference_genome_dir, my_base_filename)
    md5_file = os.path.join(reference_genome_dir, my_full_accession + "_md5checksums.txt")

    # Check if files exist and are valid
    if exists(gbff_gz_file) and exists(md5_file):
        if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
            print(f"{gbff_gz_file} SUCCEEDED.")
            return True
        else:
            os.remove(gbff_gz_file)
            os.remove(md5_file)

    # Try downloading up to 5 times
    attempts = 5
    while attempts > 0:
        try:
            await asyncio.to_thread(urllib.request.urlretrieve, gbff_ftp_path, filename=gbff_gz_file)
            await asyncio.to_thread(urllib.request.urlretrieve, md5_ftp_path, filename=md5_file)
            
            if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
                print(f"{gbff_gz_file} SUCCEEDED.")
                return True
            else:
                if exists(gbff_gz_file):
                    os.remove(gbff_gz_file)
                if exists(md5_file):
                    os.remove(md5_file)
                
        except Exception as e:
            print(f"Attempt {6-attempts} failed for {gbff_gz_file}: {str(e)}")
            if exists(gbff_gz_file):
                os.remove(gbff_gz_file)
            if exists(md5_file):
                os.remove(md5_file)
        
        attempts -= 1
    
    print(f"{gbff_gz_file} FAILED after all attempts")
    return False


async def async_download(ftp_paths, reference_genome_dir, max_concurrent=10):
    """Download genomes in parallel using asyncio task pool with rate limiting"""
    ## Create a semaphore to limit concurrent downloads
    ## Adjust based on your network capacity
    semaphore = asyncio.Semaphore(max_concurrent)
    
    # Create a rate limiter for NCBI
    rate_limiter = RateLimiter(calls_per_minute=20)  # Adjust as needed
    
    async def download_with_limits(ftp_path):
        """Download a single genome with rate and concurrency limits"""
        async with semaphore:  # Limit concurrent downloads
            async with rate_limiter.limit():  # Respect rate limits
                return await download_single_genome(ftp_path, reference_genome_dir)
    
    # Create tasks for all downloads
    tasks = []
    for ftp_path in ftp_paths:
        task = asyncio.create_task(download_with_limits(ftp_path))
        tasks.append(task)
    
    # Use tqdm to show progress
    results = []
    for task in tqdm(asyncio.as_completed(tasks), total=len(tasks)):
        result = await task
        results.append(result)
    
    return results


def fetch_reference_genomes(RunID_table_file, refseq_accession_to_ftp_path_dict, reference_genome_dir):
    """Download reference genomes for each genome in the RunID table"""
    ## Create reference genome directory if it doesn't exist
    os.makedirs(reference_genome_dir, exist_ok=True)

    ## Get RefSeq IDs from the RunID table
    with open(RunID_table_file, "r") as RunID_file_obj:
        RunID_table_lines = RunID_file_obj.read().splitlines()

    ## Remove the header from the imported data
    RunID_table_data = RunID_table_lines[1:]
    ## Get the first column to get all refseq_ids of interest
    ## Set comprehension to remove duplicates (there can be multiple SRA datasets per reference genome)
    refseq_ids = {line.split(",")[0] for line in RunID_table_data}
    
    ## Look up the FTP URLs for each refseq id
    ## IMPORTANT: we have to check to see if the ftp path exists; this is not true for suppressed entries in RefSeq,
    ## which have an ftp_path == '-'
    ftp_paths = [refseq_accession_to_ftp_path_dict[x] for x in refseq_ids if x in refseq_accession_to_ftp_path_dict]

    ## Run the async download
    asyncio.run(async_download(ftp_paths, reference_genome_dir))
    return


def get_Run_IDs_from_RunID_table(RunID_table_csv):
    """Get Run IDs from the RunID table CSV file."""
    Run_IDs = []
    with open(RunID_table_csv, 'r') as f:
        reader = csv.reader(f)
        next(reader)  ## Skip header
        for row in reader:
            ## The Run_ID column is the 3rd row (index 2)
            ## remove whitespace with strip, if needed.
            Run_IDs.append(row[2].strip())
    logging.info(f"Found {len(Run_IDs)} Run IDs in {RunID_table_csv}")    
    return Run_IDs


async def prefetch_fastq_reads(run_id, SRA_data_dir, max_retries=5):
    """prefetch FASTQ files for a run ID with retry logic."""
    logging.info(f"Prefetching {run_id} with command: prefetch --max-size 100G -O {SRA_data_dir} {run_id}")
    
    for attempt in range(max_retries):
        try:
            # First prefetch the SRA data
            logging.info(f"Prefetching {run_id} (attempt {attempt+1}/{max_retries})...")
            prefetch_cmd = ["prefetch", "--max-size", "100G", "-O", SRA_data_dir, run_id]
            prefetch_process = await asyncio.create_subprocess_exec(
                *prefetch_cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            prefetch_stdout, prefetch_stderr = await prefetch_process.communicate()
            
            if prefetch_process.returncode != 0:
                logging.error(f"Error prefetching {run_id} (attempt {attempt+1}/{max_retries}): {prefetch_stderr.decode()}")
                await asyncio.sleep(1)
                continue
            else:
                logging.info(f"Successfully prefetched {run_id}")
                return True
            
        except Exception as e:
            logging.error(f"Exception prefetching {run_id} (attempt {attempt+1}/{max_retries}): {str(e)}")
            await asyncio.sleep(1)
    
    logging.error(f"Failed to prefetch {run_id} after {max_retries} attempts")
    return False


async def prefetch_fastq_reads_parallel(SRA_data_dir, Run_IDs, max_concurrent=10, max_retries=5):
    """prefetch FASTQ reads in parallel with retry logic and integrity checks."""
    os.makedirs(SRA_data_dir, exist_ok=True)
    semaphore = asyncio.Semaphore(max_concurrent)

    async def prefetch_with_semaphore(run_id):
        async with semaphore:
            return await prefetch_fastq_reads(run_id, SRA_data_dir, max_retries)

    total_downloads = len(Run_IDs)
    completed = 0
    failed = 0

    tasks = [prefetch_with_semaphore(run_id) for run_id in Run_IDs]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    
    # Count successful downloads
    success_count = sum(1 for r in results if r is True)
    logging.info(f"Successfully prefetched {success_count}/{len(Run_IDs)} Run IDs")
    
    return success_count


async def validate_sra_download(run_id, SRA_data_dir):
    """Validate the integrity of downloaded SRA data using vdb-validate."""
    logging.info(f"Validating SRA download for {run_id}...")
    
    try:
        # Run vdb-validate on the SRA accession
        cmd = ["vdb-validate", os.path.join(SRA_data_dir, run_id)]
        
        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )
        
        stdout, stderr = await process.communicate()
        stdout_text = stdout.decode() if stdout else ""
        stderr_text = stderr.decode() if stderr else ""
        
        ## Check if validation was successful
        if process.returncode == 0 and "is consistent" in stderr_text:
            logging.info(f"SRA validation successful for {run_id}")
            return True
        else:
            error_msg = stderr_text.strip() or stdout_text.strip() or f"vdb-validate returned code {process.returncode}"
            logging.warning(f"SRA validation failed for {run_id}: {error_msg}")
    except Exception as e:
        logging.error(f"Error validating SRA download for {run_id}: {str(e)}")
    return False


async def validate_sra_download_parallel(SRA_data_dir, Run_IDs, max_concurrent=10):
    """Run vdb-validate in parallel with retry logic and integrity checks."""
    semaphore = asyncio.Semaphore(max_concurrent)

    async def validate_with_semaphore(run_id):
        async with semaphore:
            return await validate_sra_download(run_id, SRA_data_dir)

    total_downloads = len(Run_IDs)
    completed = 0
    failed = 0

    tasks = [validate_with_semaphore(run_id) for run_id in Run_IDs]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    # Count validated downloads
    success_count = sum(1 for r in results if r is True)
    logging.info(f"Successfully downloaded and validated {success_count}/{len(Run_IDs)} Run IDs")

    if success_count == len(results):
        return True ## all SRA files were validated.
    return False


async def unpack_fastq_reads(run_id, SRA_data_dir, max_retries=5):
    """Unpack prefetched sra files into FASTQ."""
    logging.info(f"Unpacking {run_id} with command: fasterq-dump --split-files --skip-technical --temp {SRA_data_dir} --outdir {SRA_data_dir} --progress {run_id}")

    for attempt in range(max_retries):
        try:
            ## Use fasterq-dump to extract FASTQ
            logging.info(f"Running fasterq-dump for {run_id}...")
            cmd = [
                "fasterq-dump",
                "--split-files",
                "--skip-technical",
                "--temp", SRA_data_dir,
                "--outdir", SRA_data_dir,
                "--progress",
                os.path.join(SRA_data_dir, run_id) ## give the path to the prefetched dirs.
            ]

            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )

            stdout, stderr = await process.communicate()

            if process.returncode != 0:
                logging.error(f"Error unpacking {run_id} (attempt {attempt+1}/{max_retries}): {stderr.decode()}")
                await asyncio.sleep(1)
                continue

            else:
                logging.info(f"Successfully unpacked {run_id} into fastq data.")
                return True

        except Exception as e:
            logging.error(f"Exception unpacking {run_id} (attempt {attempt+1}/{max_retries}): {str(e)}")
            await asyncio.sleep(1)
    
    logging.error(f"Failed to fasterq-dump {run_id}")
    return False


async def unpack_fastq_reads_parallel(SRA_data_dir, Run_IDs, max_concurrent=10, max_retries=5):
    """unpack FASTQ reads in parallel with retry logic and integrity checks."""
    os.makedirs(SRA_data_dir, exist_ok=True)
    semaphore = asyncio.Semaphore(max_concurrent)

    async def unpack_with_semaphore(run_id):
        async with semaphore:
            return await unpack_fastq_reads(run_id, SRA_data_dir)

    total_downloads = len(Run_IDs)
    completed = 0
    failed = 0

    tasks = [unpack_with_semaphore(run_id) for run_id in Run_IDs]
    results = await asyncio.gather(*tasks, return_exceptions=True)
    # Count successful downloads
    success_count = sum(1 for r in results if r is True)
    logging.info(f"Successfully unpacked {success_count}/{len(Run_IDs)} Run IDs")

    if success_count == len(results):
        return True ## all SRA files were unpacked.
    return False


def all_fastq_data_exist(run_ids, SRA_data_dir):
    """Check if all FASTQ files exist for the given run IDs.
    Does not check for corrupted data."""
    missing_run_ids = []
    
    for run_id in run_ids:
        files_found = False

        paired_files = [
            os.path.join(SRA_data_dir, f"{run_id}_1.fastq"),
            os.path.join(SRA_data_dir, f"{run_id}_2.fastq")
        ]
        single_file = os.path.join(SRA_data_dir, f"{run_id}.fastq")

        ## handle a couple weird corner cases from SRA.
        if run_id == "SRR12988984":
            paired_files = [
                os.path.join(SRA_data_dir, f"{run_id}_3.fastq"),
                os.path.join(SRA_data_dir, f"{run_id}_4.fastq")
            ]
        elif run_id == "SRR7353055":
            single_file = os.path.join(SRA_data_dir, f"{run_id}_1.fastq"),

        ## now check to see if the fastq files for this run_id exist on disk.
        if (os.path.exists(paired_files[0]) and os.path.exists(paired_files[1])) or os.path.exists(single_file):
            files_found = True
            break
            
        if not files_found:
            missing_run_ids.append(run_id)
    
    if missing_run_ids:
        logging.warning(f"Missing FASTQ files for {len(missing_run_ids)}/{len(run_ids)} Run IDs")
        logging.warning(f"Missing Run IDs: {', '.join(missing_run_ids)}...")
        return False
    return True


def generate_replicon_level_fasta_reference_for_kallisto(gbk_gz_path, outfile):
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
                ## Important: for kallisto, we need to replace spaces with underscores in the replicon annotation field.
                replicon_description = record.description.replace(" ","_")
                header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"replicon="+replicon_description])
                outfh.write(header + "\n")
                outfh.write(str(record.seq) + "\n")
    return


def make_benchmark_fasta_refs_for_kallisto(refgenomes_dir, kallisto_ref_outdir, PIRA_low_PCN_benchmark_csv_file):
    ## this function makes a genome fasta file for each genome.
    ## each genome fasta file contains fasta sequences for every replicon.

    ## make the directory if it doesn't exist yet.
    os.makedirs(kallisto_ref_outdir, exist_ok=True)
    
    ## we only run the genomes in the low PCN benchmarking list.
    benchmark_genome_df = pl.read_csv(PIRA_low_PCN_benchmark_csv_file)
    ## get rid of duplicates with this logic
    benchmark_genome_list = sorted(list(set(benchmark_genome_df["AnnotationAccession"].to_list())))
    
    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff.gz")]
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        genome_id = gzfile.split(".gbff.gz")[0].strip("_genomic")
        
        ## skip genome_ids that are not in the list of benchmark genomes.
        if genome_id not in benchmark_genome_list:
            continue

        fasta_outfile = os.path.join(kallisto_ref_outdir, genome_id+".fna")
        generate_replicon_level_fasta_reference_for_kallisto(gzpath, fasta_outfile)
    return


def make_kallisto_indices(kallisto_ref_dir, kallisto_index_dir):
    ## make the directory if it doesn't exist yet.
    os.makedirs(kallisto_index_dir, exist_ok=True)
    
    ref_fasta_filelist = [x for x in os.listdir(kallisto_ref_dir) if x.endswith(".fna")]
    for ref_fasta_file in ref_fasta_filelist:
        ref_fasta_path = os.path.join(kallisto_ref_dir, ref_fasta_file)
        genome_id = ref_fasta_file.split(".fna")[0]
        index_file = genome_id + ".idx"
        index_path = os.path.join(kallisto_index_dir, index_file)
        kallisto_index_args = ["kallisto", "index", "-i", index_path, ref_fasta_path]
        subprocess.run(kallisto_index_args)
    return


def run_kallisto_quant(RunID_table_csv, kallisto_index_dir, SRA_data_dir, results_dir):
    ## make the output directory if it doesn't exist.
    os.makedirs(results_dir, exist_ok=True)
    ## IMPORTANT: kallisto needs -l -s parameters supplied when run on single-end data.
    ## to avoid this complexity, I only process paired-end Illumina data, and skip single-end data altogether.
    RefSeq_to_SRA_RunList_dict = make_RefSeq_to_SRA_RunList_dict(RunID_table_csv)
    
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
            slurm_string = "sbatch -p scavenger --mem=16G --wrap=" + "\"" + kallisto_quant_string + "\""
            if sys.platform == "linux": ## assume that we are running on DCC
                print("sys.platform == 'linux' so we assume this script is being run on the Duke Compute Cluster")
                print(slurm_string)
                subprocess.run(slurm_string, shell=True)
            else:
                print("sys.platform != 'linux' so we assume this script is being run on a mac laptop")
                print(kallisto_quant_string)
                subprocess.run(kallisto_quant_string, shell=True)
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


def parse_replicon_metadata_in_header(target_id):
    fields = target_id.split("|")
    SeqID = fields[0].split("=")[-1]
    SeqType = fields[1].split("=")[-1]
    ## convert underscores back into spaces.
    replicon_description = fields[2].split("=")[-1].replace("_", " ")
    metadata_tuple = (SeqID, SeqType, replicon_description)
    return(metadata_tuple)


def estimate_replicon_copy_numbers(kallisto_replicon_count_tsv_path):

    chromosomal_length = 0.0
    chromosomal_est_counts = 0.0

    replicon_coverage_dict = dict()
    ## get the chromosomal gene coverage, and get the coverage for all genes
    with open(kallisto_replicon_count_tsv_path, "r") as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0: continue ## skip header
            target_id, length, eff_length, est_counts, tpm = line.split("\t")
            SeqID, SeqType, replicon_description = parse_replicon_metadata_in_header(target_id)
            coverage = float(est_counts) / float(length)
            replicon_coverage_dict[SeqID] = (SeqType, replicon_description, coverage)
            if SeqType == "chromosome":
                chromosomal_length += float(length)
                chromosomal_est_counts += float(est_counts)
    ## Return an empty dict() when nothing aligns to the chromosome.
    if chromosomal_length == 0:
        print("WARNING: no reads pseudoaligned to chromosome in file: ", kallisto_replicon_count_tsv_path)
        print("estimate_replicon_copy_numbers is returning an empty dict.")
        return(dict())
    chromosome_coverage = chromosomal_est_counts / chromosomal_length
    ## now normalize by chromosome coverage to get copy number estimates.
    replicon_copy_number_dict = dict()
    for my_SeqID, value_tuple in replicon_coverage_dict.items():
        my_SeqType, replicon_description, my_coverage = value_tuple
        my_replicon_copy_number = my_coverage / chromosome_coverage
        replicon_copy_number_dict[my_SeqID] = (my_SeqType, replicon_description, my_replicon_copy_number)
    return(replicon_copy_number_dict)


def measure_kallisto_replicon_copy_numbers(kallisto_replicon_quant_results_dir, replicon_copy_number_csv_file):
    """
    define lists to encode the following columns of the table.
    RefSeqID, SeqID, SeqType, CopyNumber
    """
    RefSeqIDVec = []
    SeqIDVec = [] ## this is for the replicon.
    SeqTypeVec = []
    RepliconDescriptionVec = []
    CopyNumberVec = []
    
    ## skip .DS_Store and any other weird files.
    genomedirectories = [x for x in os.listdir(kallisto_replicon_quant_results_dir) if x.startswith("GCF")]
    for genomedir in genomedirectories:
        refseq_id = "_".join(genomedir.split("_")[:2])
        genome_quantfile_path = os.path.join(kallisto_replicon_quant_results_dir, genomedir, "abundance.tsv")
        replicon_copy_number_dict = estimate_replicon_copy_numbers(genome_quantfile_path)
        for SeqID, value_tuple in replicon_copy_number_dict.items():
            seqtype, replicon_description, copy_number = value_tuple
            RefSeqIDVec.append(refseq_id)
            SeqIDVec.append(SeqID)
            SeqTypeVec.append(seqtype)
            RepliconDescriptionVec.append(replicon_description)
            CopyNumberVec.append(str(copy_number))

    assert len(RefSeqIDVec) == len(SeqIDVec) == len(SeqTypeVec) == len(RepliconDescriptionVec) == len(CopyNumberVec)

    ## we have to double-quote all columns-- some fields in the product column contain commas!
    RefSeqIDVec = ["\"" + x + "\"" for x in RefSeqIDVec]
    SeqIDVec = ["\"" + x + "\"" for x in SeqIDVec]
    SeqTypeVec = ["\"" + x + "\"" for x in SeqTypeVec]
    RepliconDescriptionVec = ["\"" + x + "\"" for x in RepliconDescriptionVec]
    CopyNumberVec = ["\"" + x + "\"" for x in CopyNumberVec]
    
    ## now write the replicon copy number data to file.
    with open(replicon_copy_number_csv_file, "w") as outfh:
        ## double-quote each column name in the header for consistency.
        header = "\"RefSeqID\",\"SeqID\",\"SeqType\",\"RepliconDescription\",\"CopyNumber\""
        outfh.write(header + "\n")
        for i in range(len(RefSeqIDVec)):
            outfh.write(RefSeqIDVec[i] + "," + SeqIDVec[i] + "," + SeqTypeVec[i] + "," + RepliconDescriptionVec[i] + "," + CopyNumberVec[i] + "\n")
    return


def tabulate_NCBI_replicon_lengths(refgenomes_dir, replicon_length_csv_file):
    with open(replicon_length_csv_file, 'w') as outfh:
        header = "AnnotationAccession,SeqID,SeqType,replicon_length\n"
        outfh.write(header)
        for gbk_gz in os.listdir(refgenomes_dir):
            if not gbk_gz.endswith(".gbff.gz"): continue
            annotation_accession = gbk_gz.split("_genomic.gbff")[0]
            infile = os.path.join(refgenomes_dir, gbk_gz)
            with gzip.open(infile, "rt") as genome_fh:
                for i, replicon in enumerate(SeqIO.parse(genome_fh, "gb")):
                    SeqID = replicon.id
                    if "chromosome" in replicon.description or i == 0:
                        ## IMPORTANT: we assume here that the first record is a chromosome.
                        SeqType = "chromosome"
                    elif "plasmid" in replicon.description:
                        SeqType = "plasmid"
                    else:
                        continue
                    replicon_length = str(len(replicon))
                    ## now write out the data for the replicon.
                    row = ','.join([annotation_accession, SeqID, SeqType, replicon_length])
                    outfh.write(row + "\n")
    return


def generate_replicon_fasta_references_for_themisto(gbk_gz_path, fasta_outdir):
    print("reading in as input:", gbk_gz_path)
    ## open the input reference genome file.
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
            ## replace spaces with underscores in the replicon annotation field.
            replicon_description = record.description.replace(" ","_")
            header = ">" + "|".join(["SeqID="+SeqID,"SeqType="+SeqType,"replicon="+replicon_description])
            my_replicon_fastafile = SeqID + ".fna"
            my_replicon_outfilepath = os.path.join(fasta_outdir, my_replicon_fastafile)
            with open(my_replicon_outfilepath, "w") as outfh:
                outfh.write(header + "\n")
                outfh.write(str(record.seq) + "\n")
    return


def generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir):
    """
    IMPORTANT: the replicons in the fasta reference list file need to be sorted by length.
    Directly measure the length of each of the replicons here, and then use this to sort
    the list of replicon FASTA paths, before writing to file.
    """
    genome_id = os.path.basename(fasta_outdir)
    ## IMPORTANT: exclude any FASTA file that has the genome_id as a name.
    ## (this is a file containing all the replicons, created and used in the PIRA stages downstream).
    replicon_fasta_filelist = [x for x in os.listdir(fasta_outdir) if x.endswith(".fna") and not x.startswith(genome_id)]
    replicon_fasta_pathlist = [os.path.join(fasta_outdir, x) for x in replicon_fasta_filelist]

    ## Decorate-sort-undecorate by fasta sequence length.
    decorated_replicon_fasta_pathlist = list()
    for fastapath in replicon_fasta_pathlist:
        my_replicon = SeqIO.read(fastapath, "fasta")
        ## Get the length of the replicon
        replicon_length = len(my_replicon.seq)
        ## append a tuple of (fastapath, replicon_length)
        my_tuple = (fastapath, replicon_length)
        decorated_replicon_fasta_pathlist.append(my_tuple)
    ## sort the decorated list in place by replicon_length, in descending order from
    ## from largest to smallest replicon.
    decorated_replicon_fasta_pathlist.sort(key=lambda x: x[1], reverse=True)
    ## undecorate the path list, which is now in sorted order.
    sorted_replicon_fasta_pathlist = [x[0] for x in decorated_replicon_fasta_pathlist]

    ## write the path list to file for themisto.
    replicon_listfile = os.path.join(fasta_outdir, genome_id + ".txt")
    with open(replicon_listfile, "w") as fastatxtfile_fh:
        for fasta_path in sorted_replicon_fasta_pathlist:
            fastatxtfile_fh.write(fasta_path + "\n")
    return


def make_NCBI_replicon_fasta_refs_for_themisto(refgenomes_dir, themisto_fasta_ref_outdir):
    ## this function makes a genome directory for each genome.
    ## each directory contains separate fasta files for each replicon.

    ## make the output directory if it does not exist.
    if not exists(themisto_fasta_ref_outdir):
        os.makedirs(themisto_fasta_ref_outdir)

    gzfilelist = [x for x in os.listdir(refgenomes_dir) if x.endswith("gbff.gz")]
    for gzfile in gzfilelist:
        gzpath = os.path.join(refgenomes_dir, gzfile)
        ## get rid of the useless "_genomic" suffix here.
        genome_id = gzfile.split(".gbff.gz")[0].replace("_genomic", "")
        fasta_outdir = os.path.join(themisto_fasta_ref_outdir, genome_id)
        ## make the fasta output directory if it does not exist.
        if not exists(fasta_outdir):
            os.makedirs(fasta_outdir)
        generate_replicon_fasta_references_for_themisto(gzpath, fasta_outdir)
        generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir)
    return


def run_command_and_retry_if_it_fails(command_string, tempdir=None, max_retries=3, timeout=90):
    ## This code handles a bug in themisto build-- sometimes randomly hangs, have to delete temp files
    ## and restart and then it usually works.
    retries = 0
    while retries < max_retries:
        process = subprocess.Popen(command_string, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        timer = threading.Timer(timeout, process.kill)  ## Kill process if it exceeds timeout

        try:
            timer.start()
            stdout, stderr = process.communicate()
        finally:
            timer.cancel()

        if process.returncode == 0:
            print("Command succeeded:", stdout.decode())
            return stdout.decode()
        else:
            print(f"*********COMMAND FAILED (attempt {retries + 1}):", stderr.decode())
            if tempdir is not None: ## remove temporary files from the failed run.
                print(f"removing {tempdir}")
                subprocess.run(f"rm -rf {tempdir}", shell=True)
                print(f"remaking {tempdir} before restarting")
                os.makedirs(tempdir)
            retries += 1
            time.sleep(0.1)  ## Small delay before retrying
    
    print("Command failed after maximum retries.")
    return


def make_themisto_indices(themisto_ref_dir, themisto_index_dir):
    ## make the output directory if it does not exist.
    if not exists(themisto_index_dir):
        os.makedirs(themisto_index_dir)

    ## each directory is named after the genome_id of the given genome.
    for genome_id in os.listdir(themisto_ref_dir):
        ## get the actual path for this directory
        ref_fasta_dir = os.path.join(themisto_ref_dir, genome_id)
        ## make sure that this path is real, and not an artifact of some weird file like .DS_Store in this directory
        if not os.path.isdir(ref_fasta_dir):
            continue
        index_input_filelist = os.path.join(ref_fasta_dir, genome_id + ".txt")

        ## make the directory for the index files, if it does not exist.
        genome_index_dir = os.path.join(themisto_index_dir, genome_id)
        if not exists(genome_index_dir):
            os.makedirs(genome_index_dir)
        ## set the index_prefix to write index files into the genome_index_dir.
        index_prefix = os.path.join(genome_index_dir, genome_id)
        
        ## make the temp directory if it doesn't exist.
        tempdir = os.path.join(genome_index_dir, "temp")
        if not exists(tempdir):
            os.makedirs(tempdir)
        
        themisto_build_args = ["themisto", "build", "-k","31", "-i", index_input_filelist, "--index-prefix", index_prefix, "--temp-dir", tempdir, "--mem-gigas", "8", "--n-threads", "8", "--file-colors"]

        themisto_build_string = " ".join(themisto_build_args)
        print(themisto_build_string)
        run_command_and_retry_if_it_fails(themisto_build_string, tempdir)
            
    return


def run_themisto_pseudoalign(RunID_table_csv, themisto_index_dir, SRA_data_dir, themisto_pseudoalignment_dir):
    ## make the output directory if it does not exist.
    if not exists(themisto_pseudoalignment_dir):
        os.makedirs(themisto_pseudoalignment_dir)

    ## make this data structure
    RefSeq_to_SRA_RunList_dict = make_RefSeq_to_SRA_RunList_dict(RunID_table_csv)
        
    ## each directory in themisto_index_dir is named after the genome_id of the given genome.
    genome_id_list = [x for x in os.listdir(themisto_index_dir)]
    for genome_id in genome_id_list:
        my_index_dir = os.path.join(themisto_index_dir, genome_id)
        ## make sure that this path is real, and not an artifact of some weird file in this directory
        if not os.path.isdir(my_index_dir):
            continue
        ## make the arguments relating to the paths to the index files for this genome.
        my_index_prefix = os.path.join(my_index_dir, genome_id)

        ## make the file containing the paths to the sequencing read data for this genome.
        refseq_id = "_".join(genome_id.split("_")[:2])

        ## check for inconsistencies between Run_ID_list and the genomes with themisto indices
        if refseq_id not in RefSeq_to_SRA_RunList_dict.keys():
            print(f"INCONSISTENCY ERROR: {refseq_id} not found in RunID_table.csv!")
            continue
        
        Run_ID_list = RefSeq_to_SRA_RunList_dict[refseq_id]
        ## make read_path_arg_list.
        readpath_list = list()
        for Run_ID in Run_ID_list:
            SRA_file_pattern = f"{SRA_data_dir}/{Run_ID}*.fastq"
            matched_fastq_list = sorted(glob.glob(SRA_file_pattern))
            readpath_list += matched_fastq_list

        ## if we didn't find any fastq data for this genome, then skip to the next genome.
        if not len(readpath_list):
            continue
        
        ## now, write the paths for the SRA read data to disk for themisto pseudoalign.
        SRAdata_listfile = os.path.join(themisto_pseudoalignment_dir, genome_id + "_SRAdata.txt")
        with open(SRAdata_listfile, "w") as SRAtxtfile_fh:
            for readpath in readpath_list:
                SRAtxtfile_fh.write(readpath + "\n")

        ## make the directory for the pseudoalignments.
        my_pseudoalignment_output_dir = os.path.join(themisto_pseudoalignment_dir, genome_id)
        if not exists(my_pseudoalignment_output_dir):
            os.makedirs(my_pseudoalignment_output_dir)

        ## make the temp directory if it doesn't exist.
        tempdir = os.path.join(my_pseudoalignment_output_dir, "temp")
        if not exists(tempdir):
            os.makedirs(tempdir)
        
        ## we make corresponding pseudoalignment output files for each SRA dataset.
        ## This list goes into the output listfile.
        output_listfile = os.path.join(themisto_pseudoalignment_dir, genome_id + "_pseudoalignments.txt")
        with open(output_listfile, "w") as output_listfile_fh:
            for readpath in readpath_list:
                read_filename = os.path.basename(readpath).split(".fastq")[0]
                output_filename = os.path.join(my_pseudoalignment_output_dir, read_filename + "_pseudoalignment.txt")
                output_listfile_fh.write(output_filename + "\n")
        
        ## now run themisto pseudoalign.
        themisto_pseudoalign_args = ["themisto", "pseudoalign", "--query-file-list", SRAdata_listfile, "--index-prefix", my_index_prefix, "--temp-dir", tempdir, "--out-file-list", output_listfile, "--n-threads", "4", "--threshold", "0.7"]
        themisto_pseudoalign_string = " ".join(themisto_pseudoalign_args)
        ## don't run this command with retries!
        slurm_string = "sbatch -p scavenger --cpus-per-task=4 --mem=16G --wrap=" + "\"" + themisto_pseudoalign_string + "\""
        if sys.platform == "linux": ## assume that we are running on DCC
            print("sys.platform == 'linux' so we assume this script is being run on the Duke Compute Cluster")
            print(slurm_string)
            subprocess.run(slurm_string, shell=True)
        else:
            print("sys.platform != 'linux' so we assume this script is being run on a mac laptop")
            print(themisto_pseudoalign_string)
            subprocess.run(themisto_pseudoalign_string, shell=True)
    return


def map_themisto_IDs_to_replicon_metadata(themisto_replicon_ref_dir, my_genome_dirname):
    ## now, let's map the themisto replicon ID numbers to a (SeqID, SeqType) tuple.
    themisto_ID_to_seq_metadata_dict = dict()
    my_themisto_replicon_ID_mapping_file = os.path.join(themisto_replicon_ref_dir, my_genome_dirname, my_genome_dirname + ".txt")
    with open(my_themisto_replicon_ID_mapping_file, "r") as replicon_ID_mapping_fh:
        for i,fasta_file_path in enumerate(replicon_ID_mapping_fh):
            fasta_file_path = fasta_file_path.strip()
            ## get the header from this fasta file.
            with open(fasta_file_path, "r") as fasta_fh:
                my_header = fasta_fh.readline().strip()
                my_seq = fasta_fh.readline().strip()
            fields = my_header.strip(">").split("|")
            SeqID = fields[0].split("=")[-1]
            SeqType = fields[1].split("=")[-1]
            SeqLength = len(my_seq)
            themisto_ID_to_seq_metadata_dict[i] = (SeqID, SeqType, SeqLength)
    return themisto_ID_to_seq_metadata_dict


def summarize_themisto_pseudoalignment_results(themisto_replicon_ref_dir, themisto_pseudoalignment_dir, themisto_results_csvfile_path):
    
    with open(themisto_results_csvfile_path, "w") as output_csv_fh:
        ## first, write the output header.
        output_header = "AnnotationAccession,SeqID,SeqType,ReadCount"
        output_csv_fh.write(output_header + "\n")
        
        print(output_header)
        ## Iterate over the directories in themisto_pseudoalignment_dir.
        ## These contain the pseudoalignments for each genome.
        themisto_pseudoalignment_result_dirs = [x for x in os.listdir(themisto_pseudoalignment_dir) if x.startswith("GCF") and os.path.isdir(os.path.join(themisto_pseudoalignment_dir, x))]
        for AnnotationAccession in themisto_pseudoalignment_result_dirs:
            if not AnnotationAccession.startswith("GCF"): continue ## just an additional check to remove the temp directory.
            my_cur_pseudoalignment_dir_path = os.path.join(themisto_pseudoalignment_dir, AnnotationAccession)
            
            ## initialize a dictionary to store the pseudoalignment counts.
            pseudoalignment_read_count_dict = dict()
            pseudoalignment_filepaths = [os.path.join(my_cur_pseudoalignment_dir_path, x) for x in os.listdir(my_cur_pseudoalignment_dir_path) if x.endswith("_pseudoalignment.txt")]
            ## now, summarize the read counts for this genome.
            for my_pseudoalignment_filepath in pseudoalignment_filepaths:
                with open(my_pseudoalignment_filepath) as my_pseudoalignment_fh:
                    for line in my_pseudoalignment_fh:
                        ## handle odd behavior in themisto: we need to sort the themisto replicon ID numbers ourselves.
                        ## IMPORTANT: sort the themisto ids (strings) by their numerical value.
                        replicon_set_string = " ".join(sorted(line.strip().split()[1:], key=lambda x: int(x)))
                        if replicon_set_string in pseudoalignment_read_count_dict:
                            pseudoalignment_read_count_dict[replicon_set_string] += 1
                        else:
                            pseudoalignment_read_count_dict[replicon_set_string] = 1

            ## now, let's map the themisto replicon ID numbers to a (SeqID, SeqType) tuple.
            themisto_ID_to_seq_metadata_dict = map_themisto_IDs_to_replicon_metadata(themisto_replicon_ref_dir, AnnotationAccession)
            
            ## now write the pseudoalignment counts to file.
            for replicon_set_string in sorted(pseudoalignment_read_count_dict.keys()):
                read_count = pseudoalignment_read_count_dict[replicon_set_string]
                replicon_ID_list = replicon_set_string.split()
                if len(replicon_ID_list) == 0:
                    SeqID = "NA"
                    SeqType = "NA"
                elif len(replicon_ID_list) == 1:
                    my_replicon_ID = replicon_ID_list[0]
                    ## the SeqLength parameter in themisto_ID_to_seq_metadata_dict is not used here.
                    SeqID, SeqType, _ = themisto_ID_to_seq_metadata_dict[int(my_replicon_ID)]
                else:
                    SeqID = "&".join([themisto_ID_to_seq_metadata_dict[int(replicon_ID)][0] for replicon_ID in replicon_ID_list])
                    SeqType = "multireplicon_sequence"
                ## now write to file.
                rowdata = ",".join([AnnotationAccession, SeqID, SeqType, str(read_count)])
                print(rowdata)
                output_csv_fh.write(rowdata + "\n")
    return


def naive_themisto_PCN_estimation(themisto_results_csv_file, replicon_length_csv_file, naive_themisto_PCN_csv_file):
    ## This function simply ignores multireplicon reads when estimating PCN.
    ##Also note that this function omits replicons with zero mapped reads.
    print("running naive themisto PCN estimation (ignoring multireplicon reads)")
    ## import the data as polars dataframes.
    
    replicon_length_df = pl.read_csv(replicon_length_csv_file)
    naive_themisto_read_count_df = pl.read_csv(themisto_results_csv_file).filter(
        (pl.col("SeqType") == "chromosome") | (pl.col("SeqType") == "plasmid")).join(
            replicon_length_df, on = "SeqID").with_columns(
                (pl.col("ReadCount").cast(pl.Float64) / pl.col("replicon_length").cast(pl.Float64)).alias("SequencingCoverage"))
    
    ## make a second dataframe containing just the sequencing coverage for the longest replicon for each genome.
    ## to do so, first group by AnnotationAccession and compute maximum replicon_length within each group.
    longest_replicon_df = naive_themisto_read_count_df.group_by(
        "AnnotationAccession").agg(pl.col("replicon_length").max()).join(
            ## now join with the original DataFrame to filter for rows with the maximum replicon_length
            naive_themisto_read_count_df, on=["AnnotationAccession", "replicon_length"], how="inner").select(
                pl.col("AnnotationAccession", "SequencingCoverage")).with_columns(
                    pl.col("SequencingCoverage").alias("LongestRepliconCoverage")).select(
                        pl.col("AnnotationAccession", "LongestRepliconCoverage"))

    ## now normalize SequencingCoverage by LongestRepliconCoverage for each genome to calculate PCN.
    naive_themisto_PCN_df = naive_themisto_read_count_df.join(
        longest_replicon_df, on = "AnnotationAccession", coalesce=True).with_columns(
            (pl.col("SequencingCoverage") / pl.col("LongestRepliconCoverage")).alias("CopyNumber")).select(
                ## super annoying, not sure why the AnnotationAccession_right and SeqType_right columns
                ## are kept. This select command removes these redundancies.
                ['AnnotationAccession', 'SeqID', 'SeqType', 'ReadCount', 'replicon_length', 'SequencingCoverage',
                 'LongestRepliconCoverage','CopyNumber'])

    ## now write the naive PCN estimates to file.
    naive_themisto_PCN_df.write_csv(naive_themisto_PCN_csv_file)    
    return


def assign_multireplicon_reads(genome_df):
    ## create a new data frame with just chromosome and plasmid data.
    updated_df = genome_df.filter(
        (pl.col("SeqType") == "chromosome") | (pl.col("SeqType") == "plasmid")).with_columns(
            ## cast the ReadCount to floats, and
            ## replace null values in the ReadCount columns with floating-point zeros.
            pl.col("ReadCount").cast(pl.Float64).fill_null(pl.lit(0.0)))

    ## get the multireplicon data.
    multireplicon_df = genome_df.filter(pl.col("SeqType") == "multireplicon_sequence")
    ## iterate over each multireplicon read set.
    for row_dict in multireplicon_df.iter_rows(named=True):
        seq_id_list = row_dict["SeqID"].split("&")
        read_count = row_dict["ReadCount"]
        num_reads_to_assign_to_each_replicon = float(read_count) / float(len(seq_id_list))
        ## iterate over each replicon in this multireplicon read set.
        for seq_id in seq_id_list:
            ## update the relevant values in the dataframe
            old_seqID_readcount = updated_df.filter(updated_df["SeqID"] == seq_id)["ReadCount"][0]
            new_seqID_readcount = old_seqID_readcount + num_reads_to_assign_to_each_replicon
            ## a kludgy hacky solution, but works:
            ## create a new column called temp, that has the new_seqID_readcount for the given SeqID,
            ## and zeros elsewhere in the column.
            ## then make a new column called updatedReadCount that takes the maximum of ReadCount and temp.
            ## then drop ReadCount and temp and rename updatedReadCount as ReadCount.
            updated_df = updated_df.with_columns(
                "ReadCount",
                pl.when(pl.col("SeqID") == seq_id).then(new_seqID_readcount).otherwise(pl.lit(0.0)).alias("temp")
            ).with_columns(pl.max_horizontal("ReadCount", "temp").alias("updatedReadCount")).select(
                pl.col("*").exclude(["ReadCount", "temp"])).rename({"updatedReadCount":"ReadCount"})
    return updated_df


def make_gbk_annotation_table(reference_genome_dir, gbk_annotation_file):
    with open(gbk_annotation_file,"w") as out_fh:
        header = "Annotation_Accession,host,isolation_source\n"
        out_fh.write(header)

        gbk_files = [x for x in os.listdir(reference_genome_dir) if x.endswith("_genomic.gbff.gz")]

        for x in tqdm(gbk_files):
            gbk_path = reference_genome_dir + x
            annotation_accession = x.split("_genomic.gbff.gz")[0]
            with gzip.open(gbk_path,'rt') as gbk_fh:
                host = "NA"
                isolation_source = "NA"
                ## use a buffer to store the line, to handle cases where
                ## the annotation spans multiple lines.
                in_host_field = False
                in_isolation_source_field = False
                line_buffer = []
                for line in gbk_fh:
                    line = line.strip()
                    ## We're going to delete all double-quote characters,
                    ## and replace all commas with semicolons so that they
                    ## don't clash with the csv format.
                    ## BUT: make sure that the line terminates with a double-quote--
                    ## otherwise, we need to slurp up data from multiple lines.
                    if line.startswith("/host"):
                        line_annot = line.split('=')[-1].replace('\"','').replace(',',';')
                        if line.endswith('"'):
                            host = line_annot
                        else: ## have to look at the next line too.
                            line_buffer.append(line_annot)
                            in_host_field = True
                    elif in_host_field and len(line_buffer):
                        line_annot = line.replace('\"','').replace(',',';')
                        line_buffer.append(line_annot)
                        if line.endswith('"'): ## flush the line buffer and reset the flag.
                            host = ' '.join(line_buffer)
                            line_buffer = []
                            in_host_field = False
                    elif line.startswith("/isolation_source"):
                        line_annot = line.split('=')[-1].replace('\"','').replace(',',';')
                        if line.endswith('"'):
                            isolation_source = line_annot
                        else: ## then have to look at next line too.
                            line_buffer.append(line_annot)
                            in_isolation_source_field = True
                    elif in_isolation_source_field and len(line_buffer):
                        line_annot = line.replace('\"','').replace(',',';')
                        line_buffer.append(line_annot)
                        if line.endswith('"'): ## flush the line buffer and reset flag.
                            isolation_source = ' '.join(line_buffer)
                            line_buffer = []
                            in_isolation_source_field = False
                    ## break if we have the annotation.
                    if (host != "NA") and (isolation_source != "NA"): break
                    ## also break if we're looking at gene annotation since
                    ## there's insufficient isolate annotation in this file.
                    if line.startswith("/gene"): break ## NOTE: has to be '/gene'
                    ## because we don't want to match the Genome Annotation Data in the
                    ## COMMENT field of the Genbank metadata.
                    if line.startswith("ORIGIN"): break ## break if we got to the first fasta seq
                row = ','.join([annotation_accession, host, isolation_source]) + '\n'
                out_fh.write(row)
    return


def filter_fastq_files_for_multireads(multiread_data_dir, themisto_pseudoalignment_dir, SRA_data_dir):

    ## make the directory for filtered multireads  if it does not exist.
    if not exists(multiread_data_dir):
        os.makedirs(multiread_data_dir)

    files_in_pseudoalignment_dir = [x for x in os.listdir(themisto_pseudoalignment_dir)]
    paths_in_pseudoalignment_dir = [os.path.join(themisto_pseudoalignment_dir, x) for x in files_in_pseudoalignment_dir]
    pseudoalignment_dirpaths = [x for x in paths_in_pseudoalignment_dir if os.path.isdir(x)]

    for my_pseudoalignment_results_dir in pseudoalignment_dirpaths:
        ## make sure to ignore the temp directory by filtering by filename suffixes!
        my_pseudoalignment_files = [x for x in os.listdir(my_pseudoalignment_results_dir) if x.endswith("_pseudoalignment.txt")]
        my_pseudoalignment_paths = [os.path.join(my_pseudoalignment_results_dir, x) for x in my_pseudoalignment_files]

        for cur_pseudoalignment_path in my_pseudoalignment_paths:
            ## before sorting the pseudoalignment results, check whether the filtered multireads
            ## for this pseudoalignment already exist on disk. If so, then skip--
            ## don't repeat work that has already been done.

            ## construct the path to the original fastq file.
            my_fastq_file = basename(cur_pseudoalignment_path).split("_pseudoalignment.txt")[0] + ".fastq"
            my_fastq_path = os.path.join(SRA_data_dir, my_fastq_file)

            ## construct the path to the filtered fastq file.
            my_genome = basename(my_pseudoalignment_results_dir)
            ## note: multiread_genome_dir is only created if filtered multireads exist.
            ## this means that this path only exists on disk if filtered multireads for this genome
            ## were already written to disk in a previous call of this function.
            multiread_genome_dir = os.path.join(multiread_data_dir, my_genome)
            my_filtered_fastq_file = "multireads_" + my_fastq_file
            my_filtered_fastq_path = os.path.join(multiread_genome_dir, my_filtered_fastq_file)
            ## now check to see if the filtered fastq file already exists.
            if exists(my_filtered_fastq_path): continue ## if so, then don't repeat the work.
            
            list_of_multiread_tuples = list()
            with open(cur_pseudoalignment_path, "r") as pseudoalignment_fh:
                for line in pseudoalignment_fh:
                    line = line.strip() ## remove whitespace
                    fields = line.split()
                    if len(fields) <= 2: continue ## skip reads that map to zero or one reference sequences.
                    read_index = int(fields[0])
                    matched_reference_tuple = tuple(fields[1:])
                    multiread_tuple = (read_index, matched_reference_tuple)
                    list_of_multiread_tuples.append(multiread_tuple)
            ## sort the list of multiread tuples in-place by the read_index.
            list_of_multiread_tuples.sort(key=lambda x: x[0])
            ## if there are multireads, then filter the original fastq file for multireads.
            if len(list_of_multiread_tuples):

                ## if there are multireads, then make multiread_genome_dir for this genome.
                if not exists(multiread_genome_dir):
                    os.makedirs(multiread_genome_dir)
                
                ## IMPORTANT: read indices in the themisto pseudoalignments are zero-based (first index is 0).
                multiread_indices = {x[0] for x in list_of_multiread_tuples} ## this is a set
                print(f"filtering {my_fastq_file} for multireads. multireads written into {my_filtered_fastq_path}")
                ## write out the filtered reads
                with open(my_filtered_fastq_path, "w") as filtered_fastq_fh:
                    my_fastq_reader = HTSeq.FastqReader(my_fastq_path)
                    for i, read in enumerate(my_fastq_reader):
                        ## skip reads that aren't in the set of multireads
                        if i not in multiread_indices: continue
                        read.write_to_fastq_file(filtered_fastq_fh)          
    return


def make_fasta_reference_genomes_for_minimap2(themisto_replicon_ref_dir):
    paths_in_themisto_replicon_ref_dir = [os.path.join(themisto_replicon_ref_dir, x) for x in os.listdir(themisto_replicon_ref_dir)]
    ## extra checking to make sure we are handling directories.
    dirpaths_in_themisto_replicon_ref_dir = [x for x in paths_in_themisto_replicon_ref_dir if os.path.isdir(x)]
    for my_dir in dirpaths_in_themisto_replicon_ref_dir:
        ## get the paths to each replicon FASTA file.
        ## IMPORTANT: the order of replicons in this file MATTERS.
        ## This corresponds to the zero-indexed replicon assignment done by themisto pseudoalignment.
        replicon_fasta_pathlist = []
        my_replicon_fasta_listfilepath = [os.path.join(my_dir, x) for x in os.listdir(my_dir) if x.endswith(".txt")].pop()
        with open(my_replicon_fasta_listfilepath, "r") as replicon_listfile_fh:
            for line in replicon_listfile_fh:
                line = line.strip()
                replicon_fasta_pathlist.append(line)

        my_fasta_reference_genome_outfile = basename(my_dir) + ".fna"
        my_fasta_reference_genome_outpath = os.path.join(my_dir, my_fasta_reference_genome_outfile)

        with open(my_fasta_reference_genome_outpath, "w") as fasta_outfh:
            for themisto_replicon_num, replicon_fasta_path in enumerate(replicon_fasta_pathlist):
                with open(replicon_fasta_path, "r") as my_fasta_infh:
                    for i, line in enumerate(my_fasta_infh):
                        ## add the replicon ID number assigned by themisto to the FASTA header
                        if i == 0:
                            header_string = line.lstrip(">")
                            replicon_id_string = "ThemistoRepliconID=" + str(themisto_replicon_num)
                            updated_header = ">" +  replicon_id_string + "|" + header_string
                            fasta_outfh.write(updated_header)
                        else:
                            fasta_outfh.write(line)
    return


def align_reads_for_benchmark_genomes_with_minimap2(
        PIRA_low_PCN_benchmark_csv_file, RunID_table_csv,
        themisto_replicon_ref_dir, SRA_data_dir, benchmark_alignment_dir):
    ## Example: minimap2 -x sr ref.fa read1.fq read2.fq > aln.paf
    ## Details of PAF format are here: https://github.com/lh3/miniasm/blob/master/PAF.md

    ## make the directory for read alignments  if it does not exist.
    if not exists(benchmark_alignment_dir):
        os.makedirs(benchmark_alignment_dir)

    ## get the AnnotationAccessions of interest for benchmarking.
    benchmark_genomes_df = pl.read_csv(PIRA_low_PCN_benchmark_csv_file)
    
    ## get the unique RefSeq_IDs in this dataframe
    AnnotationAccession_list = list(set(benchmark_genomes_df.get_column("AnnotationAccession").to_list()))
    ## the following is a hack to split AnnotationAccession on the second occurrence of an "_" to get the RefSeq_ID.
    AnnotationAccession_to_RefSeq_ID_dict = {x : "_".join(x.split("_", 2)[:2]) for x in AnnotationAccession_list}

    ## we will get the fastq files for each genome from this DataFrame.
    benchmark_RunID_table_df = pl.read_csv(RunID_table_csv).filter(
        pl.col("RefSeq_ID").is_in(AnnotationAccession_to_RefSeq_ID_dict.values()))
    
    for annotation_accession in AnnotationAccession_list:
        
        ## make a subdirectory for the output alignments.
        genome_alignment_dir = os.path.join(benchmark_alignment_dir, annotation_accession)

        if not exists(genome_alignment_dir):
            os.makedirs(genome_alignment_dir)

        ref_genome_fasta_file = annotation_accession + ".fna"
        reference_genome_path = os.path.join(themisto_replicon_ref_dir, annotation_accession, ref_genome_fasta_file)
        
        my_RefSeq_ID = AnnotationAccession_to_RefSeq_ID_dict[annotation_accession]
        my_RunID_df = benchmark_RunID_table_df.filter(pl.col("RefSeq_ID") == my_RefSeq_ID)
        my_RunID_list = my_RunID_df.get_column("Run_ID").to_list()

        ## Now use this list of Run IDs to pattern match for *.fastq files for this genome.
        my_read_data_pathlist = list()

        for Run_ID in my_RunID_list:
            SRA_file_pattern = f"{SRA_data_dir}/{Run_ID}*.fastq"
            matched_fastq_list = sorted(glob.glob(SRA_file_pattern))
            my_read_data_pathlist += matched_fastq_list
        my_read_data_pathlist.sort() ## sort the read data paths for this genome.

        ## run single-end read alignment on each fastq file separately.
        for my_index, cur_read_data_path in enumerate(my_read_data_pathlist):
            my_alignment_file = basename(cur_read_data_path).split(".fastq")[0] + ".paf"
            my_alignment_outpath = os.path.join(genome_alignment_dir, my_alignment_file)

            minimap2_cmd_string = " ".join(["minimap2 -x sr", reference_genome_path, cur_read_data_path, ">", my_alignment_outpath])
            print(minimap2_cmd_string)
            subprocess.run(minimap2_cmd_string, shell=True)
    return


def align_multireads_with_minimap2(themisto_replicon_ref_dir, multiread_data_dir, multiread_alignment_dir):
    ## Example: minimap2 -x sr ref.fa read1.fq read2.fq > aln.paf
    ## Details of PAF format are here: https://github.com/lh3/miniasm/blob/master/PAF.md

    ## make the directory for multiread alignments  if it does not exist.
    if not exists(multiread_alignment_dir):
        os.makedirs(multiread_alignment_dir)
    
    genomes_with_multireads = [x for x in os.listdir(multiread_data_dir) if x.startswith("GCF")]
    for my_genome in genomes_with_multireads:
        
        ## make a subdirectory for the multiread alignments.
        multiread_genome_alignment_dir = os.path.join(multiread_alignment_dir, my_genome)
        if not exists(multiread_genome_alignment_dir):
            os.makedirs(multiread_genome_alignment_dir)
        
        ref_genome_fasta_file = my_genome + ".fna"
        reference_genome_path = os.path.join(themisto_replicon_ref_dir, my_genome, ref_genome_fasta_file)
        my_multiread_data_dir = os.path.join(multiread_data_dir, my_genome)
        my_multiread_data_pathlist = [os.path.join(my_multiread_data_dir, x) for x in os.listdir(my_multiread_data_dir)]
        my_multiread_data_pathlist.sort() ## sort the filtered data.
        
        ## run single-end read alignment on each fastq file separately.
        for my_index, cur_multiread_data_path in enumerate(my_multiread_data_pathlist):
            my_alignment_file = basename(cur_multiread_data_path).split(".fastq")[0] + ".paf"
            my_multiread_alignment_outpath = os.path.join(multiread_genome_alignment_dir, my_alignment_file)

            minimap2_cmd_string = " ".join(["minimap2 -x sr", reference_genome_path, cur_multiread_data_path, ">", my_multiread_alignment_outpath])
            print(minimap2_cmd_string)
            subprocess.run(minimap2_cmd_string, shell=True)
    return


def parse_read_alignments(genome_dir):
    ## make a dictionary from reads to the multiset of replicons that the read maps to.
    paf_alignment_files = [x for x in os.listdir(genome_dir) if x.endswith(".paf")]
    paf_alignment_paths = [os.path.join(genome_dir, x) for x in paf_alignment_files]
    
    read_mapping_dict = dict()
    for paf_path in paf_alignment_paths:
        with open(paf_path, "r") as paf_fh:
            for line in paf_fh:
                line = line.strip()
                fields = line.split("\t")
                ## see PAF format documentation here:
                ## https://github.com/lh3/miniasm/blob/master/PAF.md
                read_name = fields[0]
                themisto_replicon_ID = int(fields[5].split("|")[0].split("=")[-1]) ## IMPORTANT: turn into an integer
                
                if read_name in read_mapping_dict:
                    read_mapping_dict[read_name].append(themisto_replicon_ID)
                else:
                    read_mapping_dict[read_name] = [themisto_replicon_ID]
    return read_mapping_dict


def initialize_GenomeDataFrame(themisto_ID_to_seq_metadata_dict):
    """ I use themisto_ID_to_seq_metadata_dict to initialize the dataframe for a genome.
     themisto_ID_to_seq_metadata_dict contains metadata for all replicons in the genome,
     regardless of whether any reads mapped to that replicon.
     To handle the cases that:
     1) a replicon only contains multireads, therefore it's naive PCN == 0, or
     2) a replicon does not have any multireads, therefore, it is not present in additional_replicon_reads_dict, or
     3) no replicons have any unireads mapping to them, so my_naive_themisto_PCN_df is empty.
     we use themisto_ID_to_seq_metadata_dict to initialize a basic polars dataframe containing
     all replicons in our genome.
    """
    themisto_id_col = sorted(themisto_ID_to_seq_metadata_dict.keys())
    replicon_seq_id_col = list()
    replicon_seq_type_col = list()
    replicon_length_col = list()

    for themisto_id in themisto_id_col:
        ## get the SeqID and SeqType given the Themisto ID.
        seq_id, seq_type, seq_length = themisto_ID_to_seq_metadata_dict[themisto_id]
        replicon_seq_id_col.append(seq_id)
        replicon_seq_type_col.append(seq_type)
        replicon_length_col.append(seq_length)

    genome_df = pl.DataFrame(
        {"ThemistoID" : themisto_id_col,
         "SeqID" : replicon_seq_id_col,
         "SeqType" : replicon_seq_type_col,
         "replicon_length" : replicon_length_col})

    return genome_df


def make_PIRAGenomeDataFrame(
        additional_replicon_reads_dict,
        themisto_ID_to_seq_metadata_dict,
        my_naive_themisto_PCN_df):
    """ Make the DataFrame with the data needed for PIRA on a given genome.
        We have to update the results of the Naive PCN estimates from Themisto
        (results of stage 16) by adding the additional replicon reads found by re-aligning multireads
        with minimap2.
    """

    ## Turn the additional_replicon_reads_dict into a polars dataframe.
    ## First initialize the dataframe to contain all replicons in the genome.
    additional_replicon_reads_df = initialize_GenomeDataFrame(themisto_ID_to_seq_metadata_dict)
    ## initialize the AdditionalReadCount Column with zeros for each row in the initialized DataFrame
    additional_readcount_col = [0 for i in range(additional_replicon_reads_df.shape[0])]        
    ## now update the values in additional_readcount_col using additional_replicon_reads_dict.
    for themisto_id, readcount in additional_replicon_reads_dict.items():
        additional_readcount_col[themisto_id] = readcount
    
    ## now add the AdditionalReadCount Column to the DataFrame.
    additional_readcount_series = pl.Series("AdditionalReadCount", additional_readcount_col)
    additional_replicon_reads_df = additional_replicon_reads_df.with_columns([additional_readcount_series])

    ## IMPORTANT: the previous code ensures that additional_replicon_reads_df
    ## contains ALL replicons in the genome, even if no reads pseudoaligned or aligned to that replicon.
    ## this may happen for short contigs that are annotated as plasmids or plasmid fragments in the
    ## reference genome (example: genome GCF_017654585.1_ASM1765458v1).

    ## To fill in missing values in the AnnotationAccession column after the merge with additional_replicon_reads_df,
    ## get the unique value in the AnnotationAccession column in my_naive_themisto_PCN_df, and use this
    ## to fill in the missing values.
    my_AnnotationAccession = my_naive_themisto_PCN_df.select(pl.col('AnnotationAccession').unique())
    assert len(my_AnnotationAccession) == 1 ## make sure this value is unique!
    
    ## merge the DataFrames containing the ReadCounts.
    merged_readcount_df = additional_replicon_reads_df.join(
        ## IMPORTANT: my_naive_themisto_PCN_df may not contain rows for replicons that didn't have
        ## any reads pseudoalign to it. therefore, we need to left_join my_native_themisto_PCN_df to
        ## additional_replicon_reads_df, which contains the data for ALL replicons in the genome,
        ## even if the Count is zero.
        my_naive_themisto_PCN_df, on="SeqID", how="left", coalesce=True).with_columns(
            ##fill in missing values in the AnnotationAccession column after the merge.
            pl.col("AnnotationAccession").fill_null(my_AnnotationAccession)).with_columns(
                    ## set missing values in the InitialReadCount column to 0.
                    pl.col("InitialReadCount").fill_null(strategy="zero")).with_columns(
                        ## sum those ReadCounts,
                        (pl.col("InitialReadCount") + pl.col("AdditionalReadCount")).alias("ReadCount")).with_columns(
                            ## and re-calculate SequencingCoverage,
                            (pl.col("ReadCount").cast(pl.Float64) / pl.col("replicon_length").cast(pl.Float64)).alias("SequencingCoverage"))
    
    ## The following is a hack, following code in the function native_themisto_PCN_estimation(),
    ## to recalculate LongestRepliconCoverage and CopyNumber with the additional reads.

    ## find the length of the longest replicon
    max_replicon_length = merged_readcount_df.select(pl.col("replicon_length").max()).item()

    ## filter the DataFrame to get the row for the longest replicon
    longest_replicon_row_df = merged_readcount_df.filter(
        pl.col("replicon_length") == max_replicon_length).with_columns(
            ## and define the LongestRepliconCoverage column for merging back.
            pl.col("SequencingCoverage").alias("LongestRepliconCoverage")).select(
            pl.col("AnnotationAccession", "LongestRepliconCoverage"))

    ## now normalize SequencingCoverage by LongestRepliconCoverage for each genome to calculate PCN.
    PIRAGenomeDataFrame = merged_readcount_df.join(
        longest_replicon_row_df, on = "AnnotationAccession", coalesce=True).with_columns(
        (pl.col("SequencingCoverage") / pl.col("LongestRepliconCoverage")).alias("InitialCopyNumberEstimate")).sort(
            ## and sort by the ThemistoID column.
            "ThemistoID").select(
                ## super annoying, not sure why the AnnotationAccession_right columns
                ## are kept. This select command removes these redundancies.
                ['ThemistoID', 'SeqID', 'SeqType', 'replicon_length',
                 'AdditionalReadCount', 'AnnotationAccession', 'InitialReadCount',
                 'ReadCount', 'SequencingCoverage', 'LongestRepliconCoverage', 'InitialCopyNumberEstimate']
            )
    
    """
    Stage 12 calls generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir), which
    ensures that Themisto Replicon IDs are sorted by replicon length in descending order
    (i.e., Replicon ID == 0 corresponds to the longest chromosome).
    """
    return PIRAGenomeDataFrame


def initializePIRA(multiread_mapping_dict, themisto_ID_to_seq_metadata_dict, my_naive_themisto_PCN_df):
    ## Iterate over the multiread_mapping_dict, and split into two data structures:
    ## 1) reads that map to a single replicon are counted up in a dictionary.
    additional_replicon_reads_dict = dict()
    ## 2) reads that map to multiple replicons are stored in a list of lists,
    ## that will then be turned into the numpy array to store the match matrix M.
    match_matrix_list_of_rows = list()

    for read, replicon_list in multiread_mapping_dict.items():
        replicon_set = set(replicon_list)
        if len(replicon_set) == 0: ## the read does not map to any replicons -- should not occur.
            raise AssertionError(f"READ {read} DID NOT ALIGN TO ANY REPLICON")
        elif len(replicon_set) == 1: ## the read maps to a single replicon.
            my_replicon = replicon_set.pop()
            if my_replicon in additional_replicon_reads_dict:
                additional_replicon_reads_dict[my_replicon] += 1
            else:
                additional_replicon_reads_dict[my_replicon] = 1
        else: ## the read maps to multiple replicons
            ## initialize a row of zeros, based on the number of replicons in this genome.
            match_matrix_rowlist = [0 for k in themisto_ID_to_seq_metadata_dict.keys()]
            for replicon_index in replicon_list:
                match_matrix_rowlist[replicon_index] += 1
            match_matrix_list_of_rows.append(match_matrix_rowlist)

    ## now set up the data structures for PIRA on this genome.
    MatchMatrix = np.array(match_matrix_list_of_rows)
    
    """ Generate the DataFrame containing the ReadCounts,  replicon lengths, and initial PCN estimates.
    We update the results of the Naive PCN estimates from Themisto (results of stage 16)
    by adding the additional replicon reads found by re-aligning multireads
    with minimap2.
    """
    PIRAGenomeDataFrame = make_PIRAGenomeDataFrame(
        additional_replicon_reads_dict, themisto_ID_to_seq_metadata_dict, my_naive_themisto_PCN_df)
    ## return inputs requires for PIRA for this genome.
    return (MatchMatrix, PIRAGenomeDataFrame)


def run_PIRA(M, PIRAGenomeDataFrame, epsilon = 0.00001):
    ## Run PIRA for a genome, assuming that the zero-th index of the match matrix M is the chromosome for normalization.
    print("RUNNING PIRA.")
    print(M)
    print(M.shape)
    print()

    print(PIRAGenomeDataFrame)
    print()
    print(PIRAGenomeDataFrame.glimpse(max_items_per_column=100))
    print()

    """
    Stage 12 calls generate_replicon_fasta_reference_list_file_for_themisto(fasta_outdir), which
    ensures that Themisto Replicon IDs are sorted by replicon length in descending order
    (i.e., Replicon ID == 0 corresponds to the longest chromosome).

    Therefore, we can  extract the InitialCopyNumber Column as a numpy vector as the initial PCN estimate guess for PIRA.
    """
    
    v = PIRAGenomeDataFrame["InitialCopyNumberEstimate"].to_numpy()
    readcount_vector = PIRAGenomeDataFrame["ReadCount"].to_numpy()
    replicon_length_vector = PIRAGenomeDataFrame["replicon_length"].to_numpy()
    
    """
    M may be empty, in the case that all multireads uniquely mapped to chromosome or plasmid.
    In this case, the PIRAGenomeDataFrame has already incorporated all available information,
    so we should just return v, the vector containing the InitialCopyNumberEstimate.

    I handle this logic by only running the PIRA loop if M is non-empty.
    """ 
    
    if M.shape[0] > 0: ## only run PIRA if the Match Matrix M has rows.
        convergence_error = 10000.0
        ## Iterate PIRA until error converges to zero.
        while convergence_error > epsilon:
            print(f"current convergence error: {convergence_error}")
            print(f"current PCN estimate vector: {v}")
            ## Weight M by PCN guess v-- need to turn v into a diagonal matrix first.
            diagonal_v = np.diag(v)
            weighted_M = np.matmul(M, diagonal_v)
            print(weighted_M)

            ## Normalize rows of weighted_M to sum to 1: this the probabilistic read assignment.
            ## Compute the sum of each row
            weighted_M_row_sums = weighted_M.sum(axis=1)
            ## Normalize each row by its sum
            normalized_weighted_M = weighted_M / weighted_M_row_sums[:, np.newaxis]

            print(normalized_weighted_M)
            ## sum over the rows of the normalized and weighted M matrix to generate the
            ## multiread vector.
            multiread_vector = normalized_weighted_M.sum(axis=0)
            print(multiread_vector)

            ## update the PCN estimate vector v using the multiread vector
            unnormalized_v = (multiread_vector + readcount_vector) / replicon_length_vector
            normalization_factor = unnormalized_v[0] ## coverage for the longest chromosome.
            updated_v = unnormalized_v / normalization_factor

            ## update the error estimate
            convergence_error = abs(sum(updated_v - v))
            ## and update the PCN estimate vector v
            v = updated_v
        print(f"final convergence error: {convergence_error}")

    print(f"final PCN estimate vector: {v}")
    return v


def run_PIRA_test_suite():
    ## Test 1: a simple test to check for convergence from very bad initial PCN estimates.
    Test1_DataFrame = pl.DataFrame(
        {
            "replicon_length" : [1000000, 100000, 1000], ## log 10: 6, 5, 3
            "ReadCount" : [1000000, 1000000, 1000000], ## log 10: 6, 6, 6
            "InitialCopyNumberEstimate" : [1, 1, 1] ## should converge to: 1, 10, 1000
        }
    )

    ## Test 1: the Match Matrix is negligible.
    multi_read_row1 = [1, 1, 0]
    test1_match_matrix_list_of_rows = [multi_read_row1]
    Test1_M = np.array(test1_match_matrix_list_of_rows)
    print("*"*80)
    print("PIRA TEST 1: check convergence from bad initial PCN estimates given negligible match matrix")
    run_PIRA(Test1_M, Test1_DataFrame)
    print()


    ## Test 2-- check what happens if the initial PCN vector contains zeros-- can estimates be updated
    ## properly?
    Test2_DataFrame = pl.DataFrame(
        {
            "replicon_length" : [1000000, 100000, 1000], ## log 10: 6, 5, 3
            "ReadCount" : [1000000, 1000000, 1000000], ## log 10: 6, 6, 6
            "InitialCopyNumberEstimate" : [1, 0, 0] ## should converge to: 1, 10, 1000
        }
    )

    ## Test 2: the Match Matrix is negligible, and the initial PCN estimate contains zeros
    print("*"*80)
    print("PIRA TEST 2: check convergence from bad initial PCN estimates with zeros, given negligible match matrix")
    run_PIRA(Test1_M, Test2_DataFrame)
    print()
    
    ## Test 3: Case that Match Matrix is super important for correct PCN estimation.
    Test3_DataFrame = pl.DataFrame(
        {
            "replicon_length" : [1000000, 100000, 1000], ## log 10: 6, 5, 3
            "ReadCount" : [1000000, 0, 0], ## log 10: 6, 6, 6
            "InitialCopyNumberEstimate" : [1, 1, 1] ## should converge to: 1, 10, 1000
        }
    )

    ## Test 3: the Match Matrix is super important for accurate PCN estimation.
    multi_read_row3 = [0, 1, 1]
    test3_match_matrix_list_of_rows = [multi_read_row3] * 1000000
    Test3_M = np.array(test3_match_matrix_list_of_rows)
    print("*"*80)
    print("PIRA TEST 3: check convergence when match matrix needed for accurate PCN estimation")
    run_PIRA(Test3_M, Test3_DataFrame)
    print()
    
    return


def run_PIRA_on_all_genomes(multiread_alignment_dir, themisto_replicon_ref_dir, naive_themisto_PCN_csv_file, PIRA_PCN_csv_file):
    
    ## only run PIRA on genomes with multireads.
    genomes_with_multireads = [x for x in os.listdir(multiread_alignment_dir) if x.startswith("GCF")]
    if len(genomes_with_multireads) == 0:
        print("No genomes in these data have multireads-- PIRA not run since it is not needed.")
        quit()
    
    ## import the results of the Naive PCN estimates from Themisto,
    ## and filter for rows corresponding to genomes with multireads.
    all_naive_themisto_PCN_estimates_df = pl.read_csv(naive_themisto_PCN_csv_file).filter(
        pl.col("AnnotationAccession").is_in(genomes_with_multireads))

    ## Make an empty Polars DataFrame to contain all the results.
    all_PIRA_estimates_DataFrame = pl.DataFrame()
    
    ## now populate the all_PIRA_estimates_DataFrame.
    for genome_ID in genomes_with_multireads:   
        ## get the Naive PCN estimates for this particular genome.

        ## map the themisto replicon ID numbers to a (SeqID, SeqType) tuple.
        themisto_ID_to_seq_metadata_dict = map_themisto_IDs_to_replicon_metadata(themisto_replicon_ref_dir, genome_ID)

        ## IMPORTANT: This DataFrame will NOT contain rows for replicons that didn't have
        ## any reads pseudoalign to it, and in some cases, may be completely empty.
        ##therefore, we need to join it to an initialized DataFrame
        ## which contains the data for ALL replicons in the genome, even if the Count is zero.
        my_naive_themisto_PCN_estimates_subset_df = all_naive_themisto_PCN_estimates_df.filter(
            pl.col("AnnotationAccession") == genome_ID).select(
                ## select only the columns we need.
                ['AnnotationAccession', 'SeqID', 'SeqType', 'ReadCount', 'replicon_length']
            ).rename({"ReadCount": "InitialReadCount"}) ## rename ReadCount to InitialReadCount

        ## Initialize a dataframe containing SeqID, SeqType, replicon_length metadata for all
        ## replicons in this genome, 
        my_naive_themisto_PCN_df = initialize_GenomeDataFrame(themisto_ID_to_seq_metadata_dict)
        ## Add a column for AnnotationAccession.
        annotation_accession_col = [genome_ID for i in range(my_naive_themisto_PCN_df.shape[0])]
        AnnotationAccession_series = pl.Series("AnnotationAccession", annotation_accession_col)
        my_naive_themisto_PCN_df = my_naive_themisto_PCN_df.with_columns([AnnotationAccession_series])

        ## and now add the columns with the naive themisto PCN estimate subset data.
        my_naive_themisto_PCN_df = my_naive_themisto_PCN_df.join(
            my_naive_themisto_PCN_estimates_subset_df, on="SeqID", how="left", coalesce=True).with_columns(
                ## set missing values in the InitialReadCount column to 0.
                pl.col("InitialReadCount").fill_null(strategy="zero")).select(
                ## super annoying, not sure why the AnnotationAccession_right columns
                ## are kept. This select command removes these redundancies.
                ["ThemistoID", "AnnotationAccession", "SeqID", "SeqType",
                 "replicon_length", "InitialReadCount"])
        
        ## make a dictionary mapping reads to Themisto replicon IDs.
        genome_dir = os.path.join(multiread_alignment_dir, genome_ID)
        multiread_mapping_dict = parse_read_alignments(genome_dir)
        
        ## initialize the data structures for PIRA.
        MatchMatrix, PIRAGenomeDataFrame = initializePIRA(
            multiread_mapping_dict, themisto_ID_to_seq_metadata_dict, my_naive_themisto_PCN_df)
        
        ## now run PIRA for this genome.
        PIRA_PCN_estimate_vector = run_PIRA(MatchMatrix, PIRAGenomeDataFrame)
        print(f"PIRA PCN estimate vector for genome {genome_ID} is: {PIRA_PCN_estimate_vector}")
        print("*****************************************************************************************")

        ## now add the PIRA estimates as a column to the PIRAGenomeDataFrame.
        ## First convert the NumPy array to a Polars Series
        PIRA_PCN_estimate_series = pl.Series("PIRA_CopyNumberEstimate", PIRA_PCN_estimate_vector)
        ## Then add the Polars Series of PIRA estimates  to the DataFrame with initial data
        my_PIRA_PCN_estimate_DataFrame = PIRAGenomeDataFrame.with_columns([PIRA_PCN_estimate_series])

        ## now concatenate the DataFrame for this genome to the big DataFrame for all genomes.
        all_PIRA_estimates_DataFrame = pl.concat(
            [all_PIRA_estimates_DataFrame,
             my_PIRA_PCN_estimate_DataFrame])

    ## arrange the columns of all_PIRA_estimates_DataFrame in a nice fashion.
    all_PIRA_estimates_DataFrame = all_PIRA_estimates_DataFrame.select(
        pl.col("AnnotationAccession", "SeqID", "SeqType",
               "ThemistoID", "replicon_length", "InitialReadCount", "AdditionalReadCount", "ReadCount",
               "SequencingCoverage", "LongestRepliconCoverage", "InitialCopyNumberEstimate",
               "PIRA_CopyNumberEstimate"))
    ## now save all_PIRA_estimates_DataFrame to disk.
    all_PIRA_estimates_DataFrame.write_csv(PIRA_PCN_csv_file)
    return




def get_RefSeqIDs_with_paired_end_reads(RunID_table_csv, SRA_data_dir):

    paired_end_1_pattern = f"{SRA_data_dir}/*_1.fastq"
    paired_end_2_pattern = f"{SRA_data_dir}/*_2.fastq"

    RunID_1_set = { basename(x).split("_1.fastq")[0] for x in glob.glob(paired_end_1_pattern) }
    RunID_2_set = { basename(x).split("_2.fastq")[0] for x in glob.glob(paired_end_2_pattern) }
    ## now take the intersection.
    paired_end_RunIDs = sorted(list(RunID_1_set & RunID_2_set))

    paired_end_RunID_table_df = (
        pl.read_csv(RunID_table_csv)
        .filter(pl.col("Run_ID").is_in(paired_end_RunIDs))
    )

    paired_end_RefSeqID_list = paired_end_RunID_table_df["RefSeq_ID"].to_list()
    return paired_end_RefSeqID_list


def choose_low_PCN_benchmark_genomes(PIRA_PCN_csv_file, PIRA_low_PCN_benchmark_csv_file,
                                     RunID_table_csv, SRA_data_dir):
    ## By default choose 100 genomes for the benchmarking against minimap2 PCN estimates.
    GENOME_SAMPLE_SIZE = 100
    
    ## ReadCount threshold for PCN estimates by pseudoalignment.
    ## This is only used for filtering the final set of PCN estimates by pseudoalignment.
    MIN_READ_COUNT = 10000 

    PCN_THRESHOLD = 0.8 ## threshold to define low PCN < 1 (to avoid selecting plasmids with PCN = 0.98).
    
    ## get the PIRA estimates from disk. 
    all_PIRA_estimates_DataFrame = pl.read_csv(PIRA_PCN_csv_file)

    ## Get genomes where all replicons have ReadCount > MIN_READ_COUNT,
    filtered_PIRA_estimates_DataFrame = (
        all_PIRA_estimates_DataFrame
        .with_columns(
            (pl.col("ReadCount") > MIN_READ_COUNT).alias("SufficientReadCount"))
        .group_by("AnnotationAccession")
        .agg([
            (pl.col("SufficientReadCount").all()).alias("All_Replicons_Have_SufficientReadCount")
        ])
        .filter(pl.col("All_Replicons_Have_SufficientReadCount"))
        .join(all_PIRA_estimates_DataFrame, on="AnnotationAccession", how="inner")
    ).filter(
        ## then filter for replicons with estimated PCN < PCN_THRESHOLD.
        pl.col("PIRA_CopyNumberEstimate") < PCN_THRESHOLD)

    ## get the unique annotation accessions in this dataframe
    annotation_accession_list = list(set(filtered_PIRA_estimates_DataFrame.get_column("AnnotationAccession").to_list()))

    ## only sample annotation accesions that have paired-end read data.
    paired_end_annotation_accession_list = list()
    paired_end_RefSeqID_list = get_RefSeqIDs_with_paired_end_reads(RunID_table_csv, SRA_data_dir)
    for annotation_accession in annotation_accession_list:
        my_parts = annotation_accession.split("_", 2)
        my_refseq_id = my_parts[0] + "_" + my_parts[1]
        if my_refseq_id in paired_end_RefSeqID_list:
            paired_end_annotation_accession_list.append(annotation_accession)

    ## shuffle the list
    random.shuffle(paired_end_annotation_accession_list)
    ## make sure the sample size is less than or equal to the length of the list
    assert GENOME_SAMPLE_SIZE <= len(paired_end_annotation_accession_list)
    ## and pick the first GENOME_SAMPLE_SIZE annotation accessions in the shuffled list.
    selected_annotation_accessions = paired_end_annotation_accession_list[0:GENOME_SAMPLE_SIZE]

    ## now subset the original dataframe based on these annotation accessions.
    random_subset_of_PIRA_estimates_DataFrame = all_PIRA_estimates_DataFrame.filter(
        pl.col("AnnotationAccession").is_in(selected_annotation_accessions))
    
    ## now save to disk.
    random_subset_of_PIRA_estimates_DataFrame.write_csv(PIRA_low_PCN_benchmark_csv_file)
    return


def benchmark_PCN_estimates_with_minimap2_alignments(
        PIRA_low_PCN_benchmark_csv_file, benchmark_alignment_dir,
        themisto_replicon_ref_dir, minimap2_benchmark_PIRA_PCN_csv_file):
    ## POTENTIAL TODO: refactor code as needed for here and in the function
    ## run_PIRA_on_all_genomes() as needed to reuse code and avoid duplication.
    
    ## import PIRA PCN estimates for the benchmarking genomes.
    PIRA_low_PCN_estimate_benchmark_df = pl.read_csv(PIRA_low_PCN_benchmark_csv_file)

    ## get the unique annotation accessions for the benchmark genomes.
    benchmark_genome_IDs = list(set(PIRA_low_PCN_estimate_benchmark_df.get_column("AnnotationAccession").to_list()))

    ## Make an empty Polars DataFrame to contain all the results.
    all_PIRA_estimates_DataFrame = pl.DataFrame()
    
    ## iterate over the benchmark genomes to populate all_PIRA_estimates_DataFrame.
    for my_genome_ID in benchmark_genome_IDs:

        my_genome_dir = os.path.join(benchmark_alignment_dir, my_genome_ID)
        assert os.path.isdir(my_genome_dir) ## make sure this directory exists.
        
        ## make the dictionary mapping reads to the multiset of themisto replicon IDs.
        read_mapping_dict = parse_read_alignments(my_genome_dir)

        themisto_ID_to_seq_metadata_dict = map_themisto_IDs_to_replicon_metadata(
            themisto_replicon_ref_dir, my_genome_ID)

        ## we need to run PIRA using ONLY the minimap2 alignment results.
        ## to do so, we need to pass in a data structure with some sensible default values.
        my_PCN_data_df = (
            PIRA_low_PCN_estimate_benchmark_df
            .filter(pl.col("AnnotationAccession") == my_genome_ID)
            ## only keep the metadata for this genome
            .select(pl.col("AnnotationAccession", "SeqID", "SeqType", "ThemistoID", "replicon_length"))
            ## and recreate the InitialReadCount column, set to zero.
            .with_columns(pl.lit(0).alias("InitialReadCount"))
        )
    
        ## now initialize the data structures for PIRA.
        MatchMatrix, PIRAGenomeDataFrame = initializePIRA(
            read_mapping_dict, themisto_ID_to_seq_metadata_dict, my_PCN_data_df)

        ## now run PIRA for this genome.
        PIRA_PCN_estimate_vector = run_PIRA(MatchMatrix, PIRAGenomeDataFrame)
        print(f"PIRA PCN estimate vector for genome {my_genome_ID} is: {PIRA_PCN_estimate_vector}")
        print("*****************************************************************************************")

        ## now add the PIRA estimates as a column to the PIRA_estimates_DataFrame.
        ## First convert the NumPy array to a Polars Series
        PIRA_PCN_estimate_series = pl.Series("PIRA_CopyNumberEstimate", PIRA_PCN_estimate_vector)
        ## Then add the Polars Series of PIRA estimates  to the DataFrame with initial data
        my_PIRA_PCN_estimate_DataFrame = PIRAGenomeDataFrame.with_columns([PIRA_PCN_estimate_series])

        ## now concatenate the DataFrame for this genome to the big DataFrame for all genomes.
        all_PIRA_estimates_DataFrame = pl.concat(
            [all_PIRA_estimates_DataFrame,
             my_PIRA_PCN_estimate_DataFrame])

    ## arrange the columns of all_PIRA_estimates_DataFrame in a nice fashion.
    all_PIRA_estimates_DataFrame = all_PIRA_estimates_DataFrame.select(
        pl.col(
            "AnnotationAccession", "SeqID", "SeqType",
            "ThemistoID", "replicon_length", "InitialReadCount", "AdditionalReadCount", "ReadCount",
            "SequencingCoverage", "LongestRepliconCoverage", "InitialCopyNumberEstimate",
            "PIRA_CopyNumberEstimate"))
    ## now save all_PIRA_estimates_DataFrame to disk.
    all_PIRA_estimates_DataFrame.write_csv(minimap2_benchmark_PIRA_PCN_csv_file)
    return


def benchmark_low_PCN_genomes_with_breseq(
        PIRA_low_PCN_benchmark_csv_file, RunID_table_csv,
        reference_genome_dir, SRA_data_dir, breseq_benchmark_results_dir):

    ## make the directory for breseq results  if it does not exist.
    if not exists(breseq_benchmark_results_dir):
        os.makedirs(breseq_benchmark_results_dir)

    ## get the AnnotationAccessions of interest for benchmarking.
    benchmark_genomes_df = pl.read_csv(PIRA_low_PCN_benchmark_csv_file)

    ## get the unique RefSeq_IDs in this dataframe
    AnnotationAccession_list = list(set(benchmark_genomes_df.get_column("AnnotationAccession").to_list()))
    ## the following is a hack to split AnnotationAccession on the second occurrence of an "_" to get the RefSeq_ID.
    AnnotationAccession_to_RefSeq_ID_dict = {x : "_".join(x.split("_", 2)[:2]) for x in AnnotationAccession_list}

    ## we will get the fastq files for each genome from this DataFrame.
    benchmark_RunID_table_df = pl.read_csv(RunID_table_csv).filter(
        pl.col("RefSeq_ID").is_in(AnnotationAccession_to_RefSeq_ID_dict.values()))
    
    for annotation_accession in AnnotationAccession_list:
        
        ## make a subdirectory for the breseq output.
        my_breseq_outdir = os.path.join(breseq_benchmark_results_dir, annotation_accession)
        if not exists(my_breseq_outdir):
            os.makedirs(my_breseq_outdir)

        ## IMPORTANT: breseq needs an unzipped version of this file. let's keep the original though (use -c option).
        ref_genome_gbk_gz_file = annotation_accession + "_genomic.gbff.gz"
        gz_reference_genome_path = os.path.join(reference_genome_dir, ref_genome_gbk_gz_file)
        
        ref_genome_gbk_file = annotation_accession + "_genomic.gbff"
        reference_genome_path = os.path.join(reference_genome_dir, ref_genome_gbk_file)

        gunzip_optionc_string = f"gunzip -c {gz_reference_genome_path} > {reference_genome_path}"
        subprocess.run(gunzip_optionc_string, shell=True)
        
        my_RefSeq_ID = AnnotationAccession_to_RefSeq_ID_dict[annotation_accession]
        my_RunID_df = benchmark_RunID_table_df.filter(pl.col("RefSeq_ID") == my_RefSeq_ID)
        my_RunID_list = my_RunID_df.get_column("Run_ID").to_list()

        ## Now use this list of Run IDs to pattern match for *.fastq files for this genome.
        my_read_data_pathlist = list()

        for Run_ID in my_RunID_list:
            SRA_file_pattern = f"{SRA_data_dir}/{Run_ID}*.fastq"
            matched_fastq_list = sorted(glob.glob(SRA_file_pattern))
            my_read_data_pathlist += matched_fastq_list
        my_read_data_pathlist.sort() ## sort the read data paths for this genome.

        ## now run breseq on this genome, with 8 cores.
        ncores = str(8)
        breseq_args = ["breseq", "-j", ncores, "-o", my_breseq_outdir, "-r", reference_genome_path] + my_read_data_pathlist
        breseq_string = " ".join(breseq_args)
        ## a few genomes can take more than 5 hours to run,
        ## and two genomes take more than 24h.
        sbatch_string = f"sbatch -p scavenger -t 48:00:00 --mem=16G --cpus-per-task={ncores} --wrap=\"" + breseq_string + "\""
        print(sbatch_string)
        if sys.platform == "linux": ## assume that we are running on DCC
            subprocess.run(sbatch_string, shell=True)
        else:
            print("WARNING: will only run breseq for this stage on DCC. exiting without running breseq.")
    return


def parse_breseq_results(breseq_outdir, results_csv_path):
    with open(results_csv_path, "w") as csv_fh:
        ## print a header string.
        print("AnnotationAccession,TrimmedSeqID,mean_coverage", file=csv_fh)
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
                ## breseq trims version numbers for SeqIDs. For example,
                ## "NZ_CP072604.1" is reported as "NZ_CP072604".
                TrimmedSeqID = row[2]
                mean_coverage = row[4]
                csv_string = ",".join([genome_accession, TrimmedSeqID, mean_coverage])
                print(csv_string, file=csv_fh)
    return


def run_pipeline_stage(stagenum, stage_complete_file, final_message, stage_function, *stage_function_args):
    if exists(stage_complete_file):
        print(f"{stage_complete_file} exists on disk-- skipping stage {stagenum}.")
    else:
        stage_start_time = time.time() ## Record the start time
        ## run the passed in function with its arguments for this stage.
        stage_function(*stage_function_args)
        stage_end_time = time.time() ## Record the end time.
        stage_execution_time = stage_end_time - stage_start_time
        StageTimeMessage = f"Stage {stagenum} execution time: {stage_execution_time} seconds\n"
        print(StageTimeMessage)
        logging.info(StageTimeMessage)
        with open(stage_complete_file, 'w') as stage_complete_log:
            stage_complete_log.write(StageTimeMessage)
            stage_complete_log.write(final_message)
        quit() ## exit the main program after finishing each new stage.
    return


def configure_logging(log_file):
    try:
        logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()  # Keep console output too
        ]
    )
    except IOError:
        ## If we can't write to the log file, just log to console
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler()
            ]
        )
        print("WARNING: Could not write to log file. Logging to console only.")

    logging.info("Starting the pipeline...")
    return


def create_test_subset(prokaryotes_with_plasmids_file, test_genome_count):
    """Create a subset of genomes for testing"""
    original_file = "../results/complete-prokaryotes-with-plasmids.txt"
    test_file = "../results/test-prokaryotes-with-plasmids.txt"
    
    logging.info(f"Creating test subset with {test_genome_count} genomes")
    
    ## Read the original file
    with open(original_file, "r") as f:
        lines = f.readlines()
        header = lines[0]
        data_lines = lines[1:]
        
        ## Take the first test_genome_count genomes or all if fewer
        subset_count = min(test_genome_count, len(data_lines))
        subset_lines = [header] + data_lines[:subset_count]
    
    ## Write the subset file
    with open(test_file, "w") as f:
        f.writelines(subset_lines)
        
    logging.info(f"Created test subset with {len(subset_lines)-1} genomes")
    
    return test_file


def test_pysradb_functionality():
    """Test if pysradb is working correctly with a known SRA ID."""
    test_sra_id = "SRS7822362"  
    logging.info(f"Testing pysradb with SRA ID: {test_sra_id}")
    
    try:
        # Run the command directly
        cmd = f"pysradb metadata {test_sra_id}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            logging.info(f"pysradb test successful. Output:\n{result.stdout}")
            
            # Try to parse the output
            lines = result.stdout.strip().split('\n')
            if len(lines) > 1:
                logging.info(f"Found {len(lines)-1} potential Run IDs")
                for line in lines[1:]:
                    logging.info(f"Potential Run ID line: {line}")
            else:
                logging.warning("No data rows in pysradb output")
        else:
            logging.error(f"pysradb test failed with error: {result.stderr}")
    except Exception as e:
        logging.error(f"Error testing pysradb: {e}")
    return


################################################################################
################################################################################
## THIS BLOCK OF ASYNC CODE FOR THEMISTO BUILD IS BROKEN-- DON'T USE FOR NOW
async def run_async_command_with_retries(command_string, tempdir=None, max_retries=5, timeout=20):
    """Run a shell command with retries, handling potential hangs."""
    retries = 0
    while retries < max_retries:
        process = await asyncio.create_subprocess_shell(
            command_string,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE
        )

        try:
            stdout, stderr = await asyncio.wait_for(process.communicate(), timeout=timeout)
        except asyncio.TimeoutError:
            process.kill()
            await process.wait()
            stdout, stderr = b"", b"Timeout occurred"

        if process.returncode == 0:
            logging.info("Command succeeded:", stdout.decode())
            return True
        else:
            logging.info(f"*********COMMAND\n{command_string}\nFAILED (attempt {retries + 1}):", stderr.decode())
            if tempdir is not None:
                logging.info(f"Removing {tempdir}")
                subprocess.run(f"rm -rf {tempdir}", shell=True)
                logging.info(f"Remaking {tempdir} before restarting")
                os.makedirs(tempdir)
            retries += 1
            await asyncio.sleep(0.1)
    
    logging.info(f"Command {command_string} failed after maximum retries.")
    return False


async def make_themisto_index_for_genome(genome_id, themisto_ref_dir, themisto_index_dir):
    """build a Themisto index for a single genome"""
    ref_fasta_dir = os.path.join(themisto_ref_dir, genome_id)
    if not os.path.isdir(ref_fasta_dir):
        return
    
    index_input_filelist = os.path.join(ref_fasta_dir, genome_id + ".txt")
    genome_index_dir = os.path.join(themisto_index_dir, genome_id)
    os.makedirs(genome_index_dir, exist_ok=True)
    
    index_prefix = os.path.join(genome_index_dir, genome_id)
    tempdir = os.path.join(genome_index_dir, "temp")
    os.makedirs(tempdir, exist_ok=True)
    
    themisto_build_args = [
        "themisto", "build", "-k", "31", "-i", index_input_filelist,
        "--index-prefix", index_prefix, "--temp-dir", tempdir,
        "--mem-gigas", "6", "--n-threads", "6", "--file-colors"
    ]
    themisto_build_string = " ".join(themisto_build_args)
    print(themisto_build_string)
    return await run_async_command_with_retries(themisto_build_string, tempdir)


async def make_themisto_indices_in_parallel(themisto_ref_dir, themisto_index_dir, max_concurrent=3):
    os.makedirs(themisto_index_dir, exist_ok=True)

    semaphore = asyncio.Semaphore(max_concurrent)

    async def make_themisto_index_with_semaphore(genome_id):
        async with semaphore:
            return await make_themisto_index_for_genome(genome_id, themisto_ref_dir, themisto_index_dir)

    tasks = [
        make_themisto_index_with_semaphore(genome_id)
        for genome_id in os.listdir(themisto_ref_dir)
        if os.path.isdir(os.path.join(themisto_ref_dir, genome_id))
    ]
    results = await asyncio.gather(*tasks, return_exceptions=True)

    success_count = sum(1 for r in results if r is True)
    logging.info(f"Successfully made indices for {success_count}/{len(results)} genomes")
    return success_count


def async_make_themisto_indices(themisto_ref_dir, themisto_index_dir):
    """Create Themisto indices for all genomes in the themisto reference directory."""
    asyncio.run(make_themisto_indices_in_parallel(themisto_ref_dir, themisto_index_dir))
    return
################################################################################
################################################################################
