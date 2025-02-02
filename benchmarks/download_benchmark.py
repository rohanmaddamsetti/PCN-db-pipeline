import time
import os
from os.path import basename, dirname, join
import random
from concurrent.futures import ThreadPoolExecutor, as_completed
import asyncio
import sys
import urllib.request

# Add the src directory to the Python path so we can import from PCN_pipeline
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from PCN_pipeline import (fetch_reference_genomes, create_RefSeq_SRA_RunID_table, 
                         create_refseq_accession_to_ftp_path_dict)

# Get the root directory of the project
ROOT_DIR = os.path.dirname(os.path.dirname(__file__))

# Create data directory if it doesn't exist
DATA_DIR = os.path.join(ROOT_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)

# Test directory and files
TEST_DIR = "test_downloads"
os.makedirs(TEST_DIR, exist_ok=True)

# Log file
LOG_FILE = "download_benchmark.log"

# Files
PROKARYOTES_FILE = os.path.join(DATA_DIR, "prokaryotes-with-plasmids.txt")
RUNID_TABLE = os.path.join(DATA_DIR, "RunID_table.csv")

def download_prokaryotes_file():
    """Download the prokaryotes file from NCBI"""
    if not os.path.exists(PROKARYOTES_FILE):
        print("Downloading prokaryotes file...")
        # Updated URL to the new NCBI path
        url = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
        try:
            urllib.request.urlretrieve(url, PROKARYOTES_FILE)
            print("Download complete!")
        except urllib.error.HTTPError as e:
            print(f"Error downloading file: {e}")
            print("Please check if the file exists at the URL or if NCBI's FTP structure has changed.")
            sys.exit(1)

def create_RunID_table():
    """Create RunID table from prokaryotes file"""
    print("\nCreating RunID table...")
    try:
        with open(PROKARYOTES_FILE, 'r') as f, open(RUNID_TABLE, 'w') as out:
            # Write header
            out.write("RefSeq_ID,SRA_ID,Run_ID\n")
            
            # Skip header line
            next(f)
            
            # Process each line
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 17:  # Make sure we have enough fields
                    assembly_acc = fields[17]  # Assembly Accession column
                    if assembly_acc.startswith('GCA_'):
                        # Convert GCA to GCF
                        refseq_id = 'GCF_' + assembly_acc[4:]
                        out.write(f"{refseq_id},,\n")  # SRA_ID and Run_ID will be filled later
                        print(f"Added {refseq_id} to RunID table")
    except Exception as e:
        print(f"Error creating RunID table: {e}")
        sys.exit(1)

def setup_benchmark():
    """Create necessary files and directories for benchmarking"""
    # Only download prokaryotes file if it doesn't exist
    if not os.path.exists(PROKARYOTES_FILE):
        print("Downloading prokaryotes file...")
        download_prokaryotes_file()
    else:
        print("Using existing prokaryotes file...")
    
    if not os.path.exists(RUNID_TABLE):
        print("\nCreating RunID table...")
        create_RunID_table()  # Use our simplified version instead
    
    # Create the FTP path dictionary
    print("\nCreating FTP path dictionary...")
    ftp_path_dict = create_refseq_accession_to_ftp_path_dict(PROKARYOTES_FILE)
    
    # Create test RunID table with actual data
    test_size = 6  # Number of genomes to test with
    all_refseq_ids = list(ftp_path_dict.keys())
    test_refseq_ids = random.sample(all_refseq_ids, test_size)
    
    # Create test FTP dict
    test_ftp_dict = {id: ftp_path_dict[id] for id in test_refseq_ids}
    
    # Create test RunID table
    with open(RUNID_TABLE, 'w') as out:
        out.write("RefSeq_ID,SRA_ID,Run_ID\n")
        for refseq_id in test_refseq_ids:
            out.write(f"{refseq_id},,\n")
    
    print(f"\nTest subset created with {test_size} genomes")
    print("Test RefSeq IDs:")
    for id in test_refseq_ids:
        print(f"- {id}")
    
    print("\nFTP paths for test RefSeq IDs:")
    for id in test_refseq_ids:
        print(f"- {id}: {test_ftp_dict[id]}")
    
    print("\nRunID table content:")
    with open(RUNID_TABLE, 'r') as f:
        print(f.read())
        
    return test_ftp_dict

def cleanup_test_files():
    """Remove all files in the test directory"""
    for file in os.listdir(TEST_DIR):
        try:
            os.remove(os.path.join(TEST_DIR, file))
        except Exception as e:
            print(f"Error deleting {file}: {e}")

def get_test_data():
    """Get test data from the RunID table"""
    print(f"Reading RunID table from: {RUNID_TABLE}")
    with open(RUNID_TABLE, 'r') as f:
        # Skip header
        next(f)
        # Get unique RefSeq IDs
        refseq_ids = set()
        for line in f:
            refseq_id = line.split(',')[0]
            refseq_ids.add(refseq_id)
    print(f"Found {len(refseq_ids)} RefSeq IDs")
    if len(refseq_ids) == 0:
        print("Warning: No RefSeq IDs found in the RunID table!")
    return list(refseq_ids)

def get_random_subset(data, size):
    """Get a random subset of test data"""
    return random.sample(data, min(size, len(data)))

# Current method (baseline)
def current_method(ftp_path_dict):
    start_time = time.time()
    print(f"\nFTP dictionary contains {len(ftp_path_dict)} entries:")
    for key, value in ftp_path_dict.items():
        print(f"- {key}: {value}")
        
    print("\nAttempting to fetch reference genomes...")
    try:
        fetch_reference_genomes(RUNID_TABLE, ftp_path_dict, TEST_DIR, LOG_FILE)
    except Exception as e:
        print(f"Error in fetch_reference_genomes: {e}")
        import traceback
        traceback.print_exc()
    return time.time() - start_time

def download_single_genome(ftp_path, reference_genome_dir, log_file):
    """Download a single genome and its MD5 file"""
    my_full_accession = basename(ftp_path)
    my_base_filename = my_full_accession + "_genomic.gbff.gz"
    
    # Files on the NCBI FTP site to download
    gbff_ftp_path = os.path.join(ftp_path, my_base_filename)
    md5_ftp_path = os.path.join(ftp_path, "md5checksums.txt")
    
    # Local paths
    gbff_gz_file = os.path.join(reference_genome_dir, my_base_filename)
    md5_file = os.path.join(reference_genome_dir, my_full_accession + "_md5checksums.txt")
    
    try:
        urllib.request.urlretrieve(gbff_ftp_path, filename=gbff_gz_file)
        urllib.request.urlretrieve(md5_ftp_path, filename=md5_file)
        print(f"{gbff_gz_file} SUCCEEDED.")
        return True
    except Exception as e:
        print(f"{gbff_gz_file} FAILED: {str(e)}")
        return False

# ThreadPoolExecutor method (properly parallelized)
def threadpool_download(ftp_path_dict, max_workers=3):
    start_time = time.time()
    print("\nAttempting parallel downloads with ThreadPoolExecutor...")
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for ftp_path in ftp_path_dict.values():
            future = executor.submit(download_single_genome, ftp_path, TEST_DIR, LOG_FILE)
            futures.append(future)
        
        for future in as_completed(futures):
            try:
                result = future.result()
            except Exception as e:
                print(f"Download failed: {e}")
    
    return time.time() - start_time

# Asyncio method (properly parallelized)
async def async_download_single(ftp_path):
    """Download a single genome asynchronously"""
    return await asyncio.to_thread(download_single_genome, ftp_path, TEST_DIR, LOG_FILE)

async def async_download(ftp_path_dict):
    start_time = time.time()
    print("\nAttempting parallel downloads with asyncio...")
    
    tasks = []
    for ftp_path in ftp_path_dict.values():
        task = asyncio.create_task(async_download_single(ftp_path))
        tasks.append(task)
    
    await asyncio.gather(*tasks)
    return time.time() - start_time

def run_tests():
    # Setup the benchmark environment
    print("Setting up benchmark...")
    ftp_path_dict = setup_benchmark()
    
    # Print RunID table content
    print("\nRunID table content:")
    with open(RUNID_TABLE, 'r') as f:
        print(f.read())

    print("\nRunning current method...")
    cleanup_test_files()  # Clean before each test
    current_time = current_method(ftp_path_dict)
    print(f"Current method took {current_time:.2f} seconds")

    print("\nRunning ThreadPoolExecutor method...")
    cleanup_test_files()  # Clean before each test
    threadpool_time = threadpool_download(ftp_path_dict)
    print(f"ThreadPoolExecutor method took {threadpool_time:.2f} seconds")

    print("\nRunning asyncio method...")
    cleanup_test_files()  # Clean before each test
    asyncio_time = asyncio.run(async_download(ftp_path_dict))
    print(f"Asyncio method took {asyncio_time:.2f} seconds")

    print("\nResults:")
    print(f"Current method: {current_time:.2f}s")
    print(f"ThreadPoolExecutor: {threadpool_time:.2f}s ({current_time/threadpool_time:.1f}x faster)")
    print(f"Asyncio: {asyncio_time:.2f}s ({current_time/asyncio_time:.1f}x faster)")

    # Check what was downloaded
    print("\nFiles in test directory:")
    if os.path.exists(TEST_DIR):
        files = os.listdir(TEST_DIR)
        if files:
            for file in files:
                print(f"- {file}")
        else:
            print("No files were downloaded!")

def cleanup_cached_files():
    """Remove all cached files to start fresh"""
    files_to_remove = [
        PROKARYOTES_FILE,  # prokaryotes file
        RUNID_TABLE,       # RunID table
        RUNID_TABLE + ".filtered",  # filtered RunID table
        RUNID_TABLE + ".test",      # test RunID table
        LOG_FILE,          # log file
    ]
    
    print("Cleaning up cached files...")
    for file in files_to_remove:
        try:
            if os.path.exists(file):
                os.remove(file)
                print(f"Removed: {file}")
        except Exception as e:
            print(f"Error removing {file}: {e}")
    
    # Clean test directory
    if os.path.exists(TEST_DIR):
        for file in os.listdir(TEST_DIR):
            try:
                os.remove(os.path.join(TEST_DIR, file))
                print(f"Removed: {os.path.join(TEST_DIR, file)}")
            except Exception as e:
                print(f"Error removing {os.path.join(TEST_DIR, file)}: {e}")

if __name__ == "__main__":
    # Clean up before running tests
    cleanup_cached_files()
    
    # Run the tests
    run_tests()