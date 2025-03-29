#!/usr/bin/env python

"""
PCN_pipeline.py by Rohan Maddamsetti, Maggie Wilson, and Irida Shyti.

For this pipeline to work, ncbi datasets, pysradb, kallisto, themisto, minimap2, and breseq 0.39+ must be in the $PATH.
On the Duke Compute Cluster (DCC), run the following to get these programs into the path:
conda activate PCNdb-env

IMPORTANT: this script assumes it is being run on DCC if sys.platform == "linux".
This means that users on a linux machine will need to modify a couple functions if they
are running this code locally, and cannot use slurm to submit many jobs in parallel.

"""

## Import all code from PCN_library into this namespace.
from PCN_library import *


## Test mode configuration
TEST_MODE = False  ## Set to True for testing
TEST_GENOME_COUNT = 1000  ## Number of genomes to process
TEST_DOWNLOAD_LIMIT = 50  ## Increase from 10 to 50 for better testing

################################################################################
## Main pipeline code.
################################################################################

def run_PCN_pipeline():
    
    ## Configure logging
    log_dir = "../results"  # Use the results directory which is already mounted
    log_file = f"{log_dir}/pipeline_test.log" if TEST_MODE else f"{log_dir}/pipeline.log"
    configure_logging(log_file)

    if TEST_MODE:
        initialize_test_mode()
    else:
        logging.info("Running in PRODUCTION MODE")
    
        
    ## define input and output files used in the pipeline.
    if TEST_MODE:
        prokaryotes_with_plasmids_file = "../results/test-prokaryotes-with-plasmids.txt"
        RunID_table_csv = "../results/test-RunID_table.csv"
        reference_genome_dir = "../data/test-NCBI-reference-genomes/"
        SRA_data_dir = "../data/test-SRA/"
        # Create necessary directories
        os.makedirs(reference_genome_dir, exist_ok=True)
        os.makedirs(SRA_data_dir, exist_ok=True)
        stage1_complete_file = "../results/test-stage1.done"
        stage2_complete_file = "../results/test-stage2.done"
        stage3_complete_file = "../results/test-stage3.done"
    else:
        prokaryotes_with_plasmids_file = "../results/complete-prokaryotes-with-plasmids.txt"
        RunID_table_csv = "../results/RunID_table.csv"
        reference_genome_dir = "../data/NCBI-reference-genomes/"
        SRA_data_dir = "../data/SRA/"
        os.makedirs(SRA_data_dir, exist_ok=True)
        stage1_complete_file = "../results/stage1.done"
        stage2_complete_file = "../results/stage2.done"
        stage3_complete_file = "../results/stage3.done"

    
    replicon_length_csv_file = "../results/NCBI-replicon_lengths.csv"

    ## directories for themisto inputs and outputs.
    themisto_replicon_ref_dir = "../results/themisto_replicon_references/"
    themisto_replicon_index_dir = "../results/themisto_replicon_indices/"
    themisto_pseudoalignment_dir = "../results/themisto_replicon_pseudoalignments/"
    themisto_results_csvfile_path = "../results/themisto-replicon-read-counts.csv"

    ## this file contains estimates that throw out multireplicon reads.
    naive_themisto_PCN_csv_file = "../results/naive-themisto-PCN-estimates.csv"
    ## this file contains estimates that equally apportion multireplicon reads
    ## to the relevant plasmids and chromosomes.
    simple_themisto_PCN_csv_file = "../results/simple-themisto-PCN-estimates.csv"

    gbk_annotation_file = "../results/gbk-annotation-table.csv"

    ## directory for filtered multireads.
    multiread_data_dir = "../data/filtered_multireads/"

    ## directory for multiread alignments constructed with minimap2.
    multiread_alignment_dir = "../results/multiread_alignments/"

    ## this file contains PIRA estimates for the genomes that have multireads called by themisto.
    PIRA_PCN_csv_file = "../results/PIRA-PCN-estimates.csv"

    ## this file contains a random subset of PIRA estimates for 100 genomes with at least one plasmid with
    ## PCN < 1, and ReadCount > MIN_READ_COUNT.
    PIRA_low_PCN_benchmark_csv_file = "../results/PIRA-low-PCN-benchmark-estimates.csv"

    ## directories for replicon-level copy number estimation with kallisto.
    kallisto_replicon_ref_dir = "../results/kallisto_replicon_references/"
    kallisto_replicon_index_dir = "../results/kallisto_replicon_indices/"
    kallisto_replicon_quant_results_dir = "../results/kallisto_replicon_quant/"
    kallisto_replicon_copy_number_csv_file = "../results/kallisto-replicon_copy_numbers.csv"
    
    ## directory for read alignments constructed with minimap2 for PCN estimate benchmarking.
    benchmark_alignment_dir = "../results/PIRA_benchmark_alignments/"

    minimap2_benchmark_PIRA_PCN_csv_file = "../results/minimap2-PIRA-low-PCN-benchmark-estimates.csv"

    ## directory for breseq results for PCN estimate benchmarking.
    breseq_benchmark_results_dir = "../results/breseq_benchmark_results/"
    ## summary of breseq coverage data.
    breseq_benchmark_summary_file = "../results/breseq-low-PCN-benchmark-estimates.csv"
    
    
    #####################################################################################
    ## Stage 1: get SRA IDs and Run IDs for all RefSeq bacterial genomes with chromosomes and plasmids.
    logging.info("Stage 1: getting SRA IDs and Run IDs for RefSeq bacterial genomes with chromosomes and plasmids.")

    if exists(RunID_table_csv):
        Stage1DoneMessage = f"{RunID_table_csv} exists on disk-- skipping stage 1."
        print(Stage1DoneMessage)
        logging.info(Stage1DoneMessage)
    else:
        RunID_table_start_time = time.time()  # Record the start time
        
        if TEST_MODE:
            logging.info(f"Test mode: limiting API calls to {TEST_DOWNLOAD_LIMIT} genomes")
            # Create a limited version of the input file
            with open(prokaryotes_with_plasmids_file, "r") as f:
                lines = f.readlines()
                header = lines[0]
                data_lines = lines[1:TEST_DOWNLOAD_LIMIT+1]
            
            limited_file = prokaryotes_with_plasmids_file + ".limited"
            with open(limited_file, "w") as f:
                f.write(header)
                f.writelines(data_lines)
            
            # Use the limited file for API calls
            create_RefSeq_SRA_RunID_table(limited_file, RunID_table_csv)
        else:
            create_RefSeq_SRA_RunID_table(prokaryotes_with_plasmids_file, RunID_table_csv)
            
        RunID_table_end_time = time.time()  # Record the end time
        RunID_table_execution_time = RunID_table_end_time - RunID_table_start_time
        Stage1TimeMessage = f"Stage 1 execution time: {RunID_table_execution_time} seconds"
        print(Stage1TimeMessage)
        logging.info(Stage1TimeMessage)
        
        # Display the contents of the RunID table for verification
        logging.info("Contents of the RunID table:")
        with open(RunID_table_csv, "r") as f:
            table_contents = f.read()
            logging.info(table_contents)
            
        with open(stage1_complete_file, "w") as f:
            f.write(f"Stage 1 completed in {RunID_table_execution_time:.1f} seconds\n")
            quit()
        
        if TEST_MODE:
            logging.info("Test mode: Stage 1 completed successfully")

    #####################################################################################
    ## Stage 2: download reference genomes for each of the complete bacterial genomes containing plasmids,
    ## for which we can download Illumina reads from the NCBI Short Read Archive.
    ## first, make a dictionary from RefSeq accessions to ftp paths using the
    ## prokaryotes-with-plasmids.txt file.
    ## NOTE: as of March 28 2024, 4849 reference genomes should be downloaded.

    logging.info("Stage 2: downloading gzipped genbank reference genomes for complete genomes containing plasmids.")
    if exists(stage2_complete_file):
        print(f"{stage2_complete_file} exists on disk-- skipping stage 2.")
        logging.info(f"{stage2_complete_file} exists on disk-- skipping stage 2.")
    else:
        refgenome_download_start_time = time.time()  ## Record the start time
        refseq_accession_to_ftp_path_dict = create_refseq_accession_to_ftp_path_dict(prokaryotes_with_plasmids_file)
        
        # Log the FTP paths for debugging
        logging.info(f"Found {len(refseq_accession_to_ftp_path_dict)} FTP paths")
        
        ## now download the reference genomes.
        fetch_reference_genomes(RunID_table_csv, refseq_accession_to_ftp_path_dict, reference_genome_dir)
        refgenome_download_end_time = time.time()  ## Record the end time
        refgenome_download_execution_time = refgenome_download_end_time - refgenome_download_start_time
        Stage2TimeMessage = f"Stage 2 (reference genome download) execution time: {refgenome_download_execution_time} seconds\n"
        logging.info(Stage2TimeMessage)
        
        with open(stage2_complete_file, "w") as stage2_complete_log:
            stage2_complete_log.write(Stage2TimeMessage)
            stage2_complete_log.write("reference genomes downloaded successfully.\n")
            quit()
            
        if TEST_MODE:
            logging.info("Test mode: Stage 2 completed successfully")

    
    #####################################################################################
    ## Stage 3: download Illumina reads for the genomes from the NCBI Short Read Archive (SRA).
    logging.info("Stage 3: downloading Illumina reads from the NCBI Short Read Archive (SRA).")
    if exists(stage3_complete_file):
        print(f"{stage3_complete_file} exists on disk-- skipping stage 3.")
        logging.info(f"{stage3_complete_file} exists on disk-- skipping stage 3.")
    else:
        SRA_download_start_time = time.time()
        
        try:
            ## Get Run IDs from disk.
            Run_IDs = get_Run_IDs_from_RunID_table(RunID_table_csv)
            logging.info(f"Found {len(Run_IDs)} Run IDs to download")
            
            if TEST_MODE and len(Run_IDs) > TEST_DOWNLOAD_LIMIT:
                logging.info(f"Test mode: limiting downloads to {TEST_DOWNLOAD_LIMIT} Run IDs")
                Run_IDs = Run_IDs[:TEST_DOWNLOAD_LIMIT]

            ## Run parallel prefetch pipeline with reduced concurrency and retries
            prefetch_success = asyncio.run(prefetch_fastq_reads_parallel(
                SRA_data_dir, 
                Run_IDs, 
                max_concurrent=10,  ## Reduce concurrency to avoid I/O issues
                max_retries=3      ## Add retries
            ))

            if prefetch_success:
                ## then validate each sra file to check for data corruption.
                validation_success = asyncio.run(validate_sra_download_parallel(
                    SRA_data_dir, 
                    Run_IDs, 
                    max_concurrent=10,
                ))
                
                if validation_success:
                    ## Then run fasterq-dump in parallel
                    fasterq_dump_success = asyncio.run(unpack_fastq_reads_parallel(
                        SRA_data_dir, 
                        Run_IDs, 
                        max_concurrent=10
                    ))

                    if fasterq_dump_success:
                        ## then check to see that all fastq files exist on disk
                        if all_fastq_data_exist(Run_IDs, SRA_data_dir):
                            SRA_download_end_time = time.time()
                            SRA_download_execution_time = SRA_download_end_time - SRA_download_start_time
                            Stage3TimeMessage = f"Stage 3 downloads completed successfully in {SRA_download_execution_time:.1f} seconds\n"
                            logging.info(Stage3TimeMessage)
                            try:
                                with open(stage_3_complete_file, "w") as stage3_complete_log:
                                    stage3_complete_log.write(Stage3TimeMessage)
                                    stage3_complete_log.write("reference genomes downloaded successfully.\n")
                            except IOError as e:
                                logging.error(f"Could not write to {stage_3_complete_file}: {str(e)}")
                            if TEST_MODE:
                                logging.info("Test mode: Stage 3 completed successfully")
                        else:
                            logging.warning("Warning: All sra files downloaded and validated, but some fastq files may be missing")
                    else:
                        logging.error("Error: All files prefetched and validated but fasterq-dump failed in some cases.")
                else:
                    logging.error("Error: Some prefetch downloads failed validation.. possible data corruption.")
            else:
                logging.error("All prefetch of SRA data failed.")
                if TEST_MODE:
                    logging.info("Test mode: Stage 3 completed with download failures")
        except Exception as e:
            logging.error(f"Error in Stage 3: {str(e)}")
            if TEST_MODE:
                logging.info("Test mode: Stage 3 failed with errors")

        with open(stage3_complete_file, "w") as stage3_complete_log:
            SRA_download_end_time = time.time()
            SRA_final_execution_time = SRA_download_end_time - SRA_download_start_time
            Stage3TimeMessage = f"Stage 3 completed in {SRA_final_execution_time:.1f} seconds\n"
            stage3_complete_log.write(Stage3TimeMessage)
            stage3_complete_log.write("FASTQ data download from SRA complete.\n")
            quit()

        ## Exit after stage 3 in test mode
        if TEST_MODE:
            logging.info("Test mode: All 3 stages completed. Exiting.")
            return  ## Exit the main function

    #####################################################################################
    ## Stage 4: make gbk ecological annotation file.
    stage4_complete_file = "../results/stage4.done"
    stage4_final_message = "stage 4 (gbk ecological annotation) finished successfully.\n"
    run_pipeline_stage(4, stage4_complete_file, stage4_final_message,
                       make_gbk_annotation_table,
                       reference_genome_dir, gbk_annotation_file)
    
    #####################################################################################
    ## Stage 5: tabulate the length of all chromosomes and plasmids.
    stage5_complete_file = "../results/stage5.done"
    stage5_final_message = "Stage 5 (tabulating all replicon lengths) finished successfully.\n"
    run_pipeline_stage(5, stage5_complete_file, stage5_final_message,
                       tabulate_NCBI_replicon_lengths,
                       reference_genome_dir, replicon_length_csv_file)
    
    #####################################################################################
    ## Stage 6: Make FASTA input files for Themisto.
    ## Write out separate fasta files for each replicon in each genome, in a directory for each genome.
    ## Then, write out a text file that contains the paths to the FASTA files of the genomes, one file per line.
    ## See documentation here: https://github.com/algbio/themisto.
    stage6_complete_file = "../results/stage6.done"
    stage6_final_message = "Stage 6 (making fasta references for themisto) finished successfully.\n"
    run_pipeline_stage(6, stage6_complete_file, stage6_final_message,
                       make_NCBI_replicon_fasta_refs_for_themisto,
                       reference_genome_dir, themisto_replicon_ref_dir)
    
    #####################################################################################
    ## Stage 7: Build separate Themisto indices for each genome.
    stage7_complete_file = "../results/stage7.done"
    stage7_final_message = "Stage 7 (making indices for themisto) finished successfully.\n"
    run_pipeline_stage(7, stage7_complete_file, stage7_final_message,
                       make_themisto_indices,
                       themisto_replicon_ref_dir, themisto_replicon_index_dir)

    #####################################################################################
    ## Stage 8: Pseudoalign reads for each genome against each Themisto index.
    stage8_complete_file = "../results/stage8.done"
    stage8_final_message = "Stage 8 (themisto pseudoalignment) finished successfully.\n"
    run_pipeline_stage(8, stage8_complete_file, stage8_final_message,
                       run_themisto_pseudoalign,
                       RunID_table_csv, themisto_replicon_index_dir, SRA_data_dir, themisto_pseudoalignment_dir)
    
    #####################################################################################
    ## Stage 9: generate a large CSV file summarizing the themisto pseudoalignment read counts.
    stage9_complete_file = "../results/stage9.done"
    stage9_final_message = "Stage 9 (Themisto pseudoalignment summarization) finished successfully.\n"
    run_pipeline_stage(9, stage9_complete_file, stage9_final_message,
                       summarize_themisto_pseudoalignment_results,
                       themisto_replicon_ref_dir, themisto_pseudoalignment_dir, themisto_results_csvfile_path)
    
    #####################################################################################
    ## Stage 10: estimate plasmid copy numbers using the themisto read counts.
    stage10_complete_file = "../results/stage10.done"
    stage10_final_message = "Stage 10 (naive Themisto PCN estimates) finished successfully.\n"
    run_pipeline_stage(10, stage10_complete_file, stage10_final_message,
    ## Naive PCN calculation, ignoring multireplicon reads.
                       naive_themisto_PCN_estimation,
                       themisto_results_csvfile_path, replicon_length_csv_file, naive_themisto_PCN_csv_file)

    #####################################################################################
    ## Use Probabilistic Iterative Read Assignment (PIRA) to improve PCN estimates.
    #####################################################################################
    ## The Naive PCN calculation in Stage 10 generates the initial PCN vectors and saves them on disk.
    
    #####################################################################################
    ## Stage 11: filter fastq reads for multireads.
    stage11_complete_file = "../results/stage11.done"
    stage11_final_message = "Stage 11 (fastq read filtering) finished successfully.\n"
    run_pipeline_stage(11, stage11_complete_file, stage11_final_message,
                       filter_fastq_files_for_multireads,
                       multiread_data_dir, themisto_pseudoalignment_dir, SRA_data_dir)

    #####################################################################################
    ## Stage 12: make FASTA reference genomes with Themisto Replicon IDs for multiread mapping with minimap2.
    stage12_complete_file = "../results/stage12.done"
    stage12_final_message = "Stage 12 (making FASTA reference genomes for multiread alignment) finished successfully.\n"
    run_pipeline_stage(12, stage12_complete_file, stage12_final_message,
                       make_fasta_reference_genomes_for_minimap2,
                       themisto_replicon_ref_dir)

    #####################################################################################
    ## Stage 13: for each genome, align multireads to the replicons with minimap2.
    stage13_complete_file = "../results/stage13.done"
    stage13_final_message = "Stage 13 (aligning multireads with minimap2) finished successfully.\n"
    run_pipeline_stage(13, stage13_complete_file, stage13_final_message,
                       align_multireads_with_minimap2,
                       themisto_replicon_ref_dir, multiread_data_dir, multiread_alignment_dir)

    #####################################################################################
    ## Stage 14: Run PIRA.
    ## For each genome, parse minimap2 results to form the match matrix and refine the initial PCN guesses.
    stage14_complete_file = "../results/stage14.done"
    stage14_final_message = "Stage 14 (parsing multiread alignments and running PIRA) finished successfully.\n"
    run_pipeline_stage(14, stage14_complete_file, stage14_final_message,
                       run_PIRA_on_all_genomes,
                       multiread_alignment_dir, themisto_replicon_ref_dir, naive_themisto_PCN_csv_file, PIRA_PCN_csv_file)

    #####################################################################################
    ## In order to benchmark accuracy, speed, and memory usage, estimate PCN for a subset of 100 genomes
    ## that apparently contain low PCN plasmids (PCN < 0.8). 
    #####################################################################################
    ## Stage 15: choose a set of 100 random genomes that contain plasmids with PCN < 1
    ## and ReadCount > MIN_READ_COUNT.
    
    stage15_complete_file = "../results/stage15.done"
    stage15_final_message = "Stage 15 (choosing 100 random genomes with PCN < 1) finished successfully.\n"
    run_pipeline_stage(15, stage15_complete_file, stage15_final_message,
                       ## at random, choose 100 genomes containing a plasmid with PCN < 1, paired-end read data,
                       ## and ReadCount > 10000.
                       choose_low_PCN_benchmark_genomes,
                       PIRA_PCN_csv_file, PIRA_low_PCN_benchmark_csv_file, RunID_table_csv, SRA_data_dir)

    
    #####################################################################################
    ## Benchmark PIRA estimates against kallisto.
    #####################################################################################
    ## In order to benchmark accuracy, estimate PCN for a subset of 100 genomes
    ## that apparently contain low PCN plasmids (PCN < 0.8), using kallisto.
    #####################################################################################   
    ## Stage 16: Make replicon-level FASTA reference files for copy number estimation using kallisto.
    stage16_complete_file = "../results/stage16.done"
    stage16_final_message = "Stage 16 (FASTA reference sequences for kallisto) finished successfully.\n"
    run_pipeline_stage(16, stage16_complete_file, stage16_final_message,
                       make_NCBI_replicon_fasta_refs_for_kallisto,
                       reference_genome_dir, kallisto_replicon_ref_dir)
    
    #####################################################################################
    ## Stage 17: Make replicon-level kallisto index files for each genome.
    stage17_complete_file = "../results/stage17.done"
    stage17_final_message = "Stage 17 (kallisto index file construction) finished successfully.\n"
    run_pipeline_stage(17, stage17_complete_file, stage17_final_message,
    make_NCBI_kallisto_indices,
    kallisto_replicon_ref_dir, kallisto_replicon_index_dir)
    
    #####################################################################################
    ## Stage 18: run kallisto quant on all genome data, on replicon-level indices.
    ## NOTE: right now, this only processes paired-end fastq data-- single-end fastq data is ignored.
    stage18_complete_file = "../results/stage18.done"
    stage18_final_message = "Stage 18 (kallisto quant) finished successfully.\n"
    run_pipeline_stage(18, stage18_complete_file, stage18_final_message,
                       run_kallisto_quant,
                       RunID_table_csv, kallisto_replicon_index_dir, SRA_data_dir, kallisto_replicon_quant_results_dir)
    
    #####################################################################################
    ## Stage 19: make a table of the estimated copy number for all chromosomes and plasmids using kallisto
    stage19_complete_file = "../results/stage19.done"
    stage19_final_message = "Stage 19 (tabulating replicon copy numbers estimated by kallisto) finished successfully.\n"
    run_pipeline_stage(19, stage19_complete_file, stage19_final_message,
                       measure_kallisto_replicon_copy_numbers,
                       kallisto_replicon_quant_results_dir, kallisto_replicon_copy_number_csv_file)

    #####################################################################################
    ## Benchmark PIRA estimates against traditional alignment PCN estimation with minimap2 and breseq.
    #####################################################################################
    ## Stage 20: run minimap2 on the set of 100 genomes chosen in Stage 15.
    stage20_complete_file = "../results/stage20.done"
    stage20_final_message = "stage 20 (minimap2 full alignment) finished successfully.\n"
    run_pipeline_stage(20, stage20_complete_file, stage20_final_message,
                       align_reads_for_benchmark_genomes_with_minimap2,
                       PIRA_low_PCN_benchmark_csv_file, RunID_table_csv, themisto_replicon_ref_dir,
                       SRA_data_dir, benchmark_alignment_dir)
    
    #####################################################################################
    ## Stage 21: parse minimap2 results on the set of 100 genomes chosen for benchmarking.
    stage21_complete_file = "../results/stage21.done"
    stage21_final_message = "stage 21 (minimap2 results parsing) finished successfully.\n"
    run_pipeline_stage(21, stage21_complete_file, stage21_final_message,
                       benchmark_PCN_estimates_with_minimap2_alignments,
                       PIRA_low_PCN_benchmark_csv_file, benchmark_alignment_dir,
                       themisto_replicon_ref_dir, minimap2_benchmark_PIRA_PCN_csv_file)
    
    #####################################################################################
    ## Stage 22: run breseq on the set of 100 genomes chosen for benchmarking.
    stage22_complete_file = "../results/stage22.done"
    stage22_final_message = "stage 22 (running breseq) finished successfully.\n"
    run_pipeline_stage(22, stage22_complete_file, stage22_final_message,
                       benchmark_low_PCN_genomes_with_breseq,
                       PIRA_low_PCN_benchmark_csv_file, RunID_table_csv,
                       reference_genome_dir, SRA_data_dir, breseq_benchmark_results_dir)        

    #####################################################################################
    ## Stage 23: parse breseq results on the set of 100 genomes chosen for benchmarking.
    stage23_complete_file = "../results/stage23.done"
    stage23_final_message = "stage 23 (breseq results parsing) finished successfully.\n"
    run_pipeline_stage(23, stage23_complete_file, stage23_final_message,
                       parse_breseq_results,
                       breseq_benchmark_results_dir, breseq_benchmark_summary_file)
    return


if __name__ == "__main__":
    run_PCN_pipeline()

