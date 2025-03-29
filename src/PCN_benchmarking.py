#!/usr/bin/env python

"""
PCN_benchmarking.py by Rohan Maddamsetti, Maggie Wilson, and Irida Shyti.

For this pipeline to work, ncbi datasets, pysradb, kallisto, themisto, minimap2, and breseq 0.39+ must be in the $PATH.
On the Duke Compute Cluster (DCC), run the following to get these programs into the path:
conda activate PCNdb-env

IMPORTANT: this script assumes it is being run on DCC if sys.platform == "linux".
This means that users on a linux machine will need to modify a couple functions if they
are running this code locally, and cannot use slurm to submit many jobs in parallel.

"""

## Import all code from PCN_library into this namespace.
from PCN_library import *


def main_benchmarking():

    ## Configure logging
    log_dir = "../results"  # Use the results directory which is already mounted
    log_file =  f"{log_dir}/benchmarking.log"
    configure_logging(log_file)

    ## directories for replicon-level copy number estimation with kallisto.
    kallisto_replicon_ref_dir = "../results/kallisto_replicon_references/"
    kallisto_replicon_index_dir = "../results/kallisto_replicon_indices/"
    kallisto_replicon_quant_results_dir = "../results/kallisto_replicon_quant/"
    kallisto_replicon_copy_number_csv_file = "../results/kallisto-replicon_copy_numbers.csv"

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

    ## directory for read alignments constructed with minimap2 for PCN estimate benchmarking.
    benchmark_alignment_dir = "../results/PIRA_benchmark_alignments/"

    minimap2_benchmark_PIRA_PCN_csv_file = "../results/minimap2-PIRA-low-PCN-benchmark-estimates.csv"

    ## directory for breseq results for PCN estimate benchmarking.
    breseq_benchmark_results_dir = "../results/breseq_benchmark_results/"
    ## summary of breseq coverage data.
    breseq_benchmark_summary_file = "../results/breseq-low-PCN-benchmark-estimates.csv"


    ## TODO: ONLY SELECT GENOMES WITH PAIRED-END READS FOR KALLISTO!
    print("Hello")
    quit()
    #####################################################################################
    ## In order to benchmark accuracy, speed, and memory usage, estimate PCN for a subset of 100 genomes
    ## that apparently contain low PCN plasmids (PCN < 0.8). 
    #####################################################################################
    ## Stage 15: choose a set of 100 random genomes that contain plasmids with PCN < 1
    ## and ReadCount > MIN_READ_COUNT.
    stage15_complete_file = "../results/stage15.done"
    stage15_final_message = "Stage 15 (choosing 100 random genomes with PCN < 1) finished successfully.\n"
    run_pipeline_stage(15, stage15_complete_file, stage15_final_message,
                       ## at random, choose 100 genomes containing a plasmid with PCN < 1, and ReadCount > 10000.
                       choose_low_PCN_benchmark_genomes,
                       PIRA_PCN_csv_file, PIRA_low_PCN_benchmark_csv_file)
    
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
    main_benchmarking()
