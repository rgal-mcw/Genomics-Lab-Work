# Title: Read Alignment
# Author: Ryan Gallagher (Broeckel Lab, 2024)
#
# Description: This script aligns the fastq output from the Prom to GRCh38
#              using minimap2. This output .bam file will then be sorted.


import os
import subprocess
import argparse
import logging
import sys

def setup_logging(log_file_path, log_level=logging.INFO):
    """Configures logging for the script."""
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path, mode='a'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)

def run_alignment(ref, fastq, output_name, threads=48, logger=None):

    log = logger or logging.getLogger(__name__) # Use passed logger or get default

    # Build the minimap2 command with samtools view
    minimap2_command = f"minimap2 -L -t {threads} -ax map-ont {ref} {fastq}"
    samtools_view_command = f"samtools view -bS -@ {threads}"
    combined_command = f"{minimap2_command} | {samtools_view_command} > {output_name}.bam"

    log.info("Running command: " + combined_command)
    subprocess.run(combined_command, shell=True)

    # Run Samtools sort
    sort_command = f"samtools sort -@ {threads} {output_name}.bam -o {output_name}.sorted.bam"
    log.info("Running command: " + sort_command)
    subprocess.run(sort_command, shell=True)

    # Run samtools index
    index_command = f"samtools index -@ {threads} {output_name}.sorted.bam"
    log.info("Running command: " + index_command)
    subprocess.run(index_command, shell=True)

    # Remove the unsorted bam file
    try:
        os.remove(f"{output_name}.bam")
        log.info(f"Removed intermediate file: {output_name}.bam") 
    except OSError as e:
        log.error(f"Error removing {output_name}.bam: {e}") 


def main():

    parser = argparse.ArgumentParser(description='Align reads using minimap2 and samtools')
    parser.add_argument('sample_name', type=str, help='Enter the sample name')
    parser.add_argument('to_dir', type=str, help='Enter the MAPLE directory where this should be stored')
    parser.add_argument('ref', type=str, help='Enter the reference genome')
    args = parser.parse_args()

    # --- Setup Logging (using existing args) ---
    # Construct log directory path based on to_dir
    log_dir = os.path.join(args.to_dir, 'logs')
    # Ensure log directory exists (create if it doesn't)
    os.makedirs(log_dir, exist_ok=True)
    # Define a specific log file for this script within that directory
    log_file = os.path.join(log_dir, 'alignment.log')
    # Setup logging with default level (INFO)
    # Change logging.INFO to logging.DEBUG if you want more verbose logs ("log anything")
    logger = setup_logging(log_file, log_level=logging.INFO)
    # --- End Logging Setup ---
    logger.info(f"Starting alignment script for sample: {args.sample_name}") 

    # Ensure paths are constructed safely
    # Ensure to_dir ends with a slash if needed, or use os.path.join consistently
    args.to_dir = args.to_dir if args.to_dir.endswith('/') else args.to_dir + '/'
    fastq = os.path.join(args.to_dir, 'raw_data', f'{args.sample_name}.fastq.gz')
    output_name = os.path.join(args.to_dir, args.sample_name)

    if not os.path.exists(fastq):
        logger.error(f"FASTQ file not found: {fastq}")
        sys.exit(1)

    try:
        run_alignment(args.ref, fastq, output_name, logger=logger)
        logger.info(f"Alignment completed successfully for sample: {args.sample_name}")
    except subprocess.CalledProcessError as e:
        # Log more details from the failed process
        logger.error(f"Alignment step failed with exit code {e.returncode}")
        if e.stdout:
             logger.error(f"Failed command stdout:\n{e.stdout}")
        if e.stderr:
             logger.error(f"Failed command stderr:\n{e.stderr}")
        sys.exit(1)
    except Exception as e:
        logger.exception(f"An unexpected error occurred during alignment:") 
        sys.exit(1)


if __name__ == "__main__":
    main()
