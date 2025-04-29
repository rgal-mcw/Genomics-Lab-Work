# Title: Sniffles2 for pipeline.py 
# Author: Ryan Gallagher, Broeckel Lab 2025
# Date: 04.28.2025 (Updated for CW)
#
# Description: This script runs Sniffles2 on the SV-Calling step of the pipeline.
#              This script will run the command line tool Sniffles given the parsed
#              arguments from pipeline.py..
#
#


import os
import subprocess
import argparse
import logging
import sys
import shutil

# -- Setup logging function (consistent with other scripts) ---
def setup_logging(log_file_path, log_level=logging.INFO):
    """Configures logging for the script."""
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path, mode='a'), # Append mode
            logging.StreamHandler(sys.stdout) # Log to console
        ]
    )
    return logging.getLogger(__name__)

# --- Function to check if a command exists ---
def command_exists(cmd):
    """Check if a command (like sniffles) is available in the system PATH."""
    return shutil.which(cmd) is not None

# -- Main Function ---
def main():
    # --- Argument Parser Setup ---
    parser = argparse.ArgumentParser(description='Run Sniffles2 SV Caller.')
    parser.add_argument('--input-bam', required=True, help='Path to the input BAM file (coordinate sorted)')
    parser.add_argument('--output-vcf', required=True, help='Path for the output VCF file')
    # Consider adding --output-snf for the SNF file if needed downstream
    # parser.add_argument('--output-snf', required=True, help='Path for the output SNF file')
    parser.add_argument('--ref', required=True, help='Path to the reference genome FASTA file')
    parser.add_argument('--sample-id', required=True, help='Sample ID to embed in the VCF')
    parser.add_argument('--log-base-dir', required=True, help='Base directory where the \'logs\' subdirectory is located')
    parser.add_argument('--threads', type=int, default=16, help='Number of threads for Sniffles2')
    parser.add_argument('--timeout', type=int, help='Optional timeout in seconds for the Sniffles2 process')


    args = parser.parse_args()

    # --- Setup Logging ---
    log_dir = os.path.join(args.log_base_dir, 'logs')
    # Ensure logs directory exists
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, 'sniffles2.log')
    logger = setup_logging(log_file, log_level=logging.INFO)

    logger.info("--- Starting Sniffles2 SV Calling Script ---")
    logger.info(f"Input BAM: {args.input_bam}")
    logger.info(f"Output VCF: {args.output_vcf}")
    logger.info(f"Reference: {args.ref}")
    logger.info(f"Sample ID: {args.sample_id}")
    logger.info(f"Log Base Directory: {args.log_base_dir}")
    logger.info(f"Threads: {args.threads}")
    if args.timeout:
        logger.info(f"Timeout set: {args.timeout} seconds")


    # --- Prerequisite Checks ---
    if not command_exists('sniffles'):
        logger.error("Sniffles command not found. Please ensure Sniffles2 is installed and in the PATH.")
        sys.exit(1)
    logger.info("Sniffles command found.")

    if not os.path.isfile(args.input_bam):
        logger.error(f"Input BAM file not found: {args.input_bam}")
        sys.exit(1)
    logger.info(f"Found input BAM: {args.input_bam}")

    if not os.path.isfile(args.ref):
        logger.error(f"Reference FASTA file not found: {args.ref}")
        sys.exit(1)
    logger.info(f"Found reference FASTA: {args.ref}")

    # Ensure output directory exists
    output_dir = os.path.dirname(args.output_vcf)
    if output_dir: # Handle cases where output is in current dir
        os.makedirs(output_dir, exist_ok=True)
        logger.info(f"Ensured output directory exists: {output_dir}")


    # --- Construct Sniffles2 Command ---
    # Using parameters mostly from sniffles2.sh
    # and the original sniffles2.py
    sniffles_command = [
        'sniffles',
        '--input', args.input_bam,
        '--reference', args.ref,
        '--vcf', args.output_vcf,
        # '--snf', args.output_snf, # Add if using --output-snf arg
        '--threads', str(args.threads),
        '--sample-id', args.sample_id,
        # --- Adding other parameters from sniffles2.sh/py ---
        '--minsupport', 'auto',
        '--minsvlen', '35',
        '--minsvlen-screen-ratio', '0.9',
        '--mapq', '25',
        '--qc-stdev', 'True', # Note: These boolean flags might just need presence, not 'True'
        '--qc-stdev-abs-max', '500',
        '--qc-coverage', '1',
        '--long-ins-length', '2500',
        '--long-del-length', '50000',
        '--long-del-coverage', '0.66',
        '--long-dup-length', '50000',
        '--long-dup-coverage', '1.33',
        '--max-splits-kb', '0.1',
        '--max-splits-base', '3',
        '--min-alignment-length', '1000',
        '--phase-conflict-threshold', '0.1',
        '--detect-large-ins', 'True', # Boolean flag
        '--cluster-binsize', '100',
        '--cluster-r', '2.5',
        '--cluster-repeat-h', '1.5',
        '--cluster-repeat-h-max', '1000',
        '--cluster-merge-pos', '150',
        '--cluster-merge-len', '0.33',
        '--cluster-merge-bnd', '1500',
        '--genotype-ploidy', '2',
        '--genotype-error', '0.05',
        '--allow-overwrite', # Flag
        '--symbolic', # Flag
        '--max-del-seq-len', '50000'
    ]
    # Remove boolean flags that don't take 'True' if Sniffles expects just the flag
    # e.g., replace '--qc-stdev', 'True' with just '--qc-stdev' if needed. Check Sniffles2 docs.

    # --- Execute Sniffles2 Command ---
    logger.info("Executing Sniffles2 command:")
    logger.info(' '.join(sniffles_command))

    try:
        # Execute the command, capture output, check for errors, add timeout
        result = subprocess.run(
            sniffles_command,
            check=True,
            capture_output=True,
            text=True,
            timeout=args.timeout # Pass timeout value if provided
        )

        # Log stdout and stderr from the sniffles command
        if result.stdout:
            logger.info(f"Sniffles2 stdout:\n{result.stdout.strip()}")
        if result.stderr:
            # Log stderr as warning even on success
            logger.warning(f"Sniffles2 stderr:\n{result.stderr.strip()}")

        logger.info(f"Sniffles2 command completed successfully. Output VCF: {args.output_vcf}")

    except subprocess.TimeoutExpired as e:
        # Handle the timeout specifically
        logger.warning(f"Sniffles2 command timed out after {args.timeout} seconds.")
        logger.warning("Sniffles2 process might not have terminated cleanly.")
        if e.stdout: logger.warning(f"Timeout stdout:\n{e.stdout.strip()}")
        if e.stderr: logger.warning(f"Timeout stderr:\n{e.stderr.strip()}")
        # The pipeline.py script should still handle killing any lingering process
        # Exit with a specific code or allow pipeline.py to handle it? Let's exit non-zero.
        sys.exit(1) # Indicate timeout occurred

    except subprocess.CalledProcessError as e:
        # Log detailed error information if the command fails
        logger.error(f"Sniffles2 command failed with exit code {e.returncode}")
        if e.stdout: logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
        if e.stderr: logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
        sys.exit(e.returncode) # Exit with the same error code

    except Exception as e:
        # Catch any other unexpected errors
        logger.exception("An unexpected error occurred while running Sniffles2:")
        sys.exit(1)

    logger.info("--- Sniffles2 SV Calling Script Finished ---")

if __name__ == "__main__":
    main()
