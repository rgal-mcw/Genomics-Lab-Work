# ./pipeline/pipes/spectre.py
# Description: Runs Mosdepth & Spectre CNV Caller. This is the
#              python adaptation of what was `ONT_spectre.sh`.
#              updated for new CW maple. 
# Author: Ryan Gallagher (Broeckel Lab, 2025)
# Date: 2025-04-28

import os
import subprocess
import argparse
import logging
import sys
import shutil

# -- Setup logging function (consistent with other scripts) ---
def setup_logging(log_file_path, log_level=logging.INFO):
    """Configures logging for the script."""
    # Use basicConfig to configure the root logger
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path, mode='a'), # Append mode
            logging.StreamHandler(sys.stdout) # Log to console
        ]
    )
    # Return a logger specific to this module
    return logging.getLogger(__name__)

# --- Function to check if a command exists ---
def command_exists(cmd):
    """Check if a command is available in the system PATH."""
    return shutil.which(cmd) is not None

# -- Main Function ---
def main():
    # --- Argument Parser Setup ---
    parser = argparse.ArgumentParser(description='Run Mosdepth and Spectre CNV Caller.')
    parser.add_argument('--input-dir', required=True,
                        help='Full path to the sample directory (e.g., /path/to/flowcell/sample_name/)')
    parser.add_argument('--ref-genome', required=True,
                        help='Path to the reference genome FASTA file.')
    parser.add_argument('--log-base-dir', required=True,
                        help='Base directory where the \'logs\' subdirectory should be placed (e.g., the main sample output dir)')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads to use for mosdepth.')

    args = parser.parse_args()

    # --- Setup Logging ---
    # Ensure log_base_dir exists (it should, if pipeline.py runs first)
    os.makedirs(args.log_base_dir, exist_ok=True)
    log_dir = os.path.join(args.log_base_dir, 'logs')
    os.makedirs(log_dir, exist_ok=True) # Ensure logs subdir exists
    log_file = os.path.join(log_dir, 'spectre.log')
    logger = setup_logging(log_file, log_level=logging.INFO)

    logger.info("--- Starting ONT Spectre Pipeline Script ---")
    logger.info(f"Input Directory: {args.input_dir}")
    logger.info(f"Reference Genome: {args.ref_genome}")
    logger.info(f"Log Base Directory: {args.log_base_dir}")
    logger.info(f"Threads for Mosdepth: {args.threads}")

    # --- Basic Input Validation ---
    if not os.path.isdir(args.input_dir):
        logger.error(f"Input directory does not exist: {args.input_dir}")
        sys.exit(1)

    # --- Get Sample Name ---
    sample_name = os.path.basename(os.path.normpath(args.input_dir)) # Handles trailing slashes
    logger.info(f"Derived Sample Name: {sample_name}")

    # --- Define Paths ---
    spectre_dir = os.path.join(args.input_dir, 'spectre')
    sorted_bam = os.path.join(args.input_dir, f"{sample_name}.sorted.bam")
    # Assuming PMDV ran first and created snp_indel directory
    snp_indel_dir = os.path.join(args.input_dir, 'snp_indel')
    # Standard PMDV output name (adjust if different)
    snp_vcf = os.path.join(snp_indel_dir, f"{sample_name}_PMDV_wgs.phased.vcf.gz")
    mosdepth_prefix = os.path.join(spectre_dir, sample_name)
    mosdepth_regions_bed = f"{mosdepth_prefix}.regions.bed.gz"

    # --- Create Spectre Directory ---
    try:
        os.makedirs(spectre_dir, exist_ok=True)
        logger.info(f"Ensured spectre output directory exists: {spectre_dir}")
    except OSError as e:
        logger.error(f"Failed to create spectre directory '{spectre_dir}': {e}")
        sys.exit(1)

    # --- Check for Required Commands (and implicitly, environment) ---
    required_commands = ["mosdepth", "spectre"]
    commands_ok = True
    for cmd in required_commands:
        if not command_exists(cmd):
            logger.error(f"Command '{cmd}' not found in PATH.")
            commands_ok = False
    if not commands_ok:
        logger.error("One or more required commands are missing.")
        logger.error("Please ensure the correct conda environment (e.g., 'spectre') is activated before running this script.")
        sys.exit(1)
    else:
        logger.info("Required commands (mosdepth, spectre) found in PATH.")


    # --- Check for Sorted BAM ---
    if not os.path.exists(sorted_bam):
         logger.error(f"Required sorted BAM file not found: {sorted_bam}")
         sys.exit(1)
    logger.info(f"Found input BAM: {sorted_bam}")

    # --- Run Mosdepth if needed ---
    if os.path.exists(mosdepth_regions_bed):
        logger.info(f"Mosdepth output already exists: {mosdepth_regions_bed}. Skipping mosdepth run.")
    else:
        logger.info("Mosdepth output not found. Running mosdepth...")
        mosdepth_command = [
            "mosdepth",
            "-t", str(args.threads),
            "-x",                          # Skip per-base depth calculation
            "-b", "1000",                  # BED file bin size
            "-Q", "20",                    # Mapping quality threshold
            mosdepth_prefix,               # Output prefix
            sorted_bam                     # Input BAM file
        ]
        logger.info(f"Executing Mosdepth command: {' '.join(mosdepth_command)}")
        try:
            result_mosdepth = subprocess.run(mosdepth_command, check=True, capture_output=True, text=True)
            if result_mosdepth.stdout:
                logger.info(f"Mosdepth stdout:\n{result_mosdepth.stdout.strip()}")
            if result_mosdepth.stderr:
                logger.warning(f"Mosdepth stderr:\n{result_mosdepth.stderr.strip()}") # Log stderr as warning
            logger.info("Mosdepth finished successfully.")
        except subprocess.CalledProcessError as e:
            logger.error(f"Mosdepth command failed with exit code {e.returncode}")
            if e.stdout: logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
            if e.stderr: logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
            sys.exit(e.returncode)
        except FileNotFoundError:
            logger.error("Error: 'mosdepth' command not found during execution attempt.")
            sys.exit(1)
        except Exception as e:
            logger.exception("An unexpected error occurred while running mosdepth:")
            sys.exit(1)

    # --- Check for SNV VCF from PMDV ---
    if not os.path.exists(snp_vcf):
        logger.error(f"Required SNV VCF file not found: {snp_vcf}")
        logger.error("Please ensure the PMDV step (Step 3) completed successfully.")
        sys.exit(1)
    logger.info(f"Found SNV VCF: {snp_vcf}")


    # --- Run Spectre CNVCaller ---
    logger.info("Running Spectre CNVCaller...")
    spectre_command = [
        "spectre", "CNVCaller",
        "--bin-size", "1000",
        "--coverage", spectre_dir,          # Directory containing mosdepth output
        "--sample-id", sample_name,
        "--output-dir", spectre_dir,        # Output directory for spectre results
        "--reference", args.ref_genome,
        "--snv", snp_vcf                   # Input SNV VCF path
        # Add any other spectre parameters if needed
    ]
    logger.info(f"Executing Spectre command: {' '.join(spectre_command)}")
    try:
        result_spectre = subprocess.run(spectre_command, check=True, capture_output=True, text=True)
        if result_spectre.stdout:
            logger.info(f"Spectre stdout:\n{result_spectre.stdout.strip()}")
        if result_spectre.stderr:
            logger.warning(f"Spectre stderr:\n{result_spectre.stderr.strip()}") # Log stderr as warning
        logger.info("Spectre CNVCaller finished successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Spectre command failed with exit code {e.returncode}")
        if e.stdout: logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
        if e.stderr: logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
        sys.exit(e.returncode)
    except FileNotFoundError:
        logger.error("Error: 'spectre' command not found during execution attempt.")
        sys.exit(1)
    except Exception as e:
        logger.exception("An unexpected error occurred while running spectre:")
        sys.exit(1)

    logger.info("--- ONT Spectre Pipeline Script Finished ---")

if __name__ == "__main__":
    main()
