# ./pipeline/pipes/snpeff.py 
# Description: This script runs SnpEff and SnpSift for annotation.
#              Replaces snpeff.sh for integration into pipeline.py
# Author: Ryan Gallagher (Broeckel Lab, 2025)
# Date: 2025-04-28

import os
import subprocess
import argparse
import logging
import sys
import shutil
import gzip

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
    """Check if a command (like java) is available in the system PATH."""
    return shutil.which(cmd) is not None

# --- Function to run shell commands ---
def run_command(command_list, logger, cwd=None):
    """Runs a shell command using subprocess and logs output."""
    command_str = ' '.join(command_list)
    logger.info(f"Executing command: {command_str}")
    try:
        result = subprocess.run(command_list, check=True, capture_output=True, text=True, cwd=cwd)
        if result.stdout:
            logger.info(f"Command stdout:\n{result.stdout.strip()}")
        if result.stderr:
            # Log stderr as warning even on success
            logger.warning(f"Command stderr:\n{result.stderr.strip()}")
        logger.info(f"Command finished successfully: {command_str}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}: {command_str}")
        if e.stdout:
            logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
        if e.stderr:
            logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
        return False
    except FileNotFoundError:
        logger.error(f"Error: Command not found: {command_list[0]}")
        return False
    except Exception as e:
        logger.exception(f"An unexpected error occurred while running command: {command_str}")
        return False


# -- Main Function -- #
def main():
    # --- Argument Parser Setup ---
    parser = argparse.ArgumentParser(description='Run SnpEff and SnpSift for SNP annotation.')
    parser.add_argument('--snp-indel-dir', required=True, help='Path to the directory containing the input VCF (e.g., .../snp_indel)')
    parser.add_argument('--vcf-gz', required=True, help='Filename of the input gzipped VCF (e.g., sample_PMDV_wgs.phased.vcf.gz)')
    parser.add_argument('--vcf-unzipped', required=True, help='Filename of the unzipped input VCF (e.g., sample_PMDV_wgs.phased.vcf)')
    parser.add_argument('--out-prefix', required=True, help='Output filename prefix (e.g., sample_PMDV_wgs.phased.snpeff)')
    parser.add_argument('--log-base-dir', required=True, help='Base directory where the \'logs\' subdirectory is located')
    parser.add_argument('--snpeff-base-dir', default=os.path.expanduser('/data/ref/snpEff'), help='Base directory containing snpEff.jar and SnpSift.jar')
    parser.add_argument('--ref-genome-version', default='GRCh38.105', help='SnpEff reference genome version (e.g., GRCh38.105)')
    parser.add_argument('--clinvar-vcf', default='/data/ref/ClinVar/clinvar.vcf.gz', help='Path to the ClinVar VCF file')
    # Removed DBSNP annotation step as it was commented out/unclear in the original .sh reference provided previously, but can be added back if needed.
    # parser.add_argument('--dbsnp-vcf', default='/data/ref/dbSNP/00-All.vcf.gz', help='Path to the dbSNP VCF file')


    args = parser.parse_args()

    # --- Setup Logging ---
    log_dir = os.path.join(args.log_base_dir, 'logs')
    # Ensure logs directory exists (pipeline.py should create it, but check again)
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, 'snpeff.log')
    logger = setup_logging(log_file, log_level=logging.INFO)

    logger.info("--- Starting SnpEff Annotation Script ---")
    logger.info(f"Input VCF Dir: {args.snp_indel_dir}")
    logger.info(f"Input VCF (gz): {args.vcf_gz}")
    logger.info(f"Input VCF (unzipped): {args.vcf_unzipped}")
    logger.info(f"Output Prefix: {args.out_prefix}")
    logger.info(f"Log Base Directory: {args.log_base_dir}")
    logger.info(f"SnpEff Base Directory: {args.snpeff_base_dir}")
    logger.info(f"Reference Genome Version: {args.ref_genome_version}")
    logger.info(f"ClinVar VCF: {args.clinvar_vcf}")
    # logger.info(f"dbSNP VCF: {args.dbsnp_vcf}") # If using dbSNP

    # --- Define Paths ---
    snpeff_jar = os.path.join(args.snpeff_base_dir, 'snpEff.jar')
    snpsift_jar = os.path.join(args.snpeff_base_dir, 'SnpSift.jar')
    input_vcf_gz_path = os.path.join(args.snp_indel_dir, args.vcf_gz)
    input_vcf_path = os.path.join(args.snp_indel_dir, args.vcf_unzipped)
    snpeff_output_vcf = os.path.join(args.snp_indel_dir, f"{args.out_prefix}.vcf")
    # Define final output name based on the last step in the original script
    final_output_vcf = os.path.join(args.snp_indel_dir, f"{args.out_prefix}.clinvar.vcf") # Changed from rsID.clinvar

    # --- Input Checks ---
    if not command_exists('java'):
        logger.error("Java command not found. Please ensure Java is installed and in the PATH.")
        sys.exit(1)
    logger.info("Java command found.")

    if not os.path.isfile(snpeff_jar):
        logger.error(f"SnpEff jar not found at: {snpeff_jar}")
        sys.exit(1)
    logger.info(f"Found SnpEff jar: {snpeff_jar}")

    if not os.path.isfile(snpsift_jar):
        logger.error(f"SnpSift jar not found at: {snpsift_jar}")
        sys.exit(1)
    logger.info(f"Found SnpSift jar: {snpsift_jar}")

    if not os.path.isfile(input_vcf_gz_path):
        logger.error(f"Input gzipped VCF not found: {input_vcf_gz_path}")
        sys.exit(1)
    logger.info(f"Found input gzipped VCF: {input_vcf_gz_path}")

    if not os.path.isfile(args.clinvar_vcf):
        logger.warning(f"ClinVar VCF not found: {args.clinvar_vcf}. Skipping ClinVar annotation.")
        # Set ClinVar VCF to None to skip the step later
        args.clinvar_vcf = None
    else:
         logger.info(f"Found ClinVar VCF: {args.clinvar_vcf}")


    # --- Unzip VCF if necessary ---
    if not os.path.isfile(input_vcf_path):
        logger.info(f"Unzipped VCF not found ({input_vcf_path}). Unzipping {input_vcf_gz_path}...")
        try:
            with gzip.open(input_vcf_gz_path, 'rb') as f_in:
                with open(input_vcf_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logger.info(f"Successfully unzipped VCF to: {input_vcf_path}")
        except Exception as e:
            logger.error(f"Failed to unzip {input_vcf_gz_path}: {e}")
            sys.exit(1)
    else:
        logger.info(f"Unzipped VCF already exists: {input_vcf_path}")

    # --- Run SnpEff ---
    logger.info(f"Running SnpEff on {input_vcf_path}...")
    snpeff_command = [
        "java", "-Xmx4g", "-jar", snpeff_jar, # Added memory flag, adjust if needed
        args.ref_genome_version,
        input_vcf_path
    ]
    # Run SnpEff directing stdout to the output file
    logger.info(f"Executing command (stdout to file): {' '.join(snpeff_command)} > {snpeff_output_vcf}")
    try:
        with open(snpeff_output_vcf, 'w') as outfile:
            result = subprocess.run(snpeff_command, check=True, stdout=outfile, stderr=subprocess.PIPE, text=True, cwd=args.snp_indel_dir) # Run in snp_indel_dir
        if result.stderr:
             logger.warning(f"SnpEff stderr:\n{result.stderr.strip()}")
        logger.info(f"SnpEff finished successfully. Output: {snpeff_output_vcf}")

        # Move SnpEff summary files (original script moved them from pwd)
        # Assuming SnpEff creates these in the CWD where java was called (snp_indel_dir)
        summary_html = os.path.join(args.snp_indel_dir, "snpEff_summary.html")
        genes_txt = os.path.join(args.snp_indel_dir, "snpEff_genes.txt")

        if os.path.exists(summary_html):
            logger.info(f"Found {summary_html}")
        else:
             logger.warning(f"snpEff_summary.html not found in {args.snp_indel_dir}.")

        if os.path.exists(genes_txt):
             logger.info(f"Found {genes_txt}")
        else:
             logger.warning(f"snpEff_genes.txt not found in {args.snp_indel_dir}.")

    except subprocess.CalledProcessError as e:
        logger.error(f"SnpEff command failed with exit code {e.returncode}")
        if e.stderr:
            logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
        # Attempt to read stdout from file if it exists
        if os.path.exists(snpeff_output_vcf):
             try:
                 with open(snpeff_output_vcf, 'r') as f:
                      stdout_content = f.read()
                      if stdout_content:
                           logger.error(f"Failed command partial stdout (from file):\n{stdout_content.strip()}")
             except Exception as read_err:
                  logger.error(f"Could not read partial output file {snpeff_output_vcf}: {read_err}")
        sys.exit(e.returncode)
    except Exception as e:
        logger.exception("An unexpected error occurred while running SnpEff:")
        sys.exit(1)

    # --- Run SnpSift Annotate (ClinVar) ---
    if args.clinvar_vcf and os.path.exists(snpeff_output_vcf):
        logger.info(f"Running SnpSift annotate with ClinVar: {args.clinvar_vcf}...")
        snpsift_command = [
            "java", "-Xmx4g", "-jar", snpsift_jar, "annotate", # Added memory flag
            args.clinvar_vcf,
            snpeff_output_vcf # Input is the output from SnpEff
        ]
        # Run SnpSift directing stdout to the final output file
        logger.info(f"Executing command (stdout to file): {' '.join(snpsift_command)} > {final_output_vcf}")
        try:
            with open(final_output_vcf, 'w') as outfile:
                 result = subprocess.run(snpsift_command, check=True, stdout=outfile, stderr=subprocess.PIPE, text=True, cwd=args.snp_indel_dir)
            if result.stderr:
                 logger.warning(f"SnpSift stderr:\n{result.stderr.strip()}")
            logger.info(f"SnpSift annotate finished successfully. Output: {final_output_vcf}")
            try:
                os.remove(snpeff_output_vcf)
                logger.info(f"Removed intermediate file: {snpeff_output_vcf}")
            except OSError as e:
                logger.warning(f"Could not remove intermediate file {snpeff_output_vcf}: {e}")

        except subprocess.CalledProcessError as e:
            logger.error(f"SnpSift annotate command failed with exit code {e.returncode}")
            if e.stderr:
                 logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
            # Attempt to read stdout from file
            if os.path.exists(final_output_vcf):
                 try:
                      with open(final_output_vcf, 'r') as f:
                           stdout_content = f.read()
                           if stdout_content:
                                logger.error(f"Failed command partial stdout (from file):\n{stdout_content.strip()}")
                 except Exception as read_err:
                      logger.error(f"Could not read partial output file {final_output_vcf}: {read_err}")

            sys.exit(e.returncode)
        except Exception as e:
            logger.exception("An unexpected error occurred while running SnpSift annotate:")
            sys.exit(1)
    elif not args.clinvar_vcf:
        logger.info("Skipping ClinVar annotation as VCF path was not found or provided.")
        # If ClinVar is skipped, rename the SnpEff output to be the "final" output for consistency downstream
        try:
             if os.path.exists(snpeff_output_vcf):
                  shutil.move(snpeff_output_vcf, final_output_vcf)
                  logger.info(f"Renamed {snpeff_output_vcf} to {final_output_vcf} as ClinVar step was skipped.")
             else:
                  logger.error("SnpEff output file does not exist, cannot rename.")
                  sys.exit(1) # Exit if the expected SnpEff output isn't there
        except Exception as e:
             logger.error(f"Failed to rename SnpEff output file: {e}")
             sys.exit(1)

    elif not os.path.exists(snpeff_output_vcf):
         logger.error(f"Input VCF for SnpSift annotate not found: {snpeff_output_vcf}. Cannot proceed.")
         sys.exit(1)


    logger.info("--- SnpEff Annotation Script Finished ---")

if __name__ == "__main__":
    main()
