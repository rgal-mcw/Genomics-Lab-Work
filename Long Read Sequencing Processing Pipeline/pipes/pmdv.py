# PMDV Runner
# Author: Ryan Gallagher (Broeckel Lab, 2025)
# Date: 2025-04-25
# Description: Python script to run the Pepper-Margin-DeepVariant
#              pipeline via its Docker container. Replace PMDV_pipeline.sh
#              for consistent logging and argument handling.


import os
import subprocess
import argparse
import logging
import sys

# -- Setup logging
def setup_logging(log_file_path, log_level=logging.INFO):
    """Configures logging for the script."""
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path, mode='a'), # Append mode
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


# -- Main Function
def main():
    # --- Argument Parser Setup ---
    parser = argparse.ArgumentParser(description='Run Pepper-Margin-DeepVariant via Docker.')
    parser.add_argument('--sample', required=True, help='Sample name (used for output prefix)')
    parser.add_argument('--input-dir', required=True, help='Base input directory containing the BAM file (mounted to Docker)')
    parser.add_argument('--output-dir', required=True, help='Output directory for PMDV results (e.g., .../snp_indel/)')
    parser.add_argument('--ref', required=True, help='Path to the reference FASTA file')
    parser.add_argument('--bam', required=True, help='Path to the input BAM file (relative to input-dir if mounted)')
    parser.add_argument('--chem', default='ont_r10_q20', help='ONT Chemistry/Model (e.g., ont_r10_q20, ont_r9_guppy5_sup)')
    parser.add_argument('--threads', type=int, default=32, help='Number of threads to use')
    # Argument to derive log directory location
    parser.add_argument('--log-base-dir', required=True, help='Base directory where the \'logs\' subdirectory should be placed (e.g., the main sample output dir)')

    args = parser.parse_args()

    # --- Setup Logging ---
    log_dir = os.path.join(args.log_base_dir, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, 'pmdv.log')
    logger = setup_logging(log_file, log_level=logging.INFO)
    

    # --- Communicate Parameters
    logger.info("--- Starting PMDV Pipeline Script ---")
    logger.info(f"Sample: {args.sample}")
    logger.info(f"Input Directory: {args.input_dir}")
    logger.info(f"Output Directory: {args.output_dir}")
    logger.info(f"Reference FASTA: {args.ref}")
    logger.info(f"Input BAM: {args.bam}")
    logger.info(f"Chemistry: {args.chem}")
    logger.info(f"Threads: {args.threads}")
    logger.info(f"Log File: {log_file}")

    # --- Prepare Paths and Parameters
    # Ensure Output dir exists
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        logger.info(f"Ensured output directory exists: {args.output_dir}")
    except OSError as e:
        logger.error(f"Failed to create output directory '{args.output_dir}': {e}")
        sys.exit(1)

    # Get the directory containing the reference for mounting
    ref_dir = os.path.dirname(args.ref)
    logger.info(f"Reference Directory (for mounting): {ref_dir}")

    # Define the output prefix for PMDV
    output_prefix = f"{args.sample}_PMDV_wgs"
    logger.info(f"Output Prefix: {output_prefix}")

    # Get current user and group ID for docker --user flag
    try:
        uid = os.getuid()
        gid = os.getgid()
        user_map = f"{uid}:{gid}"
        logger.info(f"Running Docker as user:group {user_map}")
    except AttributeError:
        user_map = None
        logger.warning("Could not determine UID/GID. Docker will run as default user.")

    

    # --- Construct Docker Command as a List
    docker_command = [
        "docker", "run", "--rm", # Run and remove container afterwards
        "-v", f"{args.input_dir}:{args.input_dir}",   # Mount input dir
        "-v", f"{args.output_dir}:{args.output_dir}", # Mount output dir
        "-v", f"{ref_dir}:{ref_dir}",                 # Mount reference dir
    ]

    if user_map:
        docker_command.extend(["-u", user_map])

    # Add Docker image and PMDV command/arguments
    docker_command.extend([
        "kishwars/pepper_deepvariant:r0.8",           # Docker image
        "run_pepper_margin_deepvariant", "call_variant", # PMDV command
        "-b", args.bam,                               # Input BAM path (inside container)
        "-f", args.ref,                                # Reference FASTA path (inside container)
        "-o", args.output_dir,                         # Output directory (inside container)
        "-p", output_prefix,                           # Output prefix
        "-t", str(args.threads),                       # Threads
        "--pepper_min_mapq", "15",
        "--pepper_min_snp_baseq", "15",
        "--pepper_min_indel_baseq", "15",
        "--dv_min_mapping_quality", "15",
        "--dv_min_base_quality", "15",
        "--phased_output",
        f"--{args.chem}"                               # Chemistry flag
    ])

    logger.info("Constructed docker command:")
    # Log command for debugging, quoting arguments with spaces if necessary
    logger.info(' '.join(subprocess.list2cmdline(docker_command)))

        # --- Execute Docker Command ---
    logger.info("Executing Pepper-Margin-DeepVariant Docker command...")
    try:
        # Execute the command, capture output, check for errors
        result = subprocess.run(docker_command, check=True, capture_output=True, text=True)

        # Log stdout and stderr from the docker command
        if result.stdout:
            logger.info(f"PMDV Docker stdout:\n{result.stdout.strip()}")
        if result.stderr:
            # Log stderr as warning even on success, as tools might print info here
            logger.warning(f"PMDV Docker stderr:\n{result.stderr.strip()}")

        logger.info("Pepper-Margin-DeepVariant Docker command completed successfully.")

    except subprocess.CalledProcessError as e:
        # Log detailed error information if the command fails
        logger.error(f"Pepper-Margin-DeepVariant Docker command failed with exit code {e.returncode}")
        if e.stdout:
            logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
        if e.stderr:
            logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
        sys.exit(e.returncode) # Exit with the same error code as the failed process

    except FileNotFoundError:
        logger.error("Error: 'docker' command not found. Is Docker installed and in the PATH?")
        sys.exit(1)

    except Exception as e:
        # Catch any other unexpected errors during subprocess execution
        logger.exception("An unexpected error occurred while running the Docker command:")
        sys.exit(1)

    logger.info("--- PMDV Pipeline Script Finished ---")

if __name__ == "__main__":
    main()


