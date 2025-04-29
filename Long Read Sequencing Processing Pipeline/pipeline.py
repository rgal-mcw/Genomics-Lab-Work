################################################################################
#   ____  _   _ _______   _____ _____ _____  ______ _      _____ _   _ ______  #
#  / __ \| \ | |__   __| |  __ \_   _|  __ \|  ____| |    |_   _| \ | |  ____| #
# | |  | |  \| |  | |    | |__) || | | |__) | |__  | |      | | |  \| | |__    #
# | |  | | . ` |  | |    |  ___/ | | |  ___/|  __| | |      | | | . ` |  __|   #
# | |__| | |\  |  | |    | |    _| |_| |    | |____| |____ _| |_| |\  | |____  #
#  \____/|_| \_|  |_|    |_|   |_____|_|    |______|______|_____|_| \_|______| #
#                                                                              #
################################################################################
#                                                                              #
# Broeckle Lab Bioinformatics, 2024                                            #
#                                                                              #
# This file executes the pipeline for the ONT PromethION sequencing data by    #
# executing tools at the prompt of the user. These scripts are stored in sub-  #
# directories of the main directory, and can be edited and updated at a whim.  #
#                                                                              #
# The current pipline is as follows:                                           #
#  1. Prepare data for transfer                                                # 
#    -> 1.2. Transfer data from prom to Maple                                  #
#    -> 1.3 Concatinate FASTQ files                                            #
#  2. Read Alignment (FASTQs w/ GRCh38)                                        #
#  3. SNP/INDEL via PMDV                                                       #    
#  4. Large SV via Spectre                                                     #     
#  5. SNP Annotation on .phased.vcf file                                       #
#       -> Filter snpEff VCF                                                   #
#  6. SV-Calling via SNIFFLES2                                                 #
#                                                                              #
################################################################################
################################################################################

# Notes for development:
#  - Should make a custom conda environment for this pipeline
#       -> This should check for that
#  - Sequencing summary file naming? Not *_combined.txt DONE

import os
import sys
import subprocess
import glob
import argparse
import logging

def end_in_slash(dir):
    if not dir.endswith('/'):
        dir += '/'
    return dir

def main():


        # --- Updated Input Section using argparse---
    parser = argparse.ArgumentParser(description='ONT PromethION Data Processing Pipeline.')
    parser.add_argument('sample_name', type=str, 
                        help='The name of the sample being processed. This should match a created directory.')
    parser.add_argument('data_dir', type=str, 
                        help='Full path to the directory containing the sample data on the PromethION server (e.g., /data/run_folder/sample_subfolder/).')
    parser.add_argument('to_dir', type=str, 
                        help='Full path to the MAPLE directory where output should be stored (e.g., /data2/flowcell_10.4.1/mcw_svi_.../).')
    parser.add_argument('--log-level', default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the logging level')
    
    args = parser.parse_args()

    # Assign parsed arguments to variables, keeping original variable names
    sample_name = args.sample_name
    data_dir = end_in_slash(args.data_dir) # Apply end_in_slash here
    to_dir = end_in_slash(args.to_dir)     # Apply end_in_slash here


        # --- Updating Logging --- #
    log_dir = os.path.join(to_dir, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    log_file_path = os.path.join(log_dir, 'pipeline.log')
    
    # Configure Logging
    log_level = getattr(logging, args.log_level.upper(), logging.INFO) # Get level from args
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s', # Add timestamp and level
        handlers=[
            logging.FileHandler(log_file_path, mode='w'), # Overwrite log file each run
            logging.StreamHandler(sys.stdout) # Log to console
        ]
    )
    logger = logging.getLogger()
    
    class StreamToLogger:
        """Fake file-like stream object that redirects writes to a logger instance."""
        def __init__(self, logger_instance, log_level=logging.ERROR):
            self.logger = logger_instance
            self.log_level = log_level
            self.linebuf = ''

        def write(self, buf):
            for line in buf.rstrip().splitlines():
                self.logger.log(self.log_level, line.rstrip())

        def flush(self):
            pass # Required for file-like object interface

    sys.stderr = StreamToLogger(logger, logging.ERROR)

    # Now, we have to do any print() statements to logger.info() statements
    # --- End Logging Setup ---

    ref = '/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta'
    logger.info(f'Reference Genome set to:\n{ref}')

    if not os.path.exists(to_dir):
        logger.error("The MAPLE directory does not exist. Exiting script.")
        sys.exit(1)
    logger.info('The MAPLE directory exists.')
    log_dir = os.path.join(to_dir, 'logs') # Make sure log_dir is defined


    #----------------------
    # Step 1: Prepare data for Transfer and send to Maple
    logger.info(f"Running Step 1: File Transfer for {sample_name} using file_transfer.py")
    transfer_command = ["python", "./pipes/file_transfer.py", sample_name, data_dir, to_dir]
    logger.info("Executing Command:\n" + ' '.join(transfer_command))
    try:
        # Run the script, capture output, check return code
        result_transfer = subprocess.run(transfer_command, check=True, capture_output=True, text=True)
        # Log stdout/stderr from the script for pipeline context (optional, as script logs itself)
        if result_transfer.stdout:
            logger.info(f"file_transfer.py stdout:\n{result_transfer.stdout.strip()}")
        if result_transfer.stderr:
            logger.warning(f"file_transfer.py stderr:\n{result_transfer.stderr.strip()}") # Log stderr as warning
        logger.info("File transfer script executed successfully.")

    except subprocess.CalledProcessError as e:
        logger.error(f"Execution of file_transfer.py failed with exit code {e.returncode}")
        if e.stdout: logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
        if e.stderr: logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
        sys.exit(f"File transfer step (Step 1) failed with exit code {e.returncode}") # Exit pipeline
    except FileNotFoundError:
         logger.error("Error: 'python' or './pipes/file_transfer.py' not found.")
         sys.exit("Failed to execute file transfer script.")
    except Exception as e:
         logger.exception("An unexpected error occurred during the file transfer step:")
         sys.exit("Unexpected error in file transfer step (Step 1).")

    logger.info("%%%%%%%%%%%%%\nFINISHED STEP 1: DATA TRANSFER \n%%%%%%%%%%%%%")
    #----------------------

    # Step 2: Read Alignment (FASTQs w/ GRCh38)
    # Note: alignment.py checks for the fastq file internally
    logger.info(f"Running Step 2: Read Alignment for {sample_name} using alignment.py")
    alignment_command = ["python", "./pipes/alignment.py", sample_name, to_dir, ref]
    logger.info("Executing Command:\n" + ' '.join(alignment_command))
    try:
        # Run the script, capture output, check return code
        result_alignment = subprocess.run(alignment_command, check=True, capture_output=True, text=True)
        # Log stdout/stderr from the script for pipeline context (optional)
        if result_alignment.stdout:
            logger.info(f"alignment.py stdout:\n{result_alignment.stdout.strip()}")
        if result_alignment.stderr:
            logger.warning(f"alignment.py stderr:\n{result_alignment.stderr.strip()}")
        logger.info("Alignment script executed successfully.")

    except subprocess.CalledProcessError as e:
        logger.error(f"Execution of alignment.py failed with exit code {e.returncode}")
        if e.stdout: logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
        if e.stderr: logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
        sys.exit(f"Alignment step (Step 2) failed with exit code {e.returncode}") # Exit pipeline
    except FileNotFoundError:
         logger.error("Error: 'python' or './pipes/alignment.py' not found.")
         sys.exit("Failed to execute alignment script.")
    except Exception as e:
         logger.exception("An unexpected error occurred during the alignment step:")
         sys.exit("Unexpected error in alignment step (Step 2).")

    logger.info("%%%%%%%%%%%%%\n FINISHED STEP 2: ALIGNMENT \n%%%%%%%%%%%%%")
    logger.info('Next Step: SNP/INDEL via PMDV') # Message indicating next step remains

    #----------------------
    # Step 3: SNP/INDEL via PMDV

    if os.path.exists(to_dir + f"{sample_name}.sorted.bam"):
        logger.info(f"Running PMDV on {sample_name} using pmdv.py")
        snp_indel_dir = os.path.join(to_dir, "snp_indel") # Use os.path.join for consistency
        bam_file = os.path.join(to_dir, f"{sample_name}.sorted.bam") # Use os.path.join

        # Ensure the output directory exists (pmdv.py also does this, but good practice)
        os.makedirs(snp_indel_dir, exist_ok=True)

        # Construct the command to run pmdv.py
        pmdv_command = [
            "python", "./pipes/pmdv.py",
            "--sample", sample_name,
            "--input-dir", to_dir,          # Base directory containing the BAM
            "--output-dir", snp_indel_dir,  # Specific output directory for PMDV
            "--ref", ref,
            "--bam", bam_file,
            "--log-base-dir", to_dir       # Directory where 'logs' subdir resides
            # Optional: Add --chem or --threads if you need non-default values
            # "--chem", "your_chemistry",
            # "--threads", "your_thread_count"
        ]

        # Log the command being run
        logger.info("Running Command:\n" + ' '.join(pmdv_command))

        # Execute the pmdv.py script
        # Note: pmdv.py handles its own logging to pmdv.log within the log_dir
        try:
            result_pmdv = subprocess.run(pmdv_command, check=True, capture_output=True, text=True)
            # Log stdout/stderr from pmdv.py if needed (it already logs to file)
            if result_pmdv.stdout:
                logger.info(f"pmdv.py stdout:\n{result_pmdv.stdout.strip()}")
            if result_pmdv.stderr:
                logger.warning(f"pmdv.py stderr:\n{result_pmdv.stderr.strip()}") # Log stderr as warning
            logger.info("pmdv.py script executed successfully.")

            # Index the haplotagged BAM produced by pmdv.py
            haplotagged_bam = os.path.join(snp_indel_dir, f"{sample_name}_PMDV_wgs.haplotagged.bam")
            if os.path.exists(haplotagged_bam):
                 index_command = f"samtools index -@ 8 '{haplotagged_bam}'"
                 logger.info(f"Indexing haplotagged BAM: {index_command}")
                 result_samtools = subprocess.run(index_command, shell=True, check=True, capture_output=True, text=True)
                 logger.info("Haplotagged BAM indexing completed.")
                 if result_samtools.stderr:
                     logger.warning(f"samtools index stderr:\n{result_samtools.stderr.strip()}")
            else:
                logger.warning(f"Haplotagged BAM not found for indexing: {haplotagged_bam}")

        except subprocess.CalledProcessError as e:
            logger.error(f"Execution of {' '.join(e.cmd)} failed with exit code {e.returncode}")
            if e.stdout:
                logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
            if e.stderr:
                logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
            # Decide if you want to exit the pipeline here
            # sys.exit(f"PMDV step failed with exit code {e.returncode}")
        except Exception as e:
             logger.exception("An unexpected error occurred during the PMDV step:")
             # Decide if you want to exit
             # sys.exit("Unexpected error in PMDV step.")


    else:
        # This condition remains from the original code
        logger.error(f"Input BAM file not found: {os.path.join(to_dir, f'{sample_name}.sorted.bam')}. Skipping PMDV.")
        # Decide if you want to exit here
        # sys.exit("Input BAM for PMDV not found.")


    logger.info("%%%%%%%%%%%%%\n FINISHED PMDV SNP/INDEL CALLING \n%%%%%%%%%%%%%")
    logger.info('Next Step: CNV via Spectre')
    


    #----------------------
    # Step 4: Large SV via Spectre 

    logger.info("NOTE: A dedicated spectre conda environment needs to be created for this part of the pipeline. If this environment is not created yet, do so via instructions from the ont-spectre github documentation.")

    logger.info(f"Running Spectre step for {sample_name} using ont_spectre.py")
    spectre_py_command = [
        "python", "./pipes/ont_spectre.py",
        "--input-dir", to_dir,
        "--ref-genome", ref,
        "--log-base-dir", to_dir
        # Optional: Add --threads if you want to control mosdepth threads from pipeline.py
        # "--threads", "16"
]

    logger.info("Executing Command:\n" + ' '.join(spectre_py_command))
    try:
        # Assuming the 'spectre' conda env is active before pipeline.py is run
        spectre_result = subprocess.run(spectre_py_command, check=True, capture_output=True, text=True)
        # Log output from ont_spectre.py within pipeline.log for context
        if spectre_result.stdout:
            logger.info(f"ont_spectre.py stdout:\n{spectre_result.stdout.strip()}")
        if spectre_result.stderr:
            logger.warning(f"ont_spectre.py stderr:\n{spectre_result.stderr.strip()}")
        logger.info("ont_spectre.py executed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Execution of ont_spectre.py failed with exit code {e.returncode}")
        if e.stdout: logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
        if e.stderr: logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
        # Decide if you want to exit pipeline.py here
        # sys.exit(f"Spectre step failed.")
    except Exception as e:
        logger.exception("An unexpected error occurred during the Spectre step:")
        # Decide if you want to exit
        # sys.exit("Unexpected error in Spectre step.")


    logger.info("%%%%%%%%%%%%%\n FINISHED SPECTRE \n%%%%%%%%%%%%%")
    logger.info('Next Step: SNP Annotation')
    #----------------------

    # Step 5: SNP Annotation on .phased.vcf file

        # Define paths and prefixes
    snp_indel_dir = os.path.join(to_dir, "snp_indel") # Use os.path.join for consistency
    gz = f"{sample_name}_PMDV_wgs.phased.vcf.gz"
    vcf_unzipped_relative = f"{sample_name}_PMDV_wgs.phased.vcf" # Relative name
    out_snpeff = f"{sample_name}_PMDV_wgs.phased.snpeff" # Prefix for output

    vcf_unzipped_abs_path = os.path.abspath(os.path.join(snp_indel_dir, vcf_unzipped_relative))

    # Define expected location of snpEff (adjust if different)
    # Using the path revealed in your previous error log:
    snpeff_base_dir = "/data/ref/snpEff"
    snpeff_jar_path = os.path.join(snpeff_base_dir, 'snpEff.jar')

    # Check if snpEff.jar exists before attempting to run
    if os.path.exists(snpeff_jar_path):
        logger.info(f"Running snpEff annotation for {sample_name} using snpeff.py")
        logger.info(f"Using absolute path for VCF: {vcf_unzipped_abs_path}") # Log the path being used

        # Construct the command to run snpeff.py
        snpeff_py_command = [
            "python", "./pipes/snpeff.py",
            "--snp-indel-dir", snp_indel_dir, # Pass the directory path
            "--vcf-gz", gz,                   # Pass the gzipped filename
            "--vcf-unzipped", vcf_unzipped_abs_path, # <-- PASS ABSOLUTE PATH HERE
            "--out-prefix", out_snpeff,
            "--log-base-dir", to_dir,
            "--snpeff-base-dir", snpeff_base_dir # Pass the correct snpEff location
        ]

        # Log the command being run
        logger.info("Executing Command:\n" + ' '.join(snpeff_py_command))

        # Execute the snpeff.py script
        try:
            result_snpeff_py = subprocess.run(snpeff_py_command, check=True, capture_output=True, text=True)

            if result_snpeff_py.stdout:
                logger.info(f"snpeff.py stdout:\n{result_snpeff_py.stdout.strip()}")
            if result_snpeff_py.stderr:
                logger.warning(f"snpeff.py stderr:\n{result_snpeff_py.stderr.strip()}")
            logger.info("snpeff.py script executed successfully.")

            # Define the expected final output file (adjust if dbSNP step is added back)
            final_annotated_vcf = os.path.join(snp_indel_dir, f"{out_snpeff}.clinvar.vcf")
            if not os.path.exists(final_annotated_vcf):
                 logger.warning(f"Expected final annotated VCF not found after snpeff.py run: {final_annotated_vcf}")

        except subprocess.CalledProcessError as e:
            logger.error(f"Execution of snpeff.py failed with exit code {e.returncode}")
            if e.stdout: logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
            if e.stderr: logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
            # Decide if you want to exit the pipeline here
            # sys.exit(f"SnpEff step failed with exit code {e.returncode}")
        except Exception as e:
             logger.exception("An unexpected error occurred during the SnpEff step:")
             # Decide if you want to exit
             # sys.exit("Unexpected error in SnpEff step.")

    else:
        logger.error(f"snpEff.jar not found at {snpeff_jar_path}. Skipping SnpEff annotation.")
        # Decide if you want to exit or just continue
        # sys.exit("SnpEff dependency missing.")

    logger.info("%%%%%%%%%%%%%\n FINISHED SNP ANNOTATION \n%%%%%%%%%%%%%")
    logger.info('Next Step: SV-Calling via SNIFFLES2')




    #----------------------
    # Step 6: SV-Calling via SNIFFLES2
    
     # Define input BAM (Assuming haplotagged output from PMDV is desired)
    snp_indel_dir = os.path.join(to_dir, "snp_indel")
    # *** Check if this is the correct BAM. PMDV produces *.haplotagged.bam ***
    # If PMDV haplotagging step ran successfully:
    input_bam_sniffles = os.path.join(snp_indel_dir, f"{sample_name}_PMDV_wgs.haplotagged.bam")

    # Define output VCF path
    vcf_name_sniffles = os.path.join(to_dir, f"{sample_name}.sniffles.vcf")
    sample_id = sample_name # Redundant, just for clarity

    # Check if the selected input BAM exists
    if os.path.exists(input_bam_sniffles):
        logger.info(f"Running Sniffles2 on {sample_name} using sniffles2.py")

        # Define timeout (example: 6 hours = 21600 seconds)
        # Adjust as needed based on typical run times
        sniffles_timeout_seconds = 21600

        # Construct the command to run sniffles2.py
        sniffles2_py_command = [
            "python", "./pipes/sniffles2.py",
            "--input-bam", input_bam_sniffles,
            "--output-vcf", vcf_name_sniffles,
            "--ref", ref,
            "--sample-id", sample_id,
            "--log-base-dir", to_dir,
            "--threads", "16", # Or make this configurable
            "--timeout", str(sniffles_timeout_seconds) # Pass timeout to the script
        ]

        logger.info("Executing Command:\n" + ' '.join(sniffles2_py_command))

        try:
            # Run sniffles2.py. It handles its own timeout logging now.
            result_sniffles2_py = subprocess.run(sniffles2_py_command, check=True, capture_output=True, text=True)

            # Log output from sniffles2.py within pipeline.log for context
            if result_sniffles2_py.stdout:
                logger.info(f"sniffles2.py stdout:\n{result_sniffles2_py.stdout.strip()}")
            if result_sniffles2_py.stderr:
                logger.warning(f"sniffles2.py stderr:\n{result_sniffles2_py.stderr.strip()}")
            logger.info("sniffles2.py executed successfully.")

        except subprocess.CalledProcessError as e:
            # Check if the error code indicates a timeout occurred inside sniffles2.py (exit code 1 in the script)
            if e.returncode == 1 and "Sniffles2 command timed out" in e.stderr: # Check stderr for timeout message
                 logger.warning("Sniffles2 timed out (as reported by sniffles2.py). Proceeding with cleanup.")
                 # Don't necessarily exit pipeline.py here, allow cleanup.
            else:
                 logger.error(f"Execution of sniffles2.py failed with exit code {e.returncode}")
                 if e.stdout: logger.error(f"Failed command stdout:\n{e.stdout.strip()}")
                 if e.stderr: logger.error(f"Failed command stderr:\n{e.stderr.strip()}")
                 # Decide if you want to exit the pipeline here
                 # sys.exit(f"Sniffles2 step failed.")
        except Exception as e:
             logger.exception("An unexpected error occurred during the Sniffles2 step:")
             # Decide if you want to exit
             # sys.exit("Unexpected error in Sniffles2 step.")

    else:
        logger.error(f"Input BAM file for Sniffles2 not found: {input_bam_sniffles}. Skipping Sniffles2.")
        # Decide if you want to exit here
        # sys.exit("Input BAM for Sniffles2 not found.")

    logger.info("%%%%%%%%%%%%%\n FINISHED SNIFFLES2 (or timed out) \n%%%%%%%%%%%%%")

    # --- Keep the cleanup step ---
    logger.info("Attempting to clean up any lingering Sniffles processes...")
    # Use run with capture_output=True to avoid printing errors if no process is found
    cleanup_result = subprocess.run("pgrep -x sniffles | xargs kill -9", shell=True, capture_output=True, text=True)
    if cleanup_result.returncode == 0: # pgrep found something
        logger.info("Kill command sent to potential lingering Sniffles processes.")
    elif "kill: sending signal to XXX failed: No such process" in cleanup_result.stderr or cleanup_result.returncode != 0:
         logger.info("No lingering Sniffles processes found by pgrep or kill failed (may be expected if timeout didn't occur or process exited).")
    else:
         logger.warning(f"Sniffles cleanup command stderr: {cleanup_result.stderr.strip()}")

    logger.info("Pipeline Complete.") # Assuming this is the last step

if __name__ == '__main__':
   main()
