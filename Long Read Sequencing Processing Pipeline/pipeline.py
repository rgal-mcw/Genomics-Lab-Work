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

    #----------------------

    # Step 1: Prepare data for Transfer and send to Maple
    command = f"python ./pipes/file_transfer.py {sample_name} {data_dir} {to_dir} | tee {log_dir}file_transfer.log"

    result = subprocess.run(command, shell=True)

    # Exit if subprocess fails
    if not os.listdir(to_dir):
        logger.error("File transfer failed")
        sys.exit(1)
    
    logger.info("%%%%%%%%%%%%%\nFINISHED DATA TRANSFER \n%%%%%%%%%%%%%")
    #----------------------

    # Step 2: Read Alignment (FASTQs w/ GRCh38)
    if os.path.exists(to_dir + f"raw_data/{sample_name}.fastq.gz"):
        logger.info(f"Running alignment on {sample_name}")
        #result3 = subprocess.run(["python", "./pipes/alignment.py", sample_name, to_dir,ref])
        result3 = subprocess.run(f"python ./pipes/alignment.py {sample_name} {to_dir} {ref} | tee {log_dir}alignment.log", shell=True)
    else:
        logger.error("No FASTQ file found. Exiting script.")
        sys.exit(1)

    logger.info("%%%%%%%%%%%%%\n FINISHED ALIGNMENT \n%%%%%%%%%%%%%")
    logger.info('Next Step: SNP/INDEL via PMDV')

    #----------------------
    # Step 3: SNP/INDEL via PMDV

    if os.path.exists(to_dir + f"{sample_name}.sorted.bam"):
        logger.info(f"Running PMDV on {sample_name}")
        snp_indel_dir = to_dir + "snp_indel"
        bam_file = to_dir + f"{sample_name}.sorted.bam"
        print("Running Command:\n", ["./pipes/PMDV_pipeline.sh", sample_name, to_dir, snp_indel_dir, ref, bam_file])
        # Args: sample_name, to_dir, output_dir, ref_file, bam_file, chemistry (default 10.4.1)  
        result4 = subprocess.run(f"./pipes/PMDV_pipeline.sh {sample_name} {to_dir} {snp_indel_dir} {ref} {bam_file} | tee {log_dir}PMDV.log", shell=True)
        result_samtools = subprocess.run(f"samtools index -@ 8 '{snp_indel_dir}/{sample_name}_PMDV_wgs.haplotagged.bam'", shell=True)

    logger.info("%%%%%%%%%%%%%\n FINISHED PMDV SNP/INDEL CALLING \n%%%%%%%%%%%%%")
    logger.info('Next Step: CNV via Spectre')

    #----------------------
    # Step 4: Large SV via Spectre 

    logger.info("NOTE: A dedicated spectre conda environment needs to be created for this part of the pipeline. If this environment is not created yet, do so via instructions from the ont-spectre github documentation.")

    spectre_result = subprocess.run(f"/data/svi/prom/ont_pipeline/pipes/ONT_spectre.sh {to_dir}", shell=True)

    logger.info("%%%%%%%%%%%%%\n FINISHED SPECTRE \n%%%%%%%%%%%%%")
    logger.info('Next Step: SNP Annotation')

    #----------------------

    # Step 5: SNP Annotation on .phased.vcf file

    if os.path.exists(os.path.expanduser("~/snpEff/snpEff.jar")):
        logger.info(f"Running snpEff.jar on {sample_name}_PMDV_wgs.phased.vcf")
        gz = f"{sample_name}_PMDV_wgs.phased.vcf.gz"
        vcf_snpeff = f"{sample_name}_PMDV_wgs.phased.vcf"
        out_snpeff = f"{sample_name}_PMDV_wgs.phased.snpeff"
        snp_indel_dir = to_dir + "snp_indel"
        command_snpeff = f"./pipes/snpeff.sh {snp_indel_dir} {gz} {vcf_snpeff} {out_snpeff} | tee {log_dir}snpeff.log" 
        result_snpeff = subprocess.run(command_snpeff, shell=True)

        snpeff_vcf = snp_indel_dir + f"{sample_name}_PMDV_wgs.phased.snpeff.vcf"

    else:
        logger.error("No snpEff download found (~/snpEff). Skipping to next step.")
    
    logger.info("%%%%%%%%%%%%%\n FINISHED SNP ANNOTATION \n%%%%%%%%%%%%%")
    logger.info('Next Step: SV-Calling via SNIFFLES2')


    #----------------------
    # Step 6: SV-Calling via SNIFFLES2
    
    if os.path.exists(to_dir + f"{sample_name}.sorted.bam"):
        timeout = 360
        logger.info(f"Running Sniffles2 on {sample_name}") 

        input_bam = snp_indel_dir + f"/{sample_name}_PMDV_wgs.bam"
        vcf_name = to_dir + f"{sample_name}.sniffles.vcf"
        sample_id = sample_name
        
        try:
            command4 = f"./pipes/sniffles2.sh {input_bam} {ref} {vcf_name} {sample_id} | tee {log_dir}sniffles2.log"
            result4 = subprocess.run(command4, shell=True, timeout=timeout)
        except subprocess.TimeoutExpired:
            logger.info("Sniffles2 timed out.")
            logger.info("ATTENTION: As of 04/25/2024, Sniffles2 will not terminate on its own.\nIf the final message from Sniffles2 that you see before this says:\n'Wrote x called SVs to ~dir~/*.sniffles.vcf' - then the script has completed.\n\nIf you do not see the above message, then the Sniffles2 did NOT complete and TIMEDOUT before completion. If this is the case, manually run Sniffles2!\n\n----------\n----------")

    else: 
        logger.error("No sorted BAM file found. Exiting script.")
        sys.exit(1)

    logger.info("%%%%%%%%%%%%%\n FINISHED SNIFFLES2 \n%%%%%%%%%%%%%")
    subprocess.run("pgrep -x sniffles | xargs kill -9", shell=True)
    logger.info("Pipeline Complete.")

if __name__ == '__main__':
   main()
