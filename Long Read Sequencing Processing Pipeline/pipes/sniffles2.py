# Title: Sniffles2 for pipeline.py 
# Author: Ryan Gallagher, Broeckel Lab 2024
# Date: 04-22-2024
#
# Description: This script runs Sniffles2 on the SV-Calling step of the pipeline.
#              This script will run the command line tool Sniffles given the parsed
#              arguments from pipeline.py..
#
#
# Dev Notes: 
#           - This script should check for current version of Sniffles2
#

import os
import subprocess 
import argparse


def run_sniffles2(input_bam, ref, vcf_name, sample_id):
    command = [
        'sniffles',
        '--input', input_bam,
        '--reference', ref,
        '--vcf', vcf_name,
        '--snf', 'sniffles.snf',
        '--threads', '16',
        '--minsupport', 'auto',
        '--minsvlen', '35',
        '--minsvlen-screen-ratio', '0.9',
        '--mapq', '25',
        '--qc-stdev', 'True',
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
        '--detect-large-ins', 'True',
        '--cluster-binsize', '100',
        '--cluster-r', '2.5',
        '--cluster-repeat-h', '1.5',
        '--cluster-repeat-h-max', '1000',
        '--cluster-merge-pos', '150',
        '--cluster-merge-len', '0.33',
        '--cluster-merge-bnd', '1500',
        '--genotype-ploidy', '2',
        '--genotype-error', '0.05',
        '--sample-id', sample_id,
        '--allow-overwrite',
        '--symbolic',
        '--max-del-seq-len', '50000'
    ]

    cmd_string = ' '.join(command)
    print('Running Sniffles2 with the following command: ' + cmd_string)

    result = subprocess.run(command)


def main():
    parser = argparse.ArgumentParser(description='Process inputs')
    parser.add_argument('sample_name', type=str, help='Enter the sample name')
    parser.add_argument('to_dir', type=str, help='Enter the MAPLE directory where this should be stored')
    parser.add_argument('ref', type=str, help='Enter the reference genome')
    args = parser.parse_args()

    # Set variables
    input_bam = args.to_dir + args.sample_name + ".sorted.bam"
    ref = args.ref
    vcf_name = args.to_dir + args.sample_name + ".sniffles.vcf"
    sample_id = args.sample_name
   
    print('Starting Sniffles2')  
    run_sniffles2(input_bam, ref, vcf_name, sample_id)

if __name__ == "__main__":
    main()
