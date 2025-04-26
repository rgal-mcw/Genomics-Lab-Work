#!/bin/bash
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
#


# Assign arguments to variables for better readability
input_bam=$1
ref=$2
vcf_name=$3
sample_id=$4

# Construct the command with the passed arguments
sniffles --input "$input_bam" \
         --reference "$ref" \
         --vcf "$vcf_name" \
         --threads 16 \
         --minsupport auto \
         --minsvlen 35 \
         --minsvlen-screen-ratio 0.9 \
         --mapq 25 \
         --qc-stdev True \
         --qc-stdev-abs-max 500 \
         --qc-coverage 1 \
         --long-ins-length 2500 \
         --long-del-length 50000 \
         --long-del-coverage 0.66 \
         --long-dup-length 50000 \
         --long-dup-coverage 1.33 \
         --max-splits-kb 0.1 \
         --max-splits-base 3 \
         --min-alignment-length 1000 \
         --phase-conflict-threshold 0.1 \
         --detect-large-ins True \
         --cluster-binsize 100 \
         --cluster-r 2.5 \
         --cluster-repeat-h 1.5 \
         --cluster-repeat-h-max 1000 \
         --cluster-merge-pos 150 \
         --cluster-merge-len 0.33 \
         --cluster-merge-bnd 1500 \
         --genotype-ploidy 2 \
         --genotype-error 0.05 \
         --sample-id "$sample_id" \
         --allow-overwrite \
         --symbolic \
         --max-del-seq-len 50000
