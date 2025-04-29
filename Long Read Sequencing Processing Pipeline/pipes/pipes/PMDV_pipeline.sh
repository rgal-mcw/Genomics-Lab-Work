#!/bin/bash
# --------------------------------------------------------------------------------
# Script Name: pepper🌶️_margin📈_deepvariant🧬🐳.sh
# Descripton: Implementation template for Pepper-Margin-DeepVariant SNP/INDEL caller for 
#             ONT Longread data via Docker image. This version is for pipeline integration. 
# GitHub: https://github.com/kishwarshafin/pepper/
# ================================================================================
# Author: Ryan Gallagher @ Broeckel Lab
# Date: 02-12-2023
# Usage: 
#  *Required Parameters*
#     run_pepper_margin_deepvariant call_variant \
#         -b, --bam             Path to a bam file containing mapping between reads and a reference.
#         -f, --fasta           Path to a reference file in FASTA format.
#         -o, --output_dir      Path to a output directory.
#         -t, --threads         Number of threads to use.
#         --ont_r9_guppy5_sup   Set to call variants from R9.4.1 Guppy 5+ sup Oxford Nanopore Reads.
#         --ont_r10_q20         Set to call variants from R10.4 Q20 Oxford Nanopore reads.
#
#  *Optional Parameters*
#         -r, --region          Region in [contig_name:start-end] format
#         -p, --output_prefix   Prefix for output filename. Do not include extension in prefix.
#
#   TONS OF OPTIONAL PARAMS FOR PEPPER / MARGIN / DV.
#   (https://github.com/kishwarshafin/pepper/blob/r0.8/docs/usage/usage_and_parameters.md)
# --------------------------------------------------------------------------------

#sample="23-0385"
SAMPLE="$1"
#NWGC_id="NA"

## Set Directories
INPUT_DIR="$2"
OUTPUT_DIR="$3"
## Set Reference Directory
#export REF_DIR="/data/ref/NCBI/GRCh38"
REF="${4:-/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta}"
REF_DIR="${REF%/*}"
## Set File Names (check that this is the actual name)
#export BAM="${sample}.sorted.bam"
#export BAM_i="${sample}.sorted.bam.bai"
BAM="$5"

#export HGREF='Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta'
#export HGREF_i='Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta.fai'

## Set Options 
CHEM="${6:-ont_r10_q20}"
THREADS=64 #(arbitrary value. Set to high for speed??)
REGION=''
# -r ${REGION} \  <-- add this to options if region is set
mkdir "${INPUT_DIR}snp_indel"


PREFIX="${SAMPLE}_PMDV_wgs"
# !Change chemistry option when necessary! #
echo "Sample: $SAMPLE"
echo "Input Directory: $INPUT_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo "Reference: $REF"
echo "Reference Directory: $REF_DIR"
echo "BAM: $BAM"
echo "CHEM: $CHEM"
#----------------------
## Run Pepper-Margin-DeepVariant (with phased output on 10.4 chemistry)
##                               (with minimum quality settings)
#----------------------

## Initialize Docker Container & Run Command 
docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${REF_DIR}":"${REF_DIR}" \
-u $(id -u):$(id -g) \
kishwars/pepper_deepvariant:r0.8 \
run_pepper_margin_deepvariant call_variant \
-b "${BAM}" \
-f "${REF}" \
-o "${OUTPUT_DIR}" \
-p "${PREFIX}" \
-t ${THREADS} \
--pepper_min_mapq 15 \
--pepper_min_snp_baseq 15 \
--pepper_min_indel_baseq 15 \
--dv_min_mapping_quality 15 \
--dv_min_base_quality 15 \
--phased_output \
--${CHEM}

