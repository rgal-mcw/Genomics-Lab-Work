#!/bin/bash

if [ -z "$1" ]; then
	echo "No argument provided - please provide <path/to/sample_dir>"
	exit 1
fi

INPUT_DIR=$1
sample=$(basename "$INPUT_DIR")

# Check if the conda command is available
if ! command -v conda &> /dev/null
then
    echo "conda could not be found. Please ensure conda is installed and accessible."
    exit 1
fi

# Check if the clair3 environment is active
if [[ "$CONDA_DEFAULT_ENV" != "clair3" ]]; then
    echo "The clair3 environment is not active. Activating it now..."
    # Activate the clair3 environment
    source activate clair3
    
    # Check if activation was successful
    if [[ "$CONDA_DEFAULT_ENV" != "clair3" ]]; then
        echo "Failed to activate clair3 environment."
        exit 1
    else
        echo "clair3 environment activated successfully."
    fi
else
    echo "The clair3 environment is already active."
fi

# Add your commands that need to run within the clair3 environment below
echo "Running your commands in the clair3 environment..."
# Example command


# -- FROM CLAIR3 GITHUB
# make sure channels are added in conda
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge

# create conda environment named "clair3"
# replace clair3 by clair3-illumina for using illumina data
#conda create -n clair3 -c bioconda clair3 python=3.9.0 -y
#conda activate clair3

# run clair3 like this afterward
MODEL_NAME="r1041_e82_400bps_sup_v420"         # e.g. r941_prom_hac_g360+g422
MODEL_DIR="/home/rgallagher/rerio/clair3_models"

if [ -e "$INPUT_DIR" ]; then
	echo "Input file path: $INPUT_DIR"

	OUTPUT_DIR="$INPUT_DUR/snp_indel/clair3"
	REF="/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta"
	THREADS=32

	run_clair3.sh \
	  --bam_fn="${INPUT_DIR}/${sample}.sorted.bam" \
	  --ref_fn="${REF}" \
	  --threads="${THREADS}" \
	  --platform="ont" \
	  --model_path="${MODEL_DIR}/${MODEL_NAME}" \
	  --output="${OUTPUT_DIR}" \
	  --use_whatshap_for_final_output_haplotagging \
	  --enable_phasing \

fi
# Optionally deactivate the environment if desired
# source deactivate
