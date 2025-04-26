#
#/bin/bash
## Run Spectre 



# =--=--=--=--=--= CONSOLE CHECKS =--=--=--=--=--=
# Sample to Run
if [ -z "$1" ]; then
  echo -e "Usage: ./run_spectre <directory path> \n \n NOTE: This should be the full path to the sample directory (ending with flowcell_10.4.1/--SAMPLE--)\n \n"
	exit 1
fi

DIR_PATH=$1
SAMPLE=$(basename "$DIR_PATH")

if [ ! -d "$DIR_PATH" ]; then
	echo -e "ERROR: Directory $DIR_PATH does not exist. \n \n Please check the path you've entered. \n \n"
	exit 1
else
	if [ ! -d "$DIR_PATH/spectre" ]; then
		mkdir -p "$DIR_PATH/spectre"
	fi
fi 

# ---  ---  ---  ---  ---  ---  ---  ---  --- 


# =--=--=--=--=--= CONDA CHECK =--=--=--=--=--=
# Check if the conda environment is activated
# Desired conda environment name
DESIRED_ENV="spectre"
if [[ "$CONDA_DEFAULT_ENV" != "$DESIRED_ENV" ]]; then
    echo "CONDA: Activating conda environment: $DESIRED_ENV"
    # Activate the desired conda environment
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate $DESIRED_ENV
else
    echo "CONDA: Conda environment $DESIRED_ENV is already activated."
fi
# ---  ---  ---  ---  ---  ---  ---  ---  --- 


# =--=--=--=--=--= MOSDEPTH CHECK =--=--=--=--=--=
## Check if mosdepth stuff exists
if [ -f "${DIR_PATH}/spectre/${SAMPLE}.regions.bed.gz" ]; then
	echo "MOSDEPTH already ran - starting spectre"
else
	echo "MOSDEPTH not ran"
	echo "MOSDEPTH: Beginning"
	mosdepth -t 8 -x -b 1000 -Q 20 "${DIR_PATH}/spectre/${SAMPLE}" "${DIR_PATH}/${SAMPLE}.sorted.bam"
	echo "MOSDEPTH: Finished"
fi

# ---  ---  ---  ---  ---  ---  ---  ---  --- 

# =--=--=--=--=--= POINT TO PMDV SNV VCF =--=--=--=--=--=
if [ -d "${DIR_PATH}/snp_indel/" ]; then
	echo "SNP: PMDV SNP file Detected."
else
	echo "NO snp_indel/ DIRECTORY FOUND."
	echo -e "\n\nPlease run the ONT pipeline first.\n"
	exit 1
fi


# ---  ---  ---  ---  ---  ---  ---  ---  --- 

# =--=--=--=--=--= SPECTRE =--=--=--=--=--=
MOSDEPTH="${DIR_PATH}/spectre/"

echo "SPECTRE: Beginning"
spectre CNVCaller \
  --bin-size 1000 \
  --coverage $MOSDEPTH \
  --sample-id $SAMPLE \
  --output-dir "${DIR_PATH}/spectre/" \
  --reference /data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta \
  --snv "${DIR_PATH}/snp_indel/${SAMPLE}_PMDV_wgs.phased.vcf.gz" \

echo "SPECTRE: Finished"
