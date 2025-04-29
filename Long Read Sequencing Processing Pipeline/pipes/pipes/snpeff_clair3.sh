#!/bin/bash

###################
### Title: snpEff  🦠 🛠
### Author: Ryan Gallagher 🧙 
### Date: 03-25-2024 🍀
####################

# This script relies on setup from ~/snpEff/README.txt, so go do that now.
#  -> If you don't have that file, look in /home/ryangallagher/snpEff/README.txt
#
# This script will run snpEff on a VCF file and output a new VCF file with 
# the snpEff annotations. It will then use the downloaded dbSNP file to 
# assign rsID's to SNPs
#
# The CLINVAR database file updates bi-weekly. This script will reference 
# the ftp site and download the latest version of the database to the 
# reference directory.
#
# The Ensembl annotation database is automatically updated per each version
# release when this script runs.
#
# From what I understand, dbSNP does not need to be updated.
#

# Set to run in all subdirectories of flowcell 9.4.1
current_date=$(date +%s)

if [ -z "$1" ]; then
	echo -e "No argument provided\nInput should be ./snpeff_clair.sh <path/to/sample>"
	exit 1
fi

DIR_INIT=$1
sample=$(basename "$DIR")

# Set variables
  #export DIR="/data/svi/prom/mcw_svi/flowcell_9.4.1/${subdir}/LongRead_SNP/PMDV/wgs_phased"
  #export GZ="${subdir}_PMDV_wgs.phased.vcf.gz"
  #export VCF="${subdir}_PMDV_wgs.phased.vcf"
  #export OUT="${subdir}_PMDV_wg.phased.snpeff"
  DBSNP="/data/ref/dbSNP/00-All.vcf.gz"
  CLINVAR="/data/ref/ClinVar/clinvar.vcf.gz"
  file_mod_date=$(stat -c %Y ${CLINVAR})

  DIR="${DIR_INIT}/snp_indel/clair3"
  GZ="phased_merge_output.vcf.gz"
  VCF="phased_merge_output.vcf"
  OUT="${sample}_clair3.snpeff"


  # Add annotations (BETA: How do we know that we have the most up-to-date version of the Ensembl database?)
  #
  #version=$(curl -q http://ftp.ensembl.org/pub/VERSION)
  echo "Running snpEff on ${VCF}."
  java -jar /home/rgallagher/snpEff/snpEff.jar GRCh38.105 ${DIR}/${VCF} > ${DIR}/${OUT}.vcf

  # Add rsID's via dbSNP
  java -jar /home/rgallagher/snpEff/SnpSift.jar annotate ${DBSNP} ${DIR}/${OUT}.vcf > ${DIR}/${OUT}.rsID.vcf

  # Check for ClinVar updates
 two_weeks_ago=$(date -d "2 weeks ago" +%s)
 if [[ $file_mod_date -ge $two_weeks_ago ]]
 then 
   echo "Updating ClinVar database..."
   curl -o ${CLINVAR} https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz  
   curl -o ${CLINVAR}.tbi https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
   curl -o ${CLINVAR}.md5 https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.md5
   echo "ClinVar database updated."
 else 
   echo "ClinVar database is up-to-date."
 fi


  # Add Clinvar annotations
  java -jar /home/rgallagher/snpEff/SnpSift.jar annotate ${CLINVAR} ${DIR}/${OUT}.rsID.vcf > ${DIR}/${OUT}.rsID.clinvar.vcf 

  # This script runs in pwd - move items to $DIR
  mv snpEff_summary.html ${DIR}/snpEff_summary.html
  mv snpEff_genes.txt ${DIR}/snpEff_genes.txt

  echo "Finished running snpEff on ${VCF}."


## Filter for pathogenic or unknown variants. Not interested in benign.
#### i.e. clinically relavent

## This will need to be figured out with Uli - there are likely relavent fields and non-relavent fields (see annotations.csv on local macbook (Ryan))
## As of now, with how large these annotations are, it's impossible to open the annotated .vcf files in IGV
## so we need to either reduce the amount of calls, or reduce the amount of annotations in each call.

## There are filter options in snpSift, or we can grep stuff.
## i.e. `grep ";OM;" data.txt` would grab all lines with the OMIM identifier. It's a way of filtering. 



