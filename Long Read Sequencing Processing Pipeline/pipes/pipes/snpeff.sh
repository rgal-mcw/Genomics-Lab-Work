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

# Set variables
  #export DIR="/data/svi/prom/mcw_svi/flowcell_9.4.1/${subdir}/LongRead_SNP/PMDV/wgs_phased"
  #export GZ="${subdir}_PMDV_wgs.phased.vcf.gz"
  #export VCF="${subdir}_PMDV_wgs.phased.vcf"
  #export OUT="${subdir}_PMDV_wg.phased.snpeff"
  DBSNP="/data/ref/dbSNP/00-All.vcf.gz"
  CLINVAR="/data/ref/ClinVar/clinvar.vcf.gz"
  file_mod_date=$(stat -c %Y ${CLINVAR})

  # Unzip the vcf
  if [ ! -f ${DIR}/${VCF} ]; then
    echo "Unzipping ${GZ}..."
    gunzip ${DIR}/${GZ}
  fi

  DIR=$1
  GZ=$2
  VCF=$3
  OUT=$4

if [ -f /data/ref/ClinVar/update.log ]; then
    echo "Last updated: $(head -n 1 /data/ref/ClinVar/update.log)"
    echo "Use ~~ /data/ref/ClinVar/update_clinvar.sh ~~ to update "
else
    echo "No previous updates found."
fi


  # Add annotations (BETA: How do we know that we have the most up-to-date version of the Ensembl database?)
  #
  #version=$(curl -q http://ftp.ensembl.org/pub/VERSION)
  #echo "Running snpEff on ${VCF} with Ensembl version ${version}."
  java -jar /home/rgallagher/snpEff/snpEff.jar GRCh38.105 ${DIR}/${VCF} > ${DIR}/${OUT}.vcf

  # Add Clinvar annotations
  java -jar /home/rgallagher/snpEff/SnpSift.jar annotate ${CLINVAR} ${DIR}/${OUT}.vcf > ${DIR}/${OUT}.rsID.clinvar.vcf 

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



