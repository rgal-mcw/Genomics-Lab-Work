#!/bin/bash

# Index the original VCF file if it's not already indexed
bcftools index ../../svi-0029/mcw_svi_0029_PMDV_wgs.phased.snpeff.rsID.clinvar.vcf.gz

# Subset the VCF file using the positions file
bcftools view -R ../../svi-0029/positions.txt ../../svi-0029/mcw_svi_0029_PMDV_wgs.phased.snpeff.rsiD.clinvar.vcf.gz -Oz -o ../../svi-0029/compound_hets.vcf.gz

# Index the new VCF file
bcftools index ../../svi-0029/compound_hets.vcf.gz


(
echo -e "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tGT"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t[%GT]\n' ../../svi-0029/compound_hets.vcf.gz
) > ../../svi-0029/compound_hets.tsv
