#!/bin/bash

BAM='/home/rgallagher/test/mcw_svi_9999.sorted.bam'
REF= '/data/ref/NCBI/GRCh38/Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta'
VCF='/home/rgallagher/test/mcw_svi_9999.sniffles.vcf'
SAMPLE_ID='mcw_svi_9999'


sniffles \
	--input $BAM \ 
	--reference $REF \
	--vcf $VCF \
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
  	--sample-id $SAMPLE_ID \
	--allow-overwrite \
	--symbolic \
	--max-del-seq-len 50000 
