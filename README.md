# Genomics Lab Work

This repository shows a select amount of work that I've accomplished as a Biostatistician in the Medical College of Wisconsin / Childeren's Wisconsin Pediatric Genomics Lab. My work in this lab focuses on processing and analyzing Oxford Nanopore Long Read Sequencing through a curated and optimized, in-house pipeline written in Python. I also conduct statistical analyses to aid in processes in the clincal lab space as well as in our various research projects. 

The directories in this repository feature work from these projects:

**Long Read Sequencing Processing Pipeline:**

> This directory show my code for our ONT Long-Read Sequencing Processing Pipeline that I created and implemented to (so far) process >200TB worth of sequencing data for our lab. It also features software that annotates and filters calls based on functional prediction or rarity. This pipeline has logging, data management, and modularity capabilities. 

**Batch Effect Report:**

> A document outlining my research on batch effect as it effects Differential Expression Analyses. I code through an explicit application of a function `limma::removeBatchEffect()` , and describe how it functions. I then compare this method to other potential solutions. This is done in RMarkdown.

**Clinical Statistics Reports:**

> Our clinical lab requires Quality Control reporting to assure calibration of instruments. The lab asked me if I could apply my expertise to help diagnose whether the different machines were well calibrated. There are three select reports in this directory that feature my application of Mixed Models, Hypothesis Testing, and Principal Component Analysis.

**Differential Expression Analysis:**

> This directory features R code and figures on some of the Differential Expression Analysis processes I did for one of our project. The plots in this figure show Knock-out efficiency and data availability from our experiements.

**VCF Viewer Application:**

> This code shows my attempt to code an RShiny app that could be used as a VCF viewer / filtering interface. This application intakes an annotated .VCF file, and allows the user to sift through the data in a tabular format. It allows for multiple filtering and exporting to .csv. 

