# ONT PromethION Sequencing Data Processing Pipeline

**Broeckel Lab Bioinformatics (2024-2025)**

## Overview

This repository contains a Python-based pipeline designed to process and analyze Oxford Nanopore Technologies (ONT) PromethION sequencing data. The pipeline automates several steps, starting from transferring raw data from the PromethION server to performing alignment, variant calling (SNPs, INDELs, SVs), and annotation.

## Pipeline Steps

The pipeline executes the following sequence of bioinformatics tasks:

1.  **Data Preparation & Transfer:**
    * Connects to the PromethION server via SSH (requires key-based authentication).
    * Concatenates FASTQ files from `fastq_pass/` on the server.
    * Moves report files (`.txt`, `.csv`, `.html`, etc.) on the server.
    * Transfers the concatenated FASTQ, reports, and `pod5` data to the analysis server (Maple) using `rsync`.
2.  **Read Alignment:**
    * Aligns the transferred FASTQ reads against the GRCh38 reference genome using `minimap2`.
    * Sorts and indexes the resulting BAM file using `samtools`.
3.  **SNP/INDEL Calling:**
    * Calls Single Nucleotide Polymorphisms (SNPs) and small Insertions/Deletions (INDELs) using the Pepper-Margin-DeepVariant (PMDV) Docker container.
    * Produces a phased VCF file and a haplotagged BAM file.
    * Indexes the haplotagged BAM file.
4.  **Large Structural Variant (SV) Calling (CNVs):**
    * Calculates read depth using `mosdepth` (requires 'spectre' conda environment).
    * Calls Copy Number Variants (CNVs) using `spectre` (requires 'spectre' conda environment).
5.  **SNP Annotation:**
    * Annotates the phased VCF from PMDV using `snpEff` and `SnpSift`.
    * Adds annotations from ClinVar.
    * *(Note: dbSNP annotation was present in older `.sh` scripts but appears removed/commented out in `snpeff.py`)*
    * *(Note: An optional `filter_snpeff_vcf.py` script exists but is not called by the main `pipeline.py`. It needs to be run separately)*
6.  **Structural Variant (SV) Calling (General):**
    * Calls various types of structural variants using `Sniffles2` on the haplotagged BAM from PMDV.
    * Includes a timeout mechanism and attempts cleanup of lingering processes.

## Prerequisites and Setup

Before running the pipeline, ensure the following are set up:

1.  **Software Dependencies:**
    * Python 3.x
    * `minimap2`
    * `samtools`
    * `docker` (or `podman` - may require modification for PDMV step)
    * `rsync`
    * `conda` (for managing environments)
    * `java` (for snpEff/SnpSift)
    * `mosdepth` (installed via conda, likely in 'spectre' environment)
    * `spectre` (installed via conda, in 'spectre' environment)
    * `Sniffles2` (installed, likely via conda)
    * `snpEff` & `SnpSift` (Download required, ensure `snpEff.jar` and `SnpSift.jar` are accessible. Default path expected: `/data/ref/snpEff/`)
2.  **Docker Images:**
    * Pepper-Margin-DeepVariant: `kishwars/pepper_deepvariant:r0.8` (or compatible version)
3.  **Conda Environments:**
    * A dedicated conda environment named `spectre` is required for running Mosdepth and Spectre. Activate this environment before running the main pipeline.
4.  **Reference Files:**
    * Reference Genome: GRCh38 FASTA (Default expected: `Homo_sapiens_assembly38_noALT_noHLA_noDecoy_ERCC.fasta`)
    * ClinVar VCF: (Default expected: `clinvar.vcf.gz`)
5.  **SSH Access:**
    * SSH key-based authentication must be set up from the analysis server (Maple) to the PromethION server. Ensure your public key is in the `authorized_keys` file on the PromethION server for the `prom` user.
6.  **Directory Structure:**
    * The pipeline expects a specific output directory structure on the analysis server (Maple). Create the main output directory before running. Subdirectories (`raw_data/`, `logs/`, `snp_indel/`, `spectre/`) will be created by the scripts.

## Usage

The main pipeline is executed using `pipeline.py`.

**Command:**

```bash
# Ensure the correct conda environment (e.g., spectre) is activated first.
# source activate spectre

python pipeline.py <sample_name> <promethion_data_dir> <maple_output_dir> [--log-level LEVEL]
