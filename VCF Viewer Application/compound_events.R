###
# Title: Find Compound Events
# Author: Ryan Gallagher, Broeckel Lab Bioinformatics (2024)
# 
# Description: This script takes Pepper-Margin-DeepVariant's phased
#              SNP/INDEL .vcf.gz output file and finds where there are 
#              different compounding events at the same position on
#              separate alleles. After finding those events, this script
#              outputs a position.txt file and subsets the original .vcf
#              for only those calls. That file is then converted to a .tsv
#              and displayed in our shiny app (Toucan). 
#
# To Be Added: 
#   - Need to find where different SNP calls overlap on sep allales
#   - Need to find where INSERTION + DEL occur at same pos. at diff alleles

# Libraries ----
print("Loading Packages")
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tools))
### ----

# Retrieve file path argument ----
args = commandArgs(trailingOnly = TRUE)

# Check if the required argument is provided
if (length(args) == 0) {
  stop("No arguments supplied. Please provide the path to the .phased.vcf.gz file.")
}

# Assign the first argument to the 'file' variable
file = file_path_as_absolute(args[1])
dir = dirname(file)
### ----

### Find Compound Het & Write positions.txt -----
print('Reading VCF')
vcf = readVcf(file, 'hg38')
print('Finding Compound Hets & Building positions.txt')
is_snp = isSNV(vcf)
is_indel = isIndel(vcf)

gt = geno(vcf)$GT
sample_names = colnames(gt)
sample_name = sample_names[1]  # Adjust if multiple samples

gr = rowRanges(vcf)
snp_gr = gr[is_snp]
indel_gr = gr[is_indel]

# Find overlaps between SNPs and indels
overlaps = findOverlaps(snp_gr, indel_gr)

# Get indices of overlapping SNPs and indels
snp_overlapping_idx = queryHits(overlaps)
indel_overlapping_idx = subjectHits(overlaps)

# Extract genotype strings for overlapping SNPs and indels
snp_gt = gt[is_snp, sample_name][snp_overlapping_idx]
indel_gt = gt[is_indel, sample_name][indel_overlapping_idx]

# Function to parse genotype strings
parse_genotype = function(gt_str) {
  phased = grepl("\\|", gt_str)
  alleles = strsplit(gt_str, "[/|]")
  allele1 = sapply(alleles, "[", 1)
  allele2 = sapply(alleles, "[", 2)
  het = allele1 != allele2
  return(data.frame(allele1 = allele1, allele2 = allele2, phased = phased, het = het, stringsAsFactors = FALSE))
}

# Parse the genotypes
snp_gt_parsed = parse_genotype(snp_gt)
indel_gt_parsed = parse_genotype(indel_gt)

df = data.frame(
  snp_idx = snp_overlapping_idx,
  indel_idx = indel_overlapping_idx,
  snp_allele1 = snp_gt_parsed$allele1,
  snp_allele2 = snp_gt_parsed$allele2,
  snp_phased = snp_gt_parsed$phased,
  snp_het = snp_gt_parsed$het,
  indel_allele1 = indel_gt_parsed$allele1,
  indel_allele2 = indel_gt_parsed$allele2,
  indel_phased = indel_gt_parsed$phased,
  indel_het = indel_gt_parsed$het,
  stringsAsFactors = FALSE
)

# Keep only rows where both SNP and indel are phased and heterozygous
df = df[df$snp_phased & df$indel_phased & df$snp_het & df$indel_het, ]

# Define the variant alleles based on allele codes ('0' for reference, '1' for alternative)
df$snp_variant_on_hap1 = df$snp_allele1 == '1'
df$snp_variant_on_hap2 = df$snp_allele2 == '1'
df$indel_variant_on_hap1 = df$indel_allele1 == '1'
df$indel_variant_on_hap2 = df$indel_allele2 == '1'

# Identify where variant alleles are on opposite haplotypes
df$opposite_phase = (df$snp_variant_on_hap1 != df$indel_variant_on_hap1) & 
  (df$snp_variant_on_hap2 != df$indel_variant_on_hap2)

# Extract the rows where the SNP and indel are on opposite phases
opposite_phase_df = df[df$opposite_phase, ]

# Optionally, you can extract more information or write to a file
# For example, get the positions and chromosomes
snp_positions = snp_gr[opposite_phase_df$snp_idx]
indel_positions = indel_gr[opposite_phase_df$indel_idx]

# Combine positional information with the data frame
opposite_phase_df$snp_chr = as.character(seqnames(snp_positions))
opposite_phase_df$snp_pos = start(snp_positions)
opposite_phase_df$indel_chr = as.character(seqnames(indel_positions))
opposite_phase_df$indel_pos = start(indel_positions)

# Save the opposite_phase_df to a CSV file
#write.csv(opposite_phase_df, "./svi-0029/opposite_phase_variants.csv", row.names = FALSE)

# Extract the CHROM and POS columns
positions = data.frame(
  CHROM = opposite_phase_df$snp_chr,
  POS = opposite_phase_df$snp_pos
)

print("Writing positions.txt")
# Save the positions to a tab-delimited text file without headers
write.table(positions, paste(dir, "/positions.txt", sep=''), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# View the updated data frame with positional information
#print(opposite_phase_df)

### ----

### BGZip the original .vcf ----
print("Running bcftools Commands")
# BGZip -> Index -> Subset -> ReIndex -> Write new .tsv
#system(paste(paste('bgzip -c', paste(file, '>', sep=' ')), paste(file, ".gz", sep='')))
filename_no_ext <- sub("\\.vcf\\.gz$", "", file)
file.ann = paste(filename_no_ext, '.snpeff.rsID.clinvar.vcf', sep='')
system(paste(paste('bgzip -c', paste(file.ann, '>', sep=' ')), paste(file.ann, ".gz", sep='')))
system(paste("bcftools index", paste(file.ann,".gz", sep='')))
system(paste('bcftools view -R', paste(dir, '/positions.txt', sep=''), paste(file.ann, '.gz', sep=''), '-Oz -o', paste(dir, 'compound_hets.vcf.gz', sep='/')))
system(paste('bcftools index', paste(dir, 'compound_hets.vcf.gz', sep='/')))

# Construct the command using sprintf()
print("Building compound_hets.tsv")
query.cmd <- sprintf("(
  echo -e \"CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tGT\"
  bcftools query -f '%%CHROM\\t%%POS\\t%%ID\\t%%REF\\t%%ALT\\t%%QUAL\\t%%FILTER\\t%%INFO\\t[%%GT]\\n' %s
) > %s", paste(dir, 'compound_hets.vcf.gz', sep='/'), paste(dir, 'compound_hets.tsv', sep='/'))
system(query.cmd)
### ----
