# Load libraries
library(dplyr)
library(biomaRt)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(AnnotationDbi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)  # hg38 annotations
library(org.Hs.eg.db)



setwd("~/Desktop/Project/maple/shiny_app/toucan/")
#file = "~/Desktop/Project/maple/svi-0029/mcw_svi_0029_vcf_calls.tsv"
file = "~/Desktop/Project/maple/svi-0031/test3_vcf_calls.tsv"
sample_name = paste((strsplit(basename(file), "_")[[1]])[1:3], collapse = "_")
snf = "~/Desktop/Project/maple/svi-0029/mcw_svi_0029.haplotagged.sniffles.vcf"

# ATTACH GENE BY CHROM+POS ON HG38 -----
snf_df = read.table(snf, header=F, fill=T) %>% dplyr::rename(CHROM=V1, POS=V2, ID=V3, REF=V4, ALT=V5, QUAL=V6, FILTER=V7) %>%
  separate(V10, into = c("GT", "GQ", "DR", "DV", "PS"), sep = ":") %>%
  # Extract the first entry as the Precision value
  mutate(Precision = str_extract(V8, "^[^;]+")) %>%
  # Remove the first entry from the V8 string
  mutate(V8 = str_remove(V8, "^[^;]+;")) %>%
  # Split the string by semicolons, and handle entries without '='
  separate_rows(V8, sep = ";") %>%
  mutate(V8 = ifelse(str_detect(V8, "="), V8, paste0(V8, "=", V8))) %>%
  separate(V8, into = c("key", "value"), sep = "=", fill = "right") %>%
  pivot_wider(names_from = key, values_from = value)  # Convert to wide format (columns)

snf_df = snf_df[!is.na(snf_df$POS) & !is.na(snf_df$END), ]

# Create GRanges object for your variants
snf_gr = GRanges(
  seqnames = snf_df$CHROM,
  ranges = IRanges(start = as.numeric(snf_df$POS), end = as.numeric(snf_df$END)),
  strand = "*"
)

# Load hg38 gene annotations
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
genes = genes(txdb)
# Find overlaps between variants and genes
overlaps = findOverlaps(snf_gr, genes)

# Extract overlapping gene IDs
overlap_df = as.data.frame(overlaps)
gene_ids = mcols(genes)$gene_id[subjectHits(overlaps)]

# Map Entrez IDs to gene symbols
gene_symbols = mapIds(
  org.Hs.eg.db,
  keys = as.character(gene_ids),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Create a data frame with variant indices and gene symbols
overlap_df$gene_symbol = gene_symbols

# Handle multiple overlaps by collapsing gene symbols

gene_overlap_summary = overlap_df %>%
  group_by(queryHits) %>%
  summarize(GENE_OVERLAP = paste(unique(gene_symbol), collapse = ";"))

# Initialize GENE_OVERLAP column in your DataFrame
snf_df$GENE_OVERLAP = NA

# Add gene overlap information to your DataFrame
snf_df$GENE_OVERLAP[gene_overlap_summary$queryHits] = gene_overlap_summary$GENE_OVERLAP

# ------

vcf = read.table(file, sep="\t", header=TRUE, fill=TRUE)

bby_vcf = vcf %>% filter(CLNVC == "single_nucleotide_variant")


compound = bby_vcf %>%
inner_join(snf_df, by = 'CHROM', relationship = "many-to-many") %>%
filter(as.numeric(POS.x) >= as.numeric(POS.y) & as.numeric(POS.x) <= as.numeric(END)) %>%
#filter(GT.x != "1/1", GT.y != "1/1", FILTER.x == "PASS", FILTER.y == "PASS") %>%
dplyr::select(CHROM, POS.x, POS.y, END, FILTER.y, GT.x, GT.y, ID.y, GENEINFO, GENE_OVERLAP, everything())

