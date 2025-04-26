# Ryan Gallagher
# 03.14.2025
#
# This script will check what the most Knocked-out gene is in each sample
# regardless of target. This is a small QC analysis
#
#
# Look at the differences between TRAC and WT
#
# Separate TRAC and WT as controls and make individual plots
#
# DESeq2 WT vs. TRAC for inference on control.
#
# For Amy: Identify the regions targeted by CRISP-R
#          Where do the targeted guide RNA align?
#          *Note: Multiple sites (1. Where are they cutting 
#                                 2. How efficient was the targeting)
#                 Just need the genomic coordinates
#
#         BLAST API to find these

library(tidyverse)

cts = read.csv("./data/counts.csv", row.names = 1)
counts.df = as.data.frame(cts)
meta = read.csv("./data/meta.csv")

for (gene.loop in unique(meta$gene)) {
  
  quick = counts.df %>% select(starts_with(gene.loop) | starts_with("WT") | starts_with("TRAC"))
  quick.meta = meta %>% filter(gene == gene.loop)
  plt = unique(quick.meta$plate)
  
  quick.meta = meta %>% filter(gene == gene.loop | ((gene == "WT" | gene == "TRAC") & plate == plt))
  
  quick = counts.df %>% select(all_of(quick.meta$sample))
  
  quick$gene_mean = rowMeans(quick[,1:3])
  # quick$ctrl_mean = rowMeans(quick[,4:39]) # SELECT FOR CORRECT PLATE CONTROLS
  quick = quick %>% mutate(prop = gene_mean / ctrl_mean) %>% arrange(prop) %>% dplyr::slice(1:5)
}
