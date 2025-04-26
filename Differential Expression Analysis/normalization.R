# Ryan Gallagher
# 03.12.2025
#
#
# This script will use DESeq2 to normalize the data. We will write that normalized
# data and output to send.
#
# Normalization method copied from https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
# With justification from: https://www.youtube.com/watch?v=UFB993xufUU
# 
# NOTES:
#
# I had mulled over the idea for whether we would calculate our size factors using all samples or within subsets (replicates).
# Using all samples would consider all experiments together (controls & each Knock-Downs with replicates) when calculating 
# our row medians, etc.. At first this seemed problematic. 
#
# At first, I thought it would make sense to normalize by separating our replicates into subsets,
# but when I did more reading, I found that the DESeq2 medians-of-ratios method is built to handle both:
#   1) Differences in library size (total counts)
#   2) Differences in library composition (where counts are distributed)
#
# The youtube link (above) describes this method nicely. When estimating size factors, there's more significance given to
# genes where counts are *more in common across samples* than those of outliers (like in K-O situations) within small groups of samples.
# Because of this, I now believe that there shouldn't be a problem with estimating size factors across all samples.
#
# I believe this primarily because we have so many genes within our data, and that the small amount of genes that are effected by the
# knock-down will have little impact on our size factors.
#
# For this reason, I think this code and it's output accurately normalize the data. 


library(tidyverse)
library(DESeq2)

# Read in combined counts data w. meta
cts = as.matrix(read.csv("./data/counts.csv", row.names = 1)) #colbind() of all our combined_counts.txt files.
meta = read.csv("./data/meta.csv") # Info gathered from sample name

#Check that the order is the same (requirement of DESeq2)
all(meta$sample == colnames(cts))

## Normalize
### Make meta to distinguish WT / TRAC across plates
meta.norm = meta %>%
  mutate(gene = ifelse(condition == "control", paste(gene, plate, sep = "_"), gene))

cts.norm = cts
all(meta.norm$sample == colnames(cts.norm))

# Make DESeq2 object
dds = DESeqDataSetFromMatrix(countData = cts.norm, colData = meta.norm,
                                  design = ~gene) #Design is anticipated in function by not used in size factor estimation

## Check if our DESeq object is the same as our cts object
all(counts(dds) == cts.norm)

# Estimate Size factors
dds.norm = estimateSizeFactors(dds)
sizeFactors(dds.norm)
s.facors = sizeFactors(dds.norm)
# Apply size factors to counts and output
normalized.cts = counts(dds.norm, normalized=TRUE)

vsd.norm = vst(dds.norm, blind=F)

# Output
write.csv(assay(vsd.norm), "./data/vsd_counts.csv", row.names=TRUE)
#write.csv(normalized.cts, "./data/normalized_counts.csv", row.names=TRUE)
#write.csv(as.data.frame(s.facors), "./data/sizeFactors.csv")
