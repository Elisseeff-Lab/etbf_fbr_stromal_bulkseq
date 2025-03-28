# This script performs differential gene expression analysis using DESeq2
#
# Input:
# - Processed counts matrix (output from 1_processing.R)
#
# Ouput:
# - Differential gene expression results saved as a csv

here::i_am("scripts/2_dge.R")

# Load required libraries
library(DESeq2)

# Read in counts matrix
counts = read.csv(here::here("results", "1_raw_counts.csv"), row.names = 1)

# Make metadata - no table needed since this is two groups only
metas = c(rep("ETBF_NEG", 4), rep("ETBF_POS", 4))

# Do comparisons using DESeq2
dds = DESeqDataSetFromMatrix(count = counts, col = coldata, design = ~ sample)
dds = DESeq(dds)
res = results(dds, name = "sample_ETBF_POS_vs_ETBF_NEG")
# Using shrinkage to address variance
res_shrink = lfcShrink(dds, coef = "sample_ETBF_POS_vs_ETBF_NEG", type = "apeglm")
res_shrink_sort = res_shrink[order(res_shrink$padj), ]
write.csv(as.data.frame(res_shrink_sort), here::here("results", "2_ETBF_pos_v_neg_dge.csv"))
