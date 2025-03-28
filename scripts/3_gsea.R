# This script runs ranked gene set enrichment analysis using fgsea.
#
# Input:
# - Differential gene expression results (from 2_dge.R)
# - Gene sets - downloaded from MSigDB v 7.5.1 as GMT file.
#   - REACTOME gene sets (M2 - CP - Reactome subset of CP)
#   - Gene set available in refs directory
#
# Output:
# - GSEA results saved as a csv, including top 10 leading edge genes
# A note on the outputs:
#   Due to a combination of parallelization, no recorded seed in the original analysis,
#   and use of a much older version of fgsea, these results diverge from the original analysis. 
#   The results provided in the results directory (3_reactome_gsea.csv) are from the original
#   analysis, and this script is the same as the one used to generate those original results.

here::i_am("scripts/3_gsea.R")
set.seed(42)
# Load required libraries and scripts
library(fgsea)
source(here::here("scripts", "utils.R"))

# Read in DE results and filter
de_res = read.csv(here::here("results", "2_ETBF_pos_v_neg_dge.csv"))
na_ind = which(is.na(de_res$log2FoldChange) | is.na(de_res$padj))
if (!identical(na_ind, integer(0))) de_res = de_res[-na_ind, ]

# Load gene sets
reactome = gmtPathways(here::here("refs", "m2.cp.reactome.v7.5.1.mgi.gmt"))

# Rank
rank = -log10(de_res$padj) * de_res$log2FoldChange
names(rank) = de_res$X

# Run fgsea
fgsea_res = fgsea::fgsea(pathways = reactome, stats = rank, maxSize = 500)

# Format results to include leading edge genes as well as fgsea scores
table = fgsea_res[order(fgsea_res$NES, decreasing = T), ]
leading = lapply(table$leadingEdge, function(x){x[1:10]})
if (length(leading) != 0){
    lead = matrix(unlist(leading), ncol = 10, byrow = T)
    colnames(lead) = sapply(1:10, function(x)paste0("leading_", x))
    table$leadingEdge = NULL
} else {
    lead = NULL
    table$leadingEdge = NULL
}
table = cbind(table, lead)

# Save results - note that "3_gsea_reactome.csv" is the name of the original results file
write.csv(table, here::here("results", "3_gsea_reactome_new.csv"))
