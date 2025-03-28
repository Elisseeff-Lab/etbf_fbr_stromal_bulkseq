# This script processes raw counts data:
# - Merging samples into a single matrix
# - Converting gene IDs (ENSMUSG to MGI)
# - Removing duplicate genes
# - Filtering genes with low counts
# - Creating normalized counts per million
#
# Input:
# - Raw counts data from data directory (data/*.ReadsPerGene.out.tab)
#
# Output:
# - Processed counts matrix saved as a csv for downstream analysis
# - Normalized counts per million saved as a csv

here::i_am("scripts/1_processing.R")

# Load required libraries and functions
library(DESeq2)
source(here::here("scripts", "utils.R"))

# Read in counts
files = list.files("data", ".tab$")
# Get sample names from file names
names(files) = vapply(files, function(x) {
    subset = stringr::str_sub(x, start = c(5, 13), end = c(7, 13))
    paste(subset, collapse = "")
    }, character(1))
# Read in data
reads = list()
reads = lapply(names(files), function(x) {
    magrittr::set_names(data.frame(read.table(here::here("data", files[x]), skip = 4)[, 4]), x)
})

# Make matrix
counts = purrr::list_cbind(reads)

# Add row names (genes) - same across all files
genes = read.table(here::here("data", files[1]), skip = 4)[, 1]
trim_gene = sapply(genes, function(genes) {
    split_result = strsplit(genes, spl = ".", fixed = TRUE)[[1]]
    if (length(split_result) > 0) {
        split_result[1]
    } else {
        genes
    }
}, USE.NAMES = FALSE)
rownames(counts) = trim_gene

# Convert gene names from ENSMUSG to MGI
converted = table_convert_genes(rownames(counts), from = "ENSMUSG", to = "MGI")
converted = converted[match(rownames(counts), converted[, 1]), ]
rm_id = which(is.na(converted[, 2]) | converted[, 2] == "")
if (!identical(rm_id, integer(0))) {
    counts = counts[-rm_id, ]
    converted = converted[-rm_id, ]
}
counts = as.matrix(counts)
rownames(counts) = converted[, 2]

# Find duplicates
dups = unique(converted[which(duplicated(converted[, 2])), 2])
# Combine duplicate rows
for (dup in dups) {
    rows = which(rownames(counts) == dup)
    duprows = counts[rows, ]
    counts = counts[-rows, ]
    combined = colSums(duprows)
    counts = rbind(combined, counts)
    rownames(counts)[1] = dup
    print(dup)
}

# Filter genes with no reads:
counts = counts[rowSums(counts) > 0, ]

# Write out the counts matrix
write.csv(counts, here::here("results", "1_raw_counts.csv"))

# Make metadata - no table needed since this is two groups only
metas = c(rep("ETBF_NEG", 4), rep("ETBF_POS", 4))

# Calculate normalized CPM:
coldata = data.frame("sample" = metas)
rownames(coldata) = colnames(counts)
dds = DESeqDataSetFromMatrix(count = counts, col = coldata, design = ~ sample)
dds = estimateSizeFactors(dds)
norm = counts(dds, normalized = TRUE)
write.csv(norm, here::here("results", "1_norm_cpm.csv"))
