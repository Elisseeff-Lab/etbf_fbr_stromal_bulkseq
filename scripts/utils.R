# Utility functions
here::i_am("scripts/utils.R")

#' Convert Genes Using Table
#' 
#' Takes a vector of gene inputs and returns converted gene table
#' Uses Ensembl Genes 111 with human genes (GRCh38) and Homologous
#' mouse genes (GRCm39), downloaded from Ensembl Biomart (2022)
#' 
#' @param genes The genes to convert
#' @param from  Gene symbol type of the input (ENSG, ENSMUSG, HGNC, MGI)
#' @param to    Desired gene symbol type for the output (HGNC, MGI)
#' @return Data frame of genes with original and corresponding converted symbols

table_convert_genes <- function(genes, from, to) {
    conversion_table = read.csv(here::here("refs", "gene_conv.csv"))
    if (from == 'ENSMUSG'){
        col1 = conversion_table$mm.ens 
    }
    if (from == 'ENSG'){
        col1 = conversion_table$hs.ens
    }
    if (from == 'MGI'){
        col1 = conversion_table$mgi   
    }
    if (from == 'HGNC'){
        col1 = conversion_table$hgnc
    }
    if (to == 'MGI'){
        col2 = conversion_table$mgi
    }
    if (to == 'HGNC'){
        col2 = conversion_table$hgnc
    }
    genesV2 = cbind(col1[which(col1 %in% genes)], col2[which(col1 %in% genes)])
    return(genesV2)
}

#' Convert GMT pathway genes
#' 
#' Converts the genes in a GMT file from one gene symbol type to another
#' 
#' @param fgsea_pathway The GMT file to convert
#' @param from_gene The gene symbol type of the input (HGNC, MGI)
#' @param to_gene The desired gene symbol type for the output (HGNC, MGI)
#' @return The converted list, suitable for writing to GMT file
convert_gmt_genes <- function(fgsea_pathway, from_gene, to_gene) {
    for (i in 1:length(fgsea_pathway)){
        # Convert genes
        if (from_gene != to_gene){
            tmp = table_convert_genes(fgsea_pathway[[i]], from = from_gene, to = to_gene)[ ,2]
        } else {
            tmp = fgsea_pathway[[i]]
        }
        # Remove empty genes
        tmp = tmp[tmp != ""]
        # Select distinct genes
        tmp = unique(tmp)
        fgsea_pathway[[i]] = tmp
    }
    return(fgsea_pathway)
}

#' This function will perform fast preranked gene set encrichment analysis
#' 
#' @param gene_sets     Directory to GMT file or a list of gene sets
#' @param rank_data     Data frame with gene names and data to be used in ranking
#' @param rank_by       User defined data ranking by logFC, Pvalue, or logFC_Pvalue (Default setting at
#'                      -log10(padj)*log2FC)
#' @param from_gene     "ENSG" or "ENSMUSG" or "HGNC" or "MGI"
#' @param to_gene       "MGI" or HGNC
#' @return  Outputs dataframe from fgsea function

rank_gse <- function(gene_sets, rank_data, from_gene, to_gene, 
    rank_by = 'logFC_Pvalue', proc = 4, nperm = NULL, set_min = NULL){

    if (typeof(rank_data) == 'character' || typeof(rank_data) == 'list'){
        
        if (typeof(rank_data) == 'character'){
            gse_rank = readRDS(paste0(rank_data))
        } else {
            gse_rank = rank_data
        }

        gse_sub = gse_rank[which(gse_rank$cluster == cluster_name),]

        if (rank_by == 'Pvalue'){
            rank = gse_sub$p_val_adj
        } else if(rank_by == 'sign_Pvalue'){
            rank = -log10(gse_sub$p_val_adj) * gse_sub$avg_log2FC
        } else{
            rank = gse_sub$avg_log2FC
        }

        fgsea_pathway = fgsea::gmtPathways(directory)
        lapply(fgsea_pathway, function(x) {
            # Convert genes
            if (from_gene != to_gene){
            tmp = convert_genes(x, from = from_gene, to = to_gene)
            tmp = tmp[,2]
            } else {
                tmp = x
            }
            # Remove empty gene
            tmp = tmp[tmp != ""]
            # Find all the unique gene
            tmp = unique(tmp)
            # If we want a minimal amount of overlap:
            if (!is.null(set_min)){
                # Check if enough overlap between set and genes in object:
                overlap = length(intersect(tmp, gse_sub$gene))
                if (overlap >= set_min) {
                    x = tmp
                } else {
                    x = NA
                    print(paste0("Only ", overlap, " genes overlapping."))
                }
            } else {
                x = tmp
            }
            return(x)
        })

        # Have to check for any NA values from set_min:
        fgsea_pathway = fgsea_pathway[!is.na(fgsea_pathway)]

        names(rank) <- gse_sub$gene
        
        # Eliminate the Inf and -Inf in the rank
        if (max(rank) == Inf){
            tmp = which(rank == Inf)
            rank = rank[-tmp]
        }
        if (min(rank) == -Inf){
            tmp = which(rank == -Inf)
            rank = rank[-tmp]
        }

        fgseaResult = fgsea::fgsea(pathways = fgsea_pathway, stats = rank, nproc = proc, maxSize = 500)
    }
    else{
        print("Please input the correct inputs: either directory of RDS file or dataframe")
    }

    return(fgseaResult)
}