
library("argparse")
library("ggplot2")
library("patchwork")
library("tidyverse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', '-s',
                    help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
parser$add_argument('--min_cells',  '-m', default = 5, 
                    help = 'Minimum cells per clone for pseudobulk')
args <- parser$parse_args()

# ----------------------------------------------------------------------------
# Functions ----
filter_tcr_genes <- function(markers){
    markers %>%
        dplyr::filter(! grepl("^TR[ABGD][VDJ]", gene))
} 

# Filter top genes
filter_top <- function(de, p_val = 0.001, n_per_cluster = 10, min_pct = 0.25){
    de %>%
        dplyr::rowwise() %>%
        dplyr::mutate(max_pct = max(`pct.1`, `pct.2`)) %>%
        dplyr::ungroup() %>%
        dplyr::filter(p_val_adj <= p_val,
                      max_pct >= min_pct) %>%
        dplyr::group_by(cluster) %>%
        dplyr::arrange(desc(abs(avg_log2FC))) %>%
        dplyr::slice_head(n = n_per_cluster) %>%
        ungroup()
}

main <- function(args){
    results <- file.path(args$results, "clusters")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    
    seurat_obj <- read_rds(args$seurat)
    
    # Cluster DE, per cell
    Idents(seurat_obj) <- seurat_obj$seurat_clusters
    de <- FindAllMarkers(seurat_obj)
    write_csv(de, file.path(results, "de_full.csv"))
    
    #pdf(file.path(results, "pca_genes.pdf"), width = 10, height = 10)
    #plots <- lapply(1:16, function(i){
    #    DimHeatmap(seurat_obj, dims = i, cells = 500, balanced = TRUE)
    #})
    #wrap_plots(plots)
    #dev.off()
    
    cluster_de <- filter_top(de)
    write_csv(cluster_de, file.path(results, "de_genes_in_heatmap.csv"))
    
    seurat_obj <- ScaleData(seurat_obj, features = cluster_de$gene)
    pdf(file.path(results, "heatmap.pdf"), height = 10, width = 10)
    p <- DoHeatmap(seurat_obj, features = cluster_de$gene,
                   size = 2, angle = 0) +
        theme(axis.text = element_text(size = 4))
    print(p)
    dev.off()
    
    # Cluster DE, pseudobulk
    
    n_beta_aa <- seurat_obj[[]] %>%
        dplyr::group_by(Sample, beta_aa, seurat_clusters) %>%         
        dplyr::summarise(n_clone_sample_cluster = n()) %>% 
        dplyr::filter(n_clone_sample_cluster >= args$min_cells) 

    # Check if same clone is in same cluster

    seurat_subs <- subset(seurat_obj, beta_aa %in% n_beta_aa$beta_aa)
    clones_pseudo <- AggregateExpression(seurat_subs,
                                         return.seurat = TRUE,
                                         group.by = c("Sample",
                                                      "seurat_clusters",
                                                      "beta_aa"))
    write_rds(clones_pseudo, gsub(".rds", "_pseudobulk.rds", args$seurat))
    
    md <-  seurat_subs[[]] %>%
        dplyr::select("Sample", "beta_aa", "seurat_clusters")

    clones_pseudo[[]] <- dplyr::left_join(clones_pseudo[[]], md,
                                          by = c("Sample", "beta_aa"))
    
    Idents(clones_pseudo) <- clones_pseudo$seurat_clusters
    de_pseudo <- FindAllMarkers(seurat_obj, test.use = "DESeq2")
    write_csv(de, file.path(results, "de_full_pseudobulk_sample_clone.csv"))
    
    
    pseudo_de_genes <- filter_top(de_pseudo)
    pdf(file.path(results), "heatmap_pseudobulk.pdf", height = 10)
    p <- DoHeatmap(clones_pseudo, features = pseudo_de_genes$gene)
    dev.off()
}
# ----------------------------------------------------------------------------
main(args)
