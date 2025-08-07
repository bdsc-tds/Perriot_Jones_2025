# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("tidyverse")
library("scales")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', '-s',
                    help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
args <- parser$parse_args()


# ----------------------------------------------------------------------------
# Functions ----

# match_cluster_colour ----
# Find the colour of the cluster of interest in the original UMAP
match_cluster_colour <- function(cl_names, coi){
    # The default colours
    n_clusters <- length(cl_names)
    default_pal <- structure(scales::hue_pal()(n_clusters),
                             names = cl_names)
    return(default_pal[coi])
}

# Get the cluster with the most clones of interest ----
get_cluster_coi <- function(md){
    md %>%
        dplyr::filter(is_coi == "coi") %>%
        dplyr::group_by(seurat_clusters) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::arrange(desc(n)) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::pull(seurat_clusters)
}

# get_markers ----
get_markers <- function(){
    c("CCR7", "TCF7", "SELL", "BCL2", "CD27",
      "CXCR3", "IL2RB", "PRF1", "GZMB")
}


# main ----
main <- function(args){
    results <- file.path(args$results, "umap")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    
    seurat_obj <- read_rds(args$seurat)
    
    # Figure 4b - UMAP, colour by clusters ----
    Idents(seurat_obj) <- "seurat_clusters"
    pdf(file.path(results, "fig_4b.pdf"))
    
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full") +
        labs(x = "UMAP 1", y = "UMAP 2")
    print(p)
    dev.off()
    
    
    # Figure 4c UMAP showing clone and cluster of interest highlighted ----
    
    seurat_obj$ri01_coi <- seurat_obj$coi == "Ri01"
    
    Idents(seurat_obj) <- "seurat_clusters"
    coi_cl <- get_cluster_coi(seurat_obj[[]])
    coi_cl_col <- match_cluster_colour(levels(Idents(seurat_obj)), coi_cl)
    seurat_obj[[]] <- seurat_obj[[]] %>%
        dplyr::mutate(cl_coi = case_when(coi == "Ri01" ~ "coi",
                                         seurat_clusters == coi_cl ~ "cl",
                                         TRUE ~ "Other"),
                      cl_coi = factor(cl_coi, levels = c("Other", "cl", "coi")))
    
    
    Idents(seurat_obj) <- "cl_coi"
    
    cols <- structure(c("#CCCCCC", coi_cl_col, "#e60000"),
                      names = c("Other", "cl", "coi"))
    
    pdf(file.path(results, "fig_4c.pdf"))
    
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full",
                 order = c("coi", "cl"),
                 cols = cols,
                 pt.size = 1.5,
                 alpha = 1) +
        labs(x = "UMAP 1", y = "UMAP 2") +
        guides(color = "none")
    print(p)
    dev.off()
    
    
    # Supplementary figure 1A - UMAP showing markers of interest ----
    markers <- get_markers()
    pdf(file.path(results, "supp_fig_1a.pdf"), width = 12, height = 12)
    
    p <- FeaturePlot(seurat_obj,
                     features = markers)
    print(p)
    dev.off()
    

    # Supplementary figure 2A ----
    # UMAP, coloured and split by samples, without remission sample
    seurat_subs <- subset(seurat_obj, Sample != "Ri01_5m")
    seurat_subs[[]] <- seurat_subs[[]]  %>%
        dplyr::mutate(Sample = ifelse(Sample == "Ri01_dis", "Ri01", Sample))
    
    Idents(seurat_subs) <- "Sample"
    pdf(file.path(results, "supp_fig_2a.pdf"),
        width = 12, height = 12)
    
    p <- DimPlot(seurat_subs,
                 reduction = "umap.full",
                 split.by = "Sample",
                 ncol = 4) +
        labs(x = "UMAP 1", y = "UMAP 2")
    print(p)
    dev.off()
    

}

# ----------------------------------------------------------------------------
main(args)
