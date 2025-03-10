# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("tidyverse")
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

get_markers <- function(){
    c("CCR7", "TCF7", "SELL", "BCL2", "CD27",
      "CXCR3", "IL2RB", "PRF1", "GZMB")
}

# UMAP clone of interest ----
umap_coi <- function(meta, fig_dir){
    meta$is_coi = c(2, 0.5)[as.numeric(meta$coi == "no") + 1]
    meta <- meta %>%
        dplyr::mutate(coi = factor(coi, levels = c("no", "Ri01", "Ri02"))) %>%
        dplyr::arrange(coi) %>%
        dplyr::mutate(coi = as.character(coi))
    
    coi_colours <- structure(c("lightgray","#DC050C","steelblue"),
                             names = c("no", "Ri01", "Ri02"))
    
    pdf(file.path(fig_dir, "umap_coi.pdf"),
        height = 12, width = 12)
    p <- ggplot(meta, 
                aes(x = umap_1, y = umap_2)) +
        geom_point(aes(size = is_coi, colour = coi)) +
        scale_size_identity() +
        theme_minimal(base_size = 15) +
        scale_color_manual(labels = names(coi_colours), values = coi_colours) +
        labs(x = "UMAP dimension 1", y = "UMAP dimension 2") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    print(p)
    dev.off()
}


#feature_plot <- function(out_fname, markers = get_markers(),
# width = 8, height =){
#    
#}

## UMAP of cluster markers ----
#ht <- ifelse(marker_nm == "A_bis", 7, 16)
#pdf(file.path(args$results, sprintf("marker_%s_umap.pdf", marker_nm)),
#    width = 8, height = ht)
#p <- FeaturePlot(marker_subset,
#                 features = markers$Gene,
#                 keep.scale = "all",
#                 alpha = 0.5,
#                 cols = c("#f2f2f2", "#ffae00", "#811818"),
#                 order = TRUE) +
#    plot_layout(guides = "collect") &
#    theme(axis.line = element_blank(),
#          axis.title = element_blank(),
#          axis.text = element_blank(),
#          axis.ticks = element_blank()) &
#    labs(color = "log2 Exp")
#print(p)
#dev.off()

# main ----
main <- function(args){
    results <- file.path(args$results, "umap")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    
    seurat_obj <- read_rds(args$seurat)
    seurat_obj$ri01_coi <- seurat_obj$coi == "Ri01"
    
    # UMAP all donors ----
    Idents(seurat_obj) <- "Sample"
    pdf(file.path(results, "samples_all_samples_all_cells.pdf"))
    
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full",
                 group.by = "Sample")
    print(p)
    dev.off()

    # UMAP clusters, except Ri01_5m ----
    seurat_subs <- subset(seurat_obj, Sample != "Ri01_5m")
    Idents(seurat_obj) <- "Sample"
    pdf(file.path(results, "samples_excluding_Ri01_5m.pdf"))
    
    DimPlot(seurat_obj, reduction = "umap.full")
    dev.off()
    
    # UMAP, only Ri01_dis highlighted ----
    seurat_obj$Ri01_5m <- seurat_obj$Sample == "Ri01_5m"
    Idents(seurat_obj) <- "Ri01_5m"
    pdf(file.path(results, "Ri01_5m_highlighted.pdf"))
    
    DimPlot(seurat_obj,
            reduction = "umap.full",
            order = TRUE,
            alpha = 0.5,
            cols = c("#f2f2f2", "#e60000"))
    dev.off()
    
    # UMAP, clone of interest ----
    Idents(seurat_obj) <- "ri01_coi"
    pdf(file.path(results, "ri01_clone_of_interest_all_samples.pdf"))
    
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full",
                 order = TRUE,
                 cols = c("#f2f2f2", "#e60000"))
    print(p)
    dev.off()
    
    # Feature plot, marker set A_bis ----
    markers <- get_markers()
    pdf(file.path(results, "A_bis_markers.pdf"), width = 12, height = 12)
    
    p <- FeaturePlot(seurat_obj,
                     features = markers)
    print(p)
    dev.off()
    
}

# ----------------------------------------------------------------------------
main(args)
