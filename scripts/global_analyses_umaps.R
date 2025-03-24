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

make_umaps <- function(seurat_obj,
                       results,
                       condition,
                       gray = "#CCCCCC",
                       pt.size = 1.5){
    # UMAP, colour and split by samples
    Idents(seurat_obj) <- "Sample"
    pdf(file.path(results, sprintf("samples_split_%s.pdf", condition)),
        width = 12)
    
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full",
                 split.by = "Sample") +
        labs(x = "UMAP 1", y = "UMAP 2")
    print(p)
    dev.off()
    
    # UMAP, colour by samples
    Idents(seurat_obj) <- "Sample"
    pdf(file.path(results, sprintf("samples_%s.pdf", condition)))
    
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full") +
        labs(x = "UMAP 1", y = "UMAP 2")
    print(p)
    dev.off()
    
    # UMAP, colour by clusters
    Idents(seurat_obj) <- "seurat_clusters"
    pdf(file.path(results, sprintf("clusters_%s.pdf", condition)))
    
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full") +
        labs(x = "UMAP 1", y = "UMAP 2")
    print(p)
    dev.off()
    
    # UMAP, colour by clusters, alternative palette
    cluster_alphabet <- DiscretePalette(length(unique(
        seurat_obj$seurat_clusters)))
    Idents(seurat_obj) <- "seurat_clusters"
    
    pdf(file.path(results, sprintf("clusters_%s_alphabet_pal.pdf", condition)))
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full",
                 cols = cluster_alphabet) +
        labs(x = "UMAP 1", y = "UMAP 2")
    print(p)
    dev.off()
    
    # UMAP, clone of interest 
    Idents(seurat_obj) <- "ri01_coi"
    pdf(file.path(results, sprintf("ri01_clone_of_interest_%s.pdf", condition)))
    
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full",
                 order = TRUE,
                 cols = c(gray, "#e60000"),
                 pt.size = pt.size) +
        labs(x = "UMAP 1", y = "UMAP 2") +
        guides(color = "none")
    print(p)
    dev.off()
}

# main ----
main <- function(args){
    results <- file.path(args$results, "umap")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    
    seurat_obj <- read_rds(args$seurat)
    seurat_obj$ri01_coi <- seurat_obj$coi == "Ri01"

    make_umaps(seurat_obj, results, "all_samples")
    
    # UMAPs excluding Ri01_5m 
    seurat_subs <- subset(seurat_obj, Sample != "Ri01_5m")
    make_umaps(seurat_subs, results, "excluding_Ri01_5m")
    
    # Just Ri01_dis
    seurat_dis <- subset(seurat_obj, Sample == "Ri01_dis")
    make_umaps(seurat_dis, results, "only_Ri01_dis")
    
    # UMAP, only Ri01_dis highlighted ----
    seurat_obj$Ri01_dis <- seurat_obj$Sample == "Ri01_dis"
    Idents(seurat_obj) <- "Ri01_dis"
    pdf(file.path(results, "Ri01_dis_highlighted.pdf"))
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full",
                 order = TRUE,
                 alpha = 0.5,
                 cols = c("#CCCCCC", "#e60000"))  +
        labs(x = "UMAP 1", y = "UMAP 2")
    print(p)
    dev.off()
    
    # UMAP, only Ri01_dis highlighted ----
    seurat_obj$Ri01_5m <- seurat_obj$Sample == "Ri01_5m"
    Idents(seurat_obj) <- "Ri01_5m"
    pdf(file.path(results, "Ri01_5m_highlighted.pdf"))
    p <- DimPlot(seurat_obj,
                 reduction = "umap.full",
                 order = TRUE,
                 alpha = 0.5,
                 cols = c("#CCCCCC", "#e60000")) +
        labs(x = "UMAP 1", y = "UMAP 2")
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
