# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("khroma") 
library("patchwork")
library("scales")
library("Seurat")
library("SeuratDisk")
library("tidyverse")
library("viridis")


parser <- ArgumentParser(description = "Differential expression analyses")
parser$add_argument('--seurat', '-s',
                    help = 'Filename of Seurat object')
parser$add_argument('--figures',  '-f', 
                    help = 'Directory for saving figures')
args <- parser$parse_args()

# ----------------------------------------------------------------------------
# UMAP by sample ----
umap_by_sample <- function(seurat_obj,
                           fig_dir,
                           reduction=c("umap", "umap.full"),
                           assay=c("sketch", "RNA")){
    assay <- match.arg(assay)
    DefaultAssay(seurat_obj) <- assay
    reduction <- match.arg(reduction)
    out_fname <- file.path(fig_dir, sprintf("%s_%s.pdf", assay, reduction))

    # All samples on single umap
    pdf(out_fname)
    p <- DimPlot(seurat_obj,
                 group.by = "Sample",
                 reduction = reduction)
    print(p)
    dev.off()
    
    # One plot per sample
    pdf(gsub(".pdf", "_by_sample.pdf", out_fname), height = 12, width = 9)
    p <- DimPlot(seurat_obj,
                 group.by = "seurat_clusters",
                 split.by = "Sample",
                 reduction = reduction,
                 ncol = 3)
    print(p)
    dev.off()
}

# UMAP with selected markers ----
# marker set 1 ----
marker_set_1 <- function(){
    return(c("BCL2", "CCR7", "TCF7", "PRF1", "CXCR3",
             "GZMB", "IL2RB", "CD27", "SELL"))
}

# umap_set1 ----
umap_set1 <- function(seurat_obj, fig_dir){
    marker_set_1 <- marker_set_1()
    
    pdf(file.path(fig_dir, "sketch_selected_markers_1.pdf"), width = 8)
    fp <- FeaturePlot(seurat_obj, features = marker_set_1, keep.scale = "all",
                      cols = c("lightgray", "gold", "firebrick")) +
        plot_layout(guides = "collect") &
        theme(axis.line = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank()) &
        labs(color = "log2 Exp")
    print(fp)
    dev.off()
}

# marker set 2 ----
marker_set_2 <- function(){
    return(c("ANXA1", "KIR3DL1", "IL7R", "KLRG1", "TOX", "KLRB1", "CD74", 
             "KLF2", "IFRD1", "MIF", "IRF1", "RUNX3", "ELF1", "AKNA",
             "IFNGR1", "KLRF1", "EOMES", "PDCD1", "SLAMF6", "GZMB", "GZMH",
             "GZMK", "TCF7", "XCL1"))
}

# umap_set2 ----
umap_set2 <- function(seurat_obj, fig_dir){
    marker_set_2 <- marker_set_2()
    
    pdf(file.path(fig_dir, "sketch_selected_markers_2.pdf"),
        width = 12, height = 15)
    
    fp <- FeaturePlot(seurat_obj,
                      reduction = "umap.full", 
                      features = marker_set_2,
                      keep.scale = "all",
                      cols = c("lightgray", "gold", "firebrick"),
                      alpha = 0.5) +
        plot_layout(guides = "collect") &
        theme(axis.line = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              plot.title = element_text(size = 10)) &
        labs(color = "log2 Exp")
    print(fp)
    dev.off()
}

# Dotplots for selected markers ----
dotplot_selected_markers <- function(seurat_obj, dotplot_dir, dotplot_markers){
    # Alternative method - use combine = FALSE and lapply themes 
    
    by_sample <- SplitObject(seurat_obj, split.by = "Sample")
    
    for (sample in names(by_sample)){
        
        fn_template <- file.path(dotplot_dir, "dotplot_%s_genes_on_%.pdf")
        
        pdf(sprintf(fn_template, sample, "y"), height = 5)
        p <- DotPlot(by_sample[[sample]], features = rev(dot_markers)) +
            coord_flip() + 
            scale_color_viridis_c() +
            theme(axis.text.y = element_text(size = 12),
                  axis.text.x = element_text(angle = 90, size = 14,
                                             hjust = 1, vjust = 0.5)) +
            labs(x = NULL, y = NULL)
        print(p)
        dev.off()
        
        pdf(sprintf(fn_template, sample, "x"), height = 5)
        
        p <- DotPlot(by_sample[[sample]], features = rev(dot_markers)) +
            scale_color_viridis_c() +
            theme(axis.text.y = element_text(size = 14),
                  axis.text.x = element_text(angle = 90, size = 12,
                                             hjust = 1, vjust = 0.5)) +
            labs(x = NULL, y = NULL)
        print(p)
        dev.off()
    }
}

# make_plots ----
make_plots <- function(args){
    umap_dir <- file.path(args$figures, "umap")
    dotplot_dir <- file.path(args$figures, "dotplot")
    
    if (! dir.exists(args$figures)) { dir.create(args$figures) }
    if (! dir.exists(umap_dir)) { dir.create(umap_dir) }
    #if (! dir.exists(dotplot_dir)) { dir.create(dotplot_dir) }
    
    seurat_obj <- read_rds(args$seurat)
    
    umap_by_sample(seurat_obj,
                   umap_dir,
                   reduction = "umap",
                   assay = "sketch")
    
    umap_by_sample(seurat_obj,
                   umap_dir,
                   reduction = "umap.full",
                   assay = "RNA")
    
    DefaultAssay(seurat_obj) <- "sketch"
    umap_set1(seurat_obj, umap_dir)
    umap_set2(seurat_obj, umap_dir)
    
    # need markers for dotplot
}

# ----------------------------------------------------------------------------

make_plots(args)
