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
parser$add_argument('--input', '-i',
                    help = 'Seurat object')
parser$add_argument('--figures',  '-f', 
                    help = 'Directory for saving figures')
args <- parser$parse_args()

# ----------------------------------------------------------------------------
# UMAP by sample ----
umap_by_sample_sketch <- function(seurat_obj, fig_dir,
                                  assay=c("sketch", "RNA")){
    
    assay <- match.arg(assay)
    pdf(file.path(fig_dir, sprintf("%s.pdf", assay)))
    p <- DimPlot(seurat_obj,
                 group.by = "Sample",
                 reduction = "umap")
    print(p)
    dev.off()
    
    pdf(file.path(fig_dir, "sketch_by_sample.pdf"), height = 12, width = 9)
    p <- DimPlot(seurat_obj,
                 group.by = "seurat_clusters",
                 split.by = "Sample",
                 reduction = "umap",
                 ncol = 3)
    print(p)
    dev.off()
}



# UMAP of integrated data, all samples ----
DefaultAssay(seurat_rpca) <- "RNA"

pdf(file.path(fig_dir, "umap/sketch_full.pdf"))
p <- DimPlot(seurat_rpca,
             reduction = "umap.full", 
             group.by = "Sample")
print(p)
dev.off()


pdf(file.path(fig_dir, "umap/sketch_full_cl_by_sample.pdf"),
    height = 12, width = 9)
p <- DimPlot(seurat_rpca,
             group.by = "seurat_clusters",
             split.by = "Sample",
             reduction = "umap.full",
             ncol = 3)
print(p)
dev.off()


# UMAP with selected markers ----
umap_set1 <- function(seurat_obj, fig_dir){
    marker_set_1 <- c("BCL2", "CCR7", "TCF7", "PRF1", "CXCR3",
                      "GZMB", "IL2RB", "CD27", "SELL")

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

umap_set2 <- function(seurat_obj, fig_dir){
    marker_set_2 <- 
        c("ANXA1", "KIR3DL1", "IL7R", "KLRG1", "TOX", "KLRB1", "CD74", 
          "KLF2", "IFRD1", "MIF", "IRF1", "RUNX3", "ELF1", "AKNA",
          "IFNGR1", "KLRF1", "EOMES", "PDCD1", "SLAMF6", "GZMB", "GZMH",
          "GZMK", "TCF7", "XCL1")
    
    
    pdf(file.path(fig_dir, "sketch_selected_markers_2.pdf"),
        width = 12, height = 15)
    
    fp <- FeaturePlot(seurat_rpca,
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
    print(pf)
    dev.off()
}

make_plots <- function(args){
    umap_dir <- file.path(args$figures, "umap")
    dotplot_dir <- file.path(args$figures, "dotplot")
    if (! dir.exists(umap_dir)) { dir.create(umap_dir) }
    if (! dir.exists(dotplot_dir)) { dir.create(dotplot_dir) }
    
    seurat_obj <- read_rds(args$input)
    DefaultAssay(seurat_obj) <- "sketch"
    
    umap_set1(seurat_obj, umap_dir)
    umap_set2(seurat_obj, umap_dir)
    
    umap_by_sample_sketch(seurat_obj, fig_dir)
    
}
# ----------------------------------------------------------------------------

#fig_dir <- file.path(project_dir, "figures")
#dotplot_dir <- file.path(fig_dir, "dotplots/")


# UMAP on the integrated, sketched cells






# Alternative method - use combine = FALSE and lapply themes 

# Dotplots for selected markers ----

by_sample <- SplitObject(seurat_rpca, split.by = "Sample")

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

