# Libraries and setup ----

library("ggplot2")
library("khroma") 
library("patchwork")
library("scales")
library("Seurat")
library("SeuratDisk")
library("tidyverse")
library("viridis")

project_dir <- tryCatch({Sys.getenv()[["project_dir"]]},
                        error = function(cond){return(".")})
scratch_dir <- tryCatch({Sys.getenv()[["scratch_dir"]]},
                        error = function(cond){return(".")})
fig_dir <- file.path(project_dir, "figures")
dotplot_dir <- file.path(fig_dir, "dotplots/")

marker_set_1 <- c("BCL2", "CCR7", "TCF7", "PRF1", "CXCR3",
                  "GZMB", "IL2RB", "CD27", "SELL")

marker_set_2 <- c("ANXA1", "KIR3DL1", "IL7R", "KLRG1", "TOX", "KLRB1", "CD74", 
                 "KLF2", "IFRD1", "MIF", "IRF1", "RUNX3", "ELF1", "AKNA",
                 "IFNGR1", "KLRF1", "EOMES", "PDCD1", "SLAMF6", "GZMB", "GZMH",
                 "GZMK", "TCF7", "XCL1")

seurat_rpca <- read_rds(file.path(scratch_dir, "integrated_sketch_rpca.rds"))

# UMAP on sketched data ----
# UMAP on the integrated, sketched cells
pdf(file.path(fig_dir, "umap_sketch_rpca.pdf"))
DimPlot(seurat_rpca,
        group.by = "Sample",
        reduction = "umap")
dev.off()

pdf(file.path(fig_dir, "umap_sketch_rpca_cluster_by_sample.pdf"),
    height = 12, width = 9)
DimPlot(seurat_rpca,
        group.by = "seurat_clusters",
        split.by = "Sample",
        reduction = "umap",
        ncol = 3)
dev.off()



# UMAP of integrated data, all samples ----
pdf(file.path(fig_dir, "umap_sketch_rpca_full.pdf"))
DimPlot(seurat_rpca,
        reduction = "umap.full", 
        group.by = "Sample")
dev.off()


pdf(file.path(fig_dir, "umap_sketch_rpca_full_cluster_by_sample.pdf"),
    height = 12, width = 9)
DimPlot(seurat_rpca,
        group.by = "seurat_clusters",
        split.by = "Sample",
        reduction = "umap.full",
        ncol = 3)
dev.off()


# UMAP with selected markers ----

pdf(file.path(fig_dir, "umap_sketch_rpca_selected_markers_1.pdf"), width = 8)
fp <- FeaturePlot(seurat_rpca, features = marker_set_1, keep.scale = "all",
                  cols = c("lightgray", "gold", "firebrick")) +
    plot_layout(guides = "collect") &
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) &
    labs(color = "log2 Exp")
print(fp)
dev.off()


pdf(file.path(fig_dir, "umap_sketch_rpca_selected_markers_2.pdf"),
    width = 12, height = 15)

FeaturePlot(seurat_rpca,
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
dev.off()

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

