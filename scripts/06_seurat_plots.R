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


dot_markers <- c("ANXA1", "KIR3DL1", "IL7R", "KLRG1", "TOX", "KLRB1", "CD74", 
                 "KLF2", "IFRD1", "MIF", "IRF1", "RUNX3", "ELF1", "AKNA",
                 "IFNGR1", "KLRF1", "EOMES", "PDCD1", "SLAMF6", "GZMB", "GZMH",
                 "GZMK", "TCF7", "XCL1")

seurat_rpca <- read_rds(file.path(scratch_dir, "integrated_sketch_rpca.rds"))
dotplot_dir <- file.path(fig_dir, "dotplots/")

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

