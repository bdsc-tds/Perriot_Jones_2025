# Libraries and setup ----

library("ggplot2")
library("khroma") 
library("patchwork")
library("tidyverse")
library("viridis")

project_dir <- tryCatch({Sys.getenv()[["project_dir"]]},
                        error = function(cond){return(".")})
fig_dir <- file.path(project_dir, "figures")

# Metadata with UMAP coordinates
meta <- read_csv(file.path(project_dir,
                           "results/integrated_sketch_rpca_md.csv"))


# Plot cluster by sample ----

pdf(file.path(fig_dir, "bar_sample_per_cluster_.pdf"))
ggplot(meta, aes(x = seurat_clusters, fill = Sample)) +
    geom_bar(position = "fill", color = "black") +
    scale_y_continuous(expand = expansion(c(0,0))) +
    scale_x_discrete(expand = expansion(c(0,0))) +
    coord_flip() + 
    theme_bw() + 
    labs(x = NULL, y = "Proportion of cells")
dev.off()


pdf(file.path(fig_dir, "bar_sample_per_cluster_unscaled.pdf"))
ggplot(meta, aes(x = seurat_clusters, fill = Sample)) +
    geom_bar(color = "black") +
    scale_y_continuous(expand = expansion(c(0,0))) +
    scale_x_discrete(expand = expansion(c(0,0))) +
    coord_flip() + 
    theme_bw() + 
    labs(x = NULL, y = "Number of cells")
dev.off()


pdf(file.path(fig_dir, "bar_cluster_per_sample.pdf"))
ggplot(meta, aes(fill = seurat_clusters, x = Sample)) +
    geom_bar(position = "fill", color = "black") +
    scale_y_continuous(expand = expansion(c(0,0))) +
    scale_x_discrete(expand = expansion(c(0,0))) +
    coord_flip() + 
    theme_bw() + 
    labs(x = NULL, y = "Proportion of cells")
dev.off()


# Normalise by total counts per sample
# 
# cl_pct <- seurat_rpca[[]] %>%
#     dplyr::select(Sample, seurat_clusters) %>%
#     dplyr::group_by(Sample) %>%
#     dplyr::mutate(n_sample = n()) %>%
#     dplyr::group_by(Sample, seurat_clusters, n_sample) %>%
#     dplyr::summarise(n_sample_cl = n()) %>%
#     dplyr::ungroup() %>%
#     dplyr::mutate(cl_pct = n_sample_cl/n_sample * 100) 
# 
# # Cluster proportions by sample
# 
# pdf(file.path(fig_dir, "bar_sample_ppn_per_cluster.pdf"))
# ggplot(cl_pct, aes(x = seurat_clusters, y = cl_pct, fill = Sample)) +
#     geom_bar(color = "black", stat = "identity") +
#     scale_y_continuous(expand = expansion(c(0,0))) +
#     scale_x_discrete(expand = expansion(c(0,0))) +
#     coord_flip() + 
#     theme_bw() + 
#     labs(x = NULL, y = "Proportion of cells")
# dev.off()

