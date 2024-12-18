# Libraries and setup ----

library("argparse")
library("ggplot2")
library("khroma") 
library("patchwork")
library("tidyverse")
library("viridis")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--input', '-i',
                    help = 'Metadata from seurat object')
parser$add_argument('--figures',  '-f', 
                    help = 'Directory for saving figures')

# ----------------------------------------------------------------------------
# Functions ----

# Proportion of cells from each cluster per donor ----
cells_per_donor <- function(md, fig_dir){
    pdf(file.path(fig_dir, "sample_by_cluster.pdf"), width = 9)
    p <- ggplot(md, aes(x = Sample, fill = seurat_clusters)) +
        geom_bar(position = "fill", color = "black") +
        scale_y_continuous(expand = expansion(c(0,0))) +
        scale_x_discrete(expand = expansion(c(0,0))) +
        coord_flip() + 
        scale_fill_manual(values = default_subs[levels(df$Cluster)],
                          labels = names(default_subs[levels(df$Cluster)])) +
        theme_bw() + 
        labs(x = NULL, y = "Proportion of cells")
    print(p)
    dev.off()
}

# Bar cluster per sample vertical ----
ppn_by_sample_vertical <- function(cl_pct, fig_dir){
    pdf(file.path(fig_dir, "bar_sample_pct_per_cl_vertical.pdf"),
        height = 10)
    p <- ggplot(cl_pct, aes(x = factor(seurat_clusters,
                                  levels = sort(unique(seurat_clusters))),
                       y = cl_pct, fill = Sample)) +
        geom_bar(color = "black", stat = "identity") +
        facet_wrap(~Sample, ncol = 1) + 
        scale_y_continuous(expand = expansion(c(0,0))) +
        scale_x_discrete(expand = expansion(c(0,0))) +
        theme_minimal() + 
        theme(panel.grid = element_blank(),
              axis.text.y = element_text(size = 6)) + 
        labs(x = NULL, y = "Percentage of cells") +
        guides(fill = FALSE)
    print(p)
    dev.off()
}


# Cluster proportions by sample ----
ppn_by_sample <- function(cl_pct, fig_dir){
    pdf(file.path(fig_dir, "bar_sample_pct_per_cluster.pdf"))
    p <- ggplot(cl_pct, aes(x = factor(seurat_clusters,
                                  levels = rev(sort(unique(seurat_clusters)))),
                       y = cl_pct, fill = Sample)) +
        geom_bar(color = "black", stat = "identity") +
        facet_wrap(~Sample) + 
        scale_y_continuous(expand = expansion(c(0,0))) +
        scale_x_discrete(expand = expansion(c(0,0))) +
        coord_flip() + 
        theme_minimal() + 
        theme(panel.grid = element_blank(),
              axis.text.y = element_text(size = 6)) + 
        labs(x = NULL, y = "Percentage of cells") +
        guides(fill = FALSE)
    print(p)
    dev.off()
}


# Normalise by total counts per sample ----
get_cluster_percent <- function(md){
    md %>%
        dplyr::select(Sample, seurat_clusters) %>%
        dplyr::group_by(Sample) %>%
        dplyr::mutate(n_sample = n()) %>%
        dplyr::group_by(Sample, seurat_clusters, n_sample) %>%
        dplyr::summarise(n_sample_cl = n()) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(cl_pct = n_sample_cl/n_sample * 100) 
} 

# Bar cluster per sample----
bar_cluster_per_sample <- function(md, fig_dir){
    pdf(file.path(fig_dir, "bar_cluster_per_sample.pdf"))
    p <- ggplot(md, aes(fill = seurat_clusters, x = Sample)) +
        geom_bar(position = "fill", color = "black") +
        scale_y_continuous(expand = expansion(c(0,0))) +
        scale_x_discrete(expand = expansion(c(0,0))) +
        coord_flip() + 
        theme_bw() + 
        labs(x = NULL, y = "Proportion of cells")
    print(p)
    dev.off()
}

# Plot cluster by sample unscaled ----
cluster_by_sample_unscaled <- function(md, fig_dir){
    pdf(file.path(fig_dir, "bar_sample_per_cluster_unscaled.pdf"))
    p <- ggplot(md, aes(x = seurat_clusters, fill = Sample)) +
        geom_bar(color = "black") +
        scale_y_continuous(expand = expansion(c(0,0))) +
        scale_x_discrete(expand = expansion(c(0,0))) +
        coord_flip() + 
        theme_bw() + 
        labs(x = NULL, y = "Number of cells")
    dev.off()
}

# Plot cluster by sample ----
cluster_by_sample <- function(md, fig_dir){
    pdf(file.path(fig_dir, "bar_sample_per_cluster_.pdf"))
    p <- ggplot(meta, aes(x = seurat_clusters, fill = Sample)) +
        geom_bar(position = "fill", color = "black") +
        scale_y_continuous(expand = expansion(c(0,0))) +
        scale_x_discrete(expand = expansion(c(0,0))) +
        coord_flip() + 
        theme_bw() + 
        labs(x = NULL, y = "Proportion of cells")
    print(p)
    dev.off()
}

make_plots <- function(args){
    if (! file.exists(args$figures)) { dir.create(args$figures) }
    
    # Metadata with UMAP coordinates
    md <- read_csv(args$input)
    
    cluster_by_sample(md, fig_dir)
    cluster_by_sample_unscaled(md, fig_dir)
    bar_cluster_per_sample(md, fig_dir)
    
    cl_pct <- get_cluster_percent(md)
    ppn_by_sample_vertical(cl_pct, fig_dir)
    ppn_by_sample(cl_pct, fig_dir)
    cells_per_donor(md, fig_dir)
}

# ----------------------------------------------------------------------------
make_plot(args)
