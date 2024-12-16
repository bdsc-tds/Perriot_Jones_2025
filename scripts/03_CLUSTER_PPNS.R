# Libraries and setup ----

library("ggplot2")
library("khroma") 
library("DESeq2") 
library("patchwork")
library("pheatmap")
library("scales")
library("Seurat")
library("tidyverse")
library("viridis")

set.seed(511697)
ndim <- 15
p_cutoff <- 0.05 
tcr_col <- "beta_aa"

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--input', '-i',
                    help = 'seurat object')
parser$add_argument('--figures',  '-f', 
                    help = 'Directory for saving figures')
parser$add_argument('--output', '-o', # de_dir
                    help = 'directory for writing differential expression') 
args <- parser$parse_args()

# ----------------------------------------------------------------------------
# Functions ----

# Proportion of cells from each cluster per donor ----
cells_per_donor <- function(seurat_obj, fig_dir){
    seurat_obj$seurat_clusters <-
        as.factor(as.numeric(seurat_obj$seurat_clusters))
    Idents(seurat_obj) <- seurat_obj$seurat_clusters
    
    # IS DF METADATA?
    pdf(file.path(fig_dir, "sample_by_cluster.pdf"), width = 9)
    ggplot(df, aes(x = Sample, fill = Cluster)) +
        geom_bar(position = "fill", color = "black") +
        scale_y_continuous(expand = expansion(c(0,0))) +
        scale_x_discrete(expand = expansion(c(0,0))) +
        coord_flip() + 
        scale_fill_manual(values = default_subs[levels(df$Cluster)],
                          labels = names(default_subs[levels(df$Cluster)])) +
        theme_bw() + 
        labs(x = NULL, y = "Proportion of cells")
    dev.off()
}

# run_de ----
run_de <- function(args){
    seurat_obj <- read_rds(args$input)
    
    # Setup named default and rainbow palettes
    cl_names <- levels(Idents(seurat_obj))
    n_clusters <- length(cl_names)
    default_pal <- structure(scales::hue_pal()(n_clusters),
                             names = cl_names)
    rainbow_pal <- structure(rev(color("discrete rainbow")(n_clusters)),
                             names = cl_names)
    cells_per_donor(seurat_obj, fig_dir)
    
}
