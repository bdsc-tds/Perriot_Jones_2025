# Libraries and setup ----

library("ggplot2")
library("khroma") 
library("DESeq2") 
library("patchwork")
library("pheatmap")
library("Seurat")
library("tidyverse")
library("viridis")

set.seed(511697)
ndim <- 15
p_cutoff <- 0.05 
tcr_col <- "beta_aa"

# Read integrated data
seurat_obj <- read_rds("data/processed/integrated_sketch_rpca.rds")

# Setup named default and rainbow palettes
cl_names <- levels(Idents(seurat_obj))
n_clusters <- length(cl_names)
default_pal <- structure(scales::hue_pal()(n_clusters),
                         names = cl_names)
rainbow_pal <- structure(rev(color("discrete rainbow")(n_clusters)),
                         names = cl_names)

# Add TCR information ----

combined_tcr <- read_rds(file.path(data_dir, "processed/combined_tcr.rds"))
seurat_obj <- combineExpression(combined_tcr, seurat_obj)

# Add beta chain cdr3
seurat_obj@meta.data <- seurat_obj@meta.data %>%
    tidyr::separate(CTgene,
                    into = c("TCR1", "TCR2"),
                    sep = "_", remove = FALSE) %>%
    tidyr::separate(CTaa,
                    into = c("CTaa1", "CTaa2"),
                    sep = "_", remove = FALSE) %>%
    dplyr::mutate(across(matches("TCR|CTaa"), ~na_if(.x, "NA")),
                  beta_aa = paste(TCR2, CTaa2, sep = "_")) 

# Proportion of cells from each cluster per donor ----

seurat_obj$seurat_clusters <- as.factor(as.numeric(seurat_obj$seurat_clusters))
Idents(seurat_obj) <- seurat_obj$seurat_clusters


pdf("figures/sample_by_cluster.pdf", width = 9)
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


# 2 clones for Ri01 and Ri02

# DEG for each cluster for Ri vs HD

# pseudobulk in seurat
bulk <- AggregateExpression(object,
                            return.seurat = TRUE,
                            slot = "counts",
                            assays = "RNA",
                            group.by = c("celltype.full", "sample", "disease"))
                                                                                                     "sample", "disease"))
