# Libraries and setup ----

library("ggplot2")
library("khroma") 
library("patchwork")
library("scales")
library("Seurat")
library("tidyverse")
library("viridis")

project_dir <- tryCatch({Sys.getenv()[["project_dir"]]},
                        error = function(cond){return(".")})
fig_dir <- file.path(project_dir, "figures")


# NOTE - THIS IS WITH THE ORIGINAL CLUSTERING 
seurat_rpca <- read_rds(file.path(project_dir,
                                  "data/processed/integrated_sketch_rpca.rds"))

# Get sequences for clone of interest
source(file.path(project_dir, "scripts/clones_of_interest.R"))

# UMAP of cluster 2 ----

# Names of graphs to use in FindSubCluster names(seurat_rpca@graphs)
# FindSubCluster using sketch_snn or sketch_nn doesn't assign clusters to all

cl2_subs <- subset(seurat_rpca,
                   cells = which(seurat_rpca$seurat_clusters == "2"))

DefaultAssay(cl2_subs) <- "RNA"
cl2_subs[["RNA"]] <- JoinLayers(cl2_subs[["RNA"]])

cl2 <- CreateSeuratObject(GetAssayData(cl2_subs, "RNA"),
                          meta.data = cl2_subs[[]])

cl2[["RNA"]] <- split(cl2[["RNA"]], f = cl2$Sample)
cl2 <- NormalizeData(cl2)
cl2 <- FindVariableFeatures(cl2)
cl2 <- ScaleData(cl2)
cl2 <- RunPCA(cl2)

# Integrate with RPCA ----
cl2 <- IntegrateLayers(cl2,
                       method = RPCAIntegration,
                       orig = "pca",
                       new.reduction = "integrated.rpca",
                       dims = 1:30,
                       k.anchor = 20)

cl2 <- FindNeighbors(cl2,
                     reduction = "integrated.rpca",
                     dims = 1:30)
cl2 <- FindClusters(cl2,
                    cluster.name = "subcluster")
cl2 <- RunUMAP(cl2,
               reduction = "integrated.rpca",
               dims = 1:30)

write_rds(cl2, file = file.path(project_dir,
                                "data/processed/cl2_subclusters.rds"))
