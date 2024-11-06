# Libraries and setup ----

library("ggplot2")
library("khroma") 
library("patchwork")
library("scales")
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

cl2_subs <- subset(seurat_rpca, cells = which(seurat_rpca$seurat_clusters == "2"))
cl2_subs <- FindVariableFeatures(cl2_subs)
cl2_subs <- ScaleData(cl2_subs, features = rownames(cl2_subs))
cl2_subs <- RunPCA(cl2_subs)

DefaultAssay(cl2_subs) <- "RNA"

cl2_subs <- IntegrateLayers(
    object = cl2_subs,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca.c2",
    verbose = FALSE
)

cl2_subs <- FindNeighbors(cl2_subs,
                          reduction = "integrated.rpca.c2",
                          dims = 1:30)
cl2_subs <- FindClusters(cl2_subs,
                         #resolution = 2,
                         cluster.name = "subcluster")




cl2 <- FindSubCluster(
    seurat_rpca,
    cluster = "2",
    graph.name = "sketch_nn",
    subcluster.name = "sub.cluster",
    resolution = 0.5,
    algorithm = 1
)

cl2_subs <- subset(cl2, cells = which(seurat_rpca$seurat_clusters == "2"))

DefaultAssay(cl2) <- "RNA"

#cl2 <- FindVariableFeatures(cl2)
#cl2 <- ScaleData(cl2)
#cl2 <- RunPCA(cl2)
cl2 <- RunUMAP(cl2, reduction = "integrated.rpca.full",dims = 1:30)
cl2 <- FindNeighbors(cl2, dims = 1:30,
                     reduction = "integrated.rpca.full",)
cl2 <- FindClusters(cl2)

