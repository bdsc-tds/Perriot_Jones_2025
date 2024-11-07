# Libraries and setup ----

library("ComplexHeatmap")
library("ggplot2")
library("khroma") 
library("patchwork")
library("scales")
library("Seurat")
library("tidyverse")
library("viridis")

project_dir <- tryCatch({Sys.getenv()[["project_dir"]]},
                        error = function(cond){return(".")})
fig_dir <- file.path(project_dir, "figures/subcluster_coi/")


# NOTE - THIS IS WITH THE ORIGINAL CLUSTERING 
seurat_rpca <- read_rds(file.path(project_dir,
                                  "data/processed/integrated_sketch_rpca.rds"))

# Get sequences for clone of interest
source(file.path(project_dir, "scripts/clones_of_interest.R"))

marker_set_2 <- c("ANXA1", "KIR3DL1", "IL7R", "KLRG1", "TOX", "KLRB1", "CD74", 
                  "KLF2", "IFRD1", "MIF", "IRF1", "RUNX3", "ELF1", "AKNA",
                  "IFNGR1", "KLRF1", "EOMES", "PDCD1", "SLAMF6", "GZMB", "GZMH",
                  "GZMK", "TCF7", "XCL1")

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

cl2[["RNA"]] <- JoinLayers(cl2[["RNA"]])

write_rds(cl2, file = file.path(project_dir,
                                "data/processed/cl2_subclusters.rds"))

# Add UMAP coordinates ----
um <- Embeddings(cl2[["umap"]])
identical(rownames(um), rownames(cl2[[]]))
cl2[[]] <- dplyr::bind_cols(cl2[[]], um)

# UMAP by sample and cluster ----

pdf(file.path(fig_dir, "umap_by_sample_and_cluster.pdf"),
    height = 5, width = 10)
p1 <- DimPlot(cl2, group.by = "Sample") 
p2 <- DimPlot(cl2, group.by = "subcluster") 
p1 + p2
dev.off()

# UMAP with clones of interest ----

# Are the clones of interest restricted to one subcluster ----
x <- cl2[[]] %>%
    dplyr::group_by(Sample, beta_aa, subcluster, coi) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::filter(! coi == "no") %>%
    dplyr::arrange(desc(n))


meta <- cl2[[]]

meta$is_coi = c(2, 0.5)[as.numeric(meta$coi == "no") + 1]
meta <- meta %>%
    dplyr::mutate(coi = factor(coi, levels = c("no", "Ri01", "Ri02"))) %>%
    dplyr::arrange(coi) %>%
    dplyr::mutate(coi = as.character(coi))

coi_colours <- structure(c("lightgray","#DC050C","steelblue"),
                         names = c("no", "Ri01", "Ri02"))


pdf(file.path(fig_dir, "umap_coi.pdf"),
    height = 12, width = 12)
ggplot(meta, 
       aes(x = umap_1, y = umap_2)) +
    geom_point(aes(size = is_coi, colour = coi)) +
    scale_size_identity() +
    #facet_wrap(~Sample) +
    theme_minimal(base_size = 15) +
    scale_color_manual(labels = names(coi_colours), values = coi_colours) +
    labs(x = "UMAP dimension 1", y = "UMAP dimension 2") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

dev.off()

# Heatmap expression of selected markers by sample ----

# Subset to just the clones with at least 10 cells
temp <- cl2[[]] %>%
    dplyr::group_by(Sample, beta_aa) %>%
    dplyr::mutate(top_clone = 
                    case_when(n() >= 10 & ! beta_aa == "NA_NA" ~ TRUE,
                              TRUE ~ FALSE))
cl2_top_clones <- subset(cl2, cells = Cells(cl2)[temp$top_clone])

# Aggregate across sample and clone

pseudo_cl2 <- AggregateExpression(cl2_top_clones,
                                  assays = "RNA",
                                  return.seurat = TRUE,
                                  group.by = c("Sample", "beta_aa"))

pseudo_cl2_by_sample <- SplitObject(pseudo_cl2, split.by = "Sample")

plots <- lapply(names(pseudo_cl2_by_sample), function(nm){
    DotPlot(pseudo_cl2_by_sample[[nm]],
            features = marker_set_2) +
        theme(axis.text.y = element_text(size = 8),
              axis.text.x = element_text(size = 4,
                                         angle = 90,
                                         hjust = 1,
                                         vjust = 0.5)) +
        coord_flip() + 
        scale_color_viridis_c() +
        labs(title = nm)
    
})




pdf(file.path(fig_dir, "dotplots_top_clones.pdf"),
    height = 15, width = 15)
wrap_plots(plots) + plot_layout(guides = "collect")
dev.off()


anno <- pseudo_cl2[[]][, c("Sample", "beta_aa")]
mtx <- GetAssayData(pseudo_cl2)[marker_set_2, ]

row_anno = rowAnnotation(Sample = anno$Sample)

# Raw expression
Heatmap(t(mtx),
        row_names_gp = gpar(fontsize = 2),
        right_annotation = row_anno)

# Markers scaled
pdf(file.path(fig_dir, "heatmap_cluster_2_selected_markers.pdf"),
    height = 12, width = 9)
Heatmap(scale(t(mtx)),
        row_names_gp = gpar(fontsize = 4),
        right_annotation = row_anno)

dev.off()

# Are clones specific to cluster 2 or distributed across clusters ----


# Differential expression between clone of interest and others ----

cl2$is_coi = ifelse(cl2$coi == "no", "non-coi", "coi")

Idents(cl2) <- cl2$is_coi
coi_de <- FindAllMarkers(object = cl2) %>%
    dplyr::filter(p_val_adj <= 0.01, cluster == "coi") %>%
    readr::write_csv(file = file.path(fig_dir, "de_coi_v_other_cl2.csv"))

Idents(cl2) <- cl2$subcluster
subcluster_de <- FindAllMarkers(object = cl2) %>%
    dplyr::filter(p_val_adj <= 0.01) %>%
    readr::write_csv(file = file.path(fig_dir, "de_cl2_subclusters.csv"))

# # Pseudobulk - not enough samples? Plus error from DESeq2 for non-integer
# coi_f <- gsub("-", "_", c(ri01_coi, ri02_coi))
# pseudo_cl2$is_coi <- gsub("-", "_", pseudo_cl2$beta_aa) %in% coi_f
# 
# Idents(pseudo_cl2) <- "is_coi"
# coi_de <- FindMarkers(object = pseudo_cl2, 
#                       ident.1 = TRUE, 
#                       ident.2 = FALSE,
#                       test.use = "DESeq2")