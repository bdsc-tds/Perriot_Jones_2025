# Libraries and setup ----

library("Seurat")
library("tidyverse")

# Add top neighbours of clones of interest to metadata and save ----
coi_cells <- subset(seurat_rpca,
                    cells = which(! seurat_rpca$coi == "no"))

seurat_rpca <- FindNeighbors(seurat_rpca, return.neighbor = TRUE)

nbr <- seurat_rpca[["sketch.nn"]]

# Check that all cells with clone of interest are in the sketch 
table(Cells(coi_cells) %in% Cells(nbr))

# Note - this gets the neighbours within the sketch

n_nbrs <- 10

coi_nbrs <- lapply(Cells(coi_cells), function(cell){ 
    TopNeighbors(nbr, cell, n = n_nbrs)
})

coi_nbrs <- tibble(cell = rep(Cells(coi_cells), each = n_nbrs),
                   neighbour = unlist(coi_nbrs),
                   nbr_coi = neighbour %in% Cells(coi_cells)) %>%
    dplyr::mutate(beta_aa = seurat_rpca[[]][cell, "beta_aa"],
                  nbr_beta_aa = seurat_rpca[[]][neighbour, "beta_aa"])

write_csv(coi_nbrs,
          file = file.path(project_dir, "results/coi_10_neighbours.csv"))

coi_nbrs %>%
    dplyr::group_by(beta_aa, nbr_beta_aa) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::arrange(beta_aa, desc(n)) %>%
    write_csv(file = file.path(project_dir, "results/coi_10_nbr_counts.csv"))


# Add column to metadata indicating if a cell is a coi neighbour
# (without specifying which cell)

seurat_rpca$coi_nbr <- rownames(seurat_rpca[[]]) %in% coi_nbrs$neighbour 

write_csv(seurat_rpca[[]],
          "integrated_sketch_rpca_res_0_5_md.csv")

# ---------------
meta <- seurat_rpca[[]] 
# ---------------
# Differential expression ----

# Add condition * cluster to metadata
seurat_rpca$Condition <- as.factor(gsub("[0-9]+(-dis)?$", "", seurat_rpca$Sample))
seurat_rpca$Cond_Cl <- paste(seurat_rpca$Condition,
                             seurat_rpca$seurat_clusters,
                             sep = "_")

DefaultAssay(seurat_rpca) <- "RNA"
seurat_rpca[["RNA"]] <- JoinLayers(seurat_rpca[["RNA"]])

# Pseudobulk
pseudo_rpca <- AggregateExpression(seurat_rpca,
                                   assays = "RNA",
                                   return.seurat = TRUE,
                                   group.by = c("Sample", "seurat_clusters"))

# Differential expression between clusters ----
Idents(pseudo_rpca) <- "seurat_clusters"
bulk_cl_markers <- FindAllMarkers(object = pseudo_rpca,
                                  min.pct = 0,
                                  test.use = "DESeq2")

write_csv(bulk_cl_markers,
          file.path(de_results, "cluster_markers_all_samples.csv"))

features <- bulk_cl_markers %>%
    dplyr::filter(! grepl("^TR[AB]", gene),
                  p_val_adj <= 0.01,
                  abs(avg_log2FC) > 0.5) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = 5) %>%
    pull(gene) %>%
    unique()

# Cluster heatmap for pseudobulked data
pdf(file.path(fig_dir, "sketch_rpca_cluster_heatmap_pseudobulk.pdf"),
    width = 12, height = 12)
DoHeatmap(pseudo_rpca,
          features = features) +
    theme(axis.text = element_text(size = 8))
dev.off()

#seurat_rpca <- ScaleData(seurat_rpca)

# Differential expression between clusters, without Ri01-5m ----

Idents(pseudo_rpca) <- pseudo_rpca$Sample
pseudo_subs <- subset(pseudo_rpca,
                      idents = setdiff(unique(pseudo_rpca$Sample), c("Ri01-5m")))


# Differential expression between HD and CL, per cluster ----

pseudo_subs$Condition <- as.factor(gsub("[0-9]+(-dis)?$", "", pseudo_subs$Sample))

Idents(pseudo_subs) <- pseudo_subs$Condition
pseudo_rpca_by_cl <- SplitObject(pseudo_subs,
                                 split.by = "seurat_clusters")


hd_vs_ri <- lapply(seq_along(pseudo_rpca_by_cl), function(i){
    mk <- FindAllMarkers(object = pseudo_rpca_by_cl[[i]],
                         min.pct = 0,
                         test.use = "DESeq2") 
    if (nrow(mk) > 0){
        mk <- mk %>%
            # Filter so logFC always refers to HD vs Ri
            dplyr::filter(cluster == "HD") %>%
            dplyr::mutate(cluster = i)
    }
    mk
})

hd_vs_ri <- bind_rows(hd_vs_ri)

write_csv(hd_vs_ri,
          file.path(de_results, "cluster_markers_hd_v_ri_without_Ri01_5m.csv"))

hd_vs_ri_sig <- hd_vs_ri %>%
    dplyr::filter(p_val_adj <= 0.01,
                  abs(avg_log2FC) > 0.5)

write_csv(hd_vs_ri_sig,
          file.path(de_results,
                    "cluster_markers_hd_v_ri_without_Ri01_5m_sig.csv"))



# Run scRepertoire positional entropy

clusters_w_roi = seurat_rpca[[]] %>%
    dplyr::filter(! coi == "no") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n >= 5) %>%
    dplyr::pull(seurat_clusters) %>%
    unique()

coi_cl_subs <- subset(seurat_rpca,
                      subset = seurat_clusters %in% clusters_w_roi)
coi_cl_subs <- SplitObject(coi_cl_subs, split.by = "seurat_clusters")

dummy <- lapply(names(coi_cl_subs), function(nm){
    pdf(file.path(fig_dir, 
                  sprintf("positional_entropy_cl_%s.pdf", nm)))
    scRepertoire::positionalEntropy(coi_cl_subs)
    dev.off()
})

# Error
#  Please provide rownames to the matrix before converting to a Graph
# WORKING HERE 


# 




# TO DO - adjust TCR for cellbender ----

