# Libraries and setup ----

library("ComplexHeatmap")
library("ggplot2")
library("ggrepel")
library("janitor")
library("khroma") 
library("patchwork")
library("PCAtools")
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
cl2 <- ScaleData(cl2, features = Features(cl2))
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

# Write table of top clones
cl2_top_clones[[]] %>%
    dplyr::group_by(Sample, beta_aa) %>%
    dplyr::summarise(n_cells = n()) %>%
    readr::write_csv(file = file.path(fig_dir, "cluster_2_top_clones.csv"))

# Aggregate across sample and clone
# Note that Seurat appears to use "average" by default -
# AggregateExpression wraps PseudobulkExpression which has default method
# "average"
pseudo_cl2 <- AggregateExpression(cl2_top_clones,
                                  assays = "RNA",
                                  return.seurat = TRUE,
                                  group.by = c("Sample", "beta_aa"))

write_rds(pseudo_cl2, file = file.path(project_dir,
                                "data/processed/cl2_pseudobulk.rds"))



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


# Differential expression between clone of interest and others ----

cl2$is_coi = ifelse(cl2$coi == "no", "non-coi", "coi")

Idents(cl2) <- cl2$Sample
cl2_ri01_02 <- subset(cl2, idents = c("Ri01_dis", "Ri02"))
Idents(cl2_ri01_02) <- cl2_ri01_02$is_coi

coi_v_other_cl2_ri_samples <- FindAllMarkers(object = cl2_ri01_02) %>%
    dplyr::filter(p_val_adj <= 0.05, cluster == "coi") %>%
    dplyr::select(-cluster) %>%
    dplyr::relocate(gene)

readr::write_csv(coi_v_other_cl2_ri_samples,
                 file = file.path(fig_dir, "de_coi_v_other_cl2_ri01_dis_ri02.csv"))

# Differential expression between Ri and Healthy samples ----

cl2$Condition <- gsub("[0-9]+(_dis|_5m)?$", "", cl2$Sample)
Idents(cl2) <- cl2$Condition

ri_v_healthy <- FindAllMarkers(object = cl2) %>%
    dplyr::filter(p_val_adj <= 0.05, cluster == "Ri") %>%
    dplyr::select(-cluster) %>%
    dplyr::relocate(gene) %>%
    write_csv(file = file.path(fig_dir, "de_Ri_v_healthy_cl_2.csv"))

# Heatmap of significant genes, grouped by clone ----

# TCR by cell and sample

# Top clones
ri01_02_top_clones <- cl2_ri01_02[[]] %>%
    dplyr::semi_join(cl2_top_clones[[]])

# Check whether any clones are top in one sample but not in the other
cl2_ri01_02[[]] %>%
    dplyr::filter(beta_aa %in% ri01_02_top_clones$beta_aa) %>%
    dplyr::anti_join(ri01_02_top_clones)
    
Idents(cl2_ri01_02) <- cl2_ri01_02$beta_aa
cl2_ri01_02_top_clones <- subset(cl2_ri01_02,
                                 cells = Cells(cl2_ri01_02)[cl2_ri01_02$beta_aa %in% 
                                     ri01_02_top_clones$beta_aa],
                                 features = coi_v_other_cl2_ri_samples$gene)


Idents(cl2_ri01_02_top_clones) <- cl2_ri01_02_top_clones$beta_aa

pdf(file.path(fig_dir, "heatmap_cl2_ri01_02_top_clones.pdf"),
    width = 12)
DoHeatmap(cl2_ri01_02_top_clones,
               features = coi_v_other_cl2_ri_samples$gene,
          angle = 90, size = 2) +
    guides(color = FALSE) +
    theme(axis.text.y = element_text(size = 6))
dev.off()


# # Replicate heatmap
# scale_data <- GetAssayData(cl2_ri01_02_top_clones, layer = "scale.data")
# colnames(scale_data) <- NULL
# 
# column_ha = HeatmapAnnotation(clone = cl2_ri01_02_top_clones$beta_aa)
# dend1 = cluster_between_groups(scale_data, cl2_ri01_02_top_clones$beta_aa)
# 
# Heatmap(scale_data,
#         cluster_columns = dend1,
#         show_row_dend = FALSE,
#         show_column_dend = FALSE,
#         row_names_gp = gpar(fontsize = 8),
#         top_annotation = column_ha)

# By clone of interest by cell and sample
Idents(cl2_ri01_02) <- cl2_ri01_02$coi
pdf(file.path(fig_dir, "heatmap_coi_v_other_cl2_ri01_dis_ri02.pdf"))
DoHeatmap(cl2_ri01_02,
          features = coi_v_other_cl2_ri_samples$gene) 
dev.off()

# Differential expression coi per sample ----

# Ri01_dis

Idents(cl2) <- "Sample"
cl2_ri01_dis <- subset(cl2, idents = c("Ri01_dis"))
Idents(cl2_ri01_dis) <- cl2_ri01_dis$is_coi

coi_v_other_cl2_ri01_dis <- FindAllMarkers(object = cl2_ri01_dis) %>%
    dplyr::filter(p_val_adj <= 0.05, cluster == "coi") %>%
    dplyr::select(-cluster) %>%
    dplyr::relocate(gene)

readr::write_csv(coi_v_other_cl2_ri01_dis,
                 file = file.path(fig_dir, "de_coi_v_other_cl2_ri01_dis.csv"))


# Ri01_5m

Idents(cl2) <- "Sample"
cl2_ri01_5m <- subset(cl2, idents = c("Ri01_5m"))
Idents(cl2_ri01_5m) <- cl2_ri01_5m$is_coi

coi_v_other_cl2_ri01_5m <- FindAllMarkers(object = cl2_ri01_5m) %>%
    dplyr::filter(p_val_adj <= 0.05, cluster == "coi") %>%
    dplyr::select(-cluster) %>%
    dplyr::relocate(gene)

readr::write_csv(coi_v_other_cl2_ri01_dis,
                 file = file.path(fig_dir, "de_coi_v_other_cl2_ri01_5m.csv"))

# Clone of interest, _dis v _5m 

Idents(cl2) <- "Sample"
cl2_ri01 <- subset(cl2, idents = c("Ri01_5m", "Ri01_dis"))
Idents(cl2_ri01) <- "coi"
cl2_ri01_coi <- subset(cl2_ri01, idents = c("Ri01"))
Idents(cl2_ri01_coi) <- "Sample"

ri01_dis_v_5m <- FindAllMarkers(object = cl2_ri01_coi) %>%
    dplyr::filter(p_val_adj <= 0.1, cluster == "Ri01_dis") %>%
    dplyr::select(-cluster) %>%
    dplyr::relocate(gene)

readr::write_csv(ri01_dis_v_5m,
                 file = file.path(fig_dir, "de_ri01_dis_v_5m.csv"))



# Other DE analyses ----
Idents(cl2) <- cl2$is_coi
coi_de <- FindAllMarkers(object = cl2) %>%
    dplyr::filter(p_val_adj <= 0.01, cluster == "coi") %>%
    readr::write_csv(file = file.path(fig_dir, "de_coi_v_other_cl2.csv"))

Idents(cl2) <- cl2$subcluster
subcluster_de <- FindAllMarkers(object = cl2) %>%
    dplyr::filter(p_val_adj <= 0.01) %>%
    readr::write_csv(file = file.path(fig_dir, "de_cl2_subclusters.csv"))

# Ri01_dis coi versus Ri01_dis other clones ----

# Are the clones in cluster 2 unique to cluster 2?

meta <- read_csv(file.path(project_dir,
                           "results/default_res/integrated_sketch_rpca_md.csv"))

meta_top_clones <- meta %>%
    # Cluster proportions
    dplyr::group_by(seurat_clusters) %>%
    dplyr::mutate(n_cluster = n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(cluster_ppn = n_cluster/n()) %>%
    
    # Keep clones of interest
    dplyr::semi_join(cl2_top_clones[[]] %>%
                          dplyr::select(Sample, beta_aa)) %>%
    
    # Cluster distribution for clones of interest
    dplyr::group_by(Sample, beta_aa) %>%
    dplyr::mutate(n_sample_clone = n()) %>%
    dplyr::group_by(Sample, beta_aa, seurat_clusters,
                    n_cluster, n_sample_clone, cluster_ppn) %>%
    dplyr::summarise(n_sample_clone_cluster = n()) %>%
    dplyr::mutate(pct_sample_clone_cluster = 
                      n_sample_clone_cluster/n_sample_clone * 100,
                  
                  # Adjust for different cluster sizes 
                  adj_ppn = pct_sample_clone_cluster / cluster_ppn) %>%
    dplyr::ungroup()


for (sample in unique(meta_top_clones$Sample)){ 

    pdf(file.path(fig_dir, "clones_by_cluster", 
                  sprintf("top_clones_by_cluster_%s.pdf", sample)),
        height = 10, width = 10)
    
    p <- ggplot(meta_top_clones %>% filter(Sample == sample),
           aes(x = pct_sample_clone_cluster,
               y = beta_aa,
               fill = as.factor(seurat_clusters))) + 
        geom_bar(stat = "identity", color = "black") + 
        theme_bw() +
        facet_wrap(~Sample, scales = "free_y")
    
    print(p)
    
    dev.off()
}



for (sample in unique(meta_top_clones$Sample)){ 
    
    pdf(file.path(fig_dir, "clones_by_cluster_adj", 
                  sprintf("top_clones_by_cluster_adj_%s.pdf", sample)),
        height = 10, width = 10)
    
    p <- ggplot(meta_top_clones %>% filter(Sample == sample),
                aes(x = adj_ppn,
                    y = beta_aa,
                    fill = as.factor(seurat_clusters))) + 
        geom_bar(stat = "identity", 
                 position = "fill",
                 color = "black") + 
        theme_bw() +
        facet_wrap(~Sample, scales = "free_y") +
        labs(y = NULL, y = "Proportion, adjusted for cluster size")
    
    print(p)
    
    dev.off()
}

meta_top_clones %>%
    janitor::clean_names() %>%
    dplyr::select(- n_cluster) %>%
    dplyr:: mutate(across(where(is.numeric),
                          \(x) tibble::num(x, digits = 2))) %>%
    write_csv(file.path(fig_dir, "top_clone_cluster_proportions.csv"))


# PCA with markers of interest ----

pca_markers <- c("KIR3DL1", "ANXA1", "IL7R", "CD74", "TOX")
cl2_pca1 <- subset(cl2_top_clones, features = pca_markers)

cl2_pca1 <- RunPCA(cl2_pca1, features = pca_markers)

pdf(file.path(fig_dir, 
              "subcluster_coi/pca_cl2_top_clones_markers_by_cell.pdf"), height = 10)
FeaturePlot(cl2_pca1, reduction = "pca", features = pca_markers)
dev.off()


cl2_pca3 <- RunPCA(cl2_pca1, features = setdiff(pca_markers, "TOX"))

pdf(file.path(fig_dir, 
              "subcluster_coi/pca_cl2_top_clones_markers_no_tox_by_cell.pdf"))
FeaturePlot(cl2_pca3, reduction = "pca", features = setdiff(pca_markers, "TOX"))
dev.off()

# Group similar clones
cl2_pca3[[]] <- cl2_pca3[[]] %>%
    dplyr::mutate(tg = gsub("(TRBV[0-9]+-)[0-9]", "\\1\\*", TCR2),
                  tg = gsub("(.*TRBJ[0-9]+-)[0-9]", "\\1\\*", tg)) %>%
    tidyr::separate_wider_delim(tg, delim = ".", names_sep = "_" ) %>%
    dplyr::mutate(gp_nm = sprintf("%s.*.%s.%s", tg_1, tg_3, tg_4)) %>%
    dplyr::select(-matches("^tg_[0-9]")) %>%
    dplyr::arrange(gp_nm)

pdf(file.path(fig_dir, 
              "subcluster_coi/pca_cl2_top_clones_by_clone_no_tox_by_cell.pdf"),
    width = 8, height = 10)

# Colour by cluster name
Idents(cl2_pca3) <- cl2_pca3$gp_nm
DimPlot(cl2_pca3, reduction = "pca", pt.size = 1) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 8)) + 
    guides(color = guide_legend(ncol = 3,
                                override.aes = list(size = 3)))
dev.off()


# PCA with expression aggregated by clone ----

Idents(cl2_pca1) <- cl2_pca1$beta_aa
cl2_agg <- AggregateExpression(cl2_pca1, return.seurat = TRUE)
cl2_agg_no_tox <- RunPCA(cl2_agg, features = setdiff(pca_markers, "TOX"))

pdf(file.path(fig_dir, 
              "subcluster_coi/pca_cl2_top_clones_markers_no_tox.pdf"))
FeaturePlot(cl2_agg_no_tox,
            reduction = "pca",
            features = setdiff(pca_markers, "TOX"))
dev.off()

cl2_agg_tox <- RunPCA(cl2_agg, features = pca_markers)

pdf(file.path(fig_dir, 
              "subcluster_coi/pca_cl2_top_clones_markers.pdf"), height = 10)
FeaturePlot(cl2_agg_tox, reduction = "pca", features = pca_markers)
dev.off()

pca_coords <- Embeddings(cl2_agg_tox) %>%
    tibble::as_tibble(rownames = "clone") %>%
    dplyr::mutate(coi = clone %in% gsub("_", "-", c(ri01_coi, ri02_coi)))


highlight_points = pca_coords %>% 
    filter((PC_1 >= -0.25 & PC_1 <= 0.5 & PC_2 >= -0.75 & PC_2 <= 0.0) |
              (PC_1 >= 0 & PC_1 <= 0.5 &  PC_2 >= 1.50 & PC_2 <= 2))
    

pdf(file.path(fig_dir, "subcluster_coi/pca_cl2_top_clones_with_tox.pdf"),
    width = 10, height = 10)
ggplot(pca_coords,
       aes(x = PC_1, y = PC_2, color = coi)) +
    geom_point(size = 2) +
    theme_minimal() +
    theme(axis.line = element_line(size = 0.5,
                                     linetype = "solid",
                                     colour = "black")) +
    labs(x = "PC 1", y = "PC 2") +
    scale_color_manual(values = c("lightgray", "red")) +
    
    geom_text_repel(data = highlight_points,
                    aes(label = clone), size = 2.5,
                    max.overlaps = Inf,
                    min.segment.length = 0,
                    color = "black",
                    nudge_x = .15,
                    nudge_y = .25) +
    # Get rid of an outlier for better viewing
    coord_cartesian(ylim = c(-2.5, 2.5))
dev.off()


# Heatmap all donors, top clones ----

# (note that Seurat DoHeatmap doesn't cluster rows or columns)

agg_dat <- GetAssayData(pseudo_cl2,
                        layer = "scale.data")[coi_v_other_cl2_ri_samples$gene, ]

raw_dat <- GetAssayData(pseudo_cl2,
                        layer = "data")[coi_v_other_cl2_ri_samples$gene, ]

pseudo_cl2[[]] <- pseudo_cl2[[]] %>%
    dplyr::mutate(coi = case_when(beta_aa == gsub("_", "-", ri01_coi) ~ "Ri01",
                  beta_aa == gsub("_", "-", ri02_coi) ~ "Ri02",
                  TRUE ~ "no"))

column_ha = HeatmapAnnotation(Clone = pseudo_cl2$coi,
                              Sample = pseudo_cl2$Sample)

pdf(file.path(fig_dir, "heatmap_cl2_clones_scale_data.pdf"),
    width = 12, height = 10)
Heatmap(agg_dat,
        top_annotation = column_ha,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 4),
        show_row_dend = FALSE)#,
        #col = PurpleAndYellow()) # purple and yellow looks mostly purple

dev.off()


pdf(file.path(fig_dir, "heatmap_cl2_clones_log_count_data.pdf"),
    width = 12, height = 10)
Heatmap(raw_dat,
        top_annotation = column_ha,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 4),
        show_row_dend = FALSE)#,
#col = PurpleAndYellow()) # purple and yellow looks mostly purple

dev.off()
