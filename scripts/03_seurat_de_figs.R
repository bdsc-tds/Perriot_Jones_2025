# -------------------------------------------------------------------
# Read data, add clonotypes, filter
# -------------------------------------------------------------------
# Libraries and setup ----

library("scRepertoire")
library("Seurat")
library("SeuratDisk")
library("patchwork")
library("tidyverse")
library("khroma")
library("ggplot2")
library("org.Hs.eg.db")
library("pheatmap")
library("viridis")

set.seed(511697)
ndim <- 15
p_cutoff <- 0.05 
tcr_col <- "beta_aa"

data_dir <- "data"

# Read data ----

combined_tcr <- read_rds(file.path(data_dir, "processed/combined_tcr.rds"))

exp_dirs <- list.files(file.path(data_dir, "raw/GEX"), full.names = TRUE)
exp_dirs <- file.path(exp_dirs, "filtered_feature_bc_matrix")
names(exp_dirs) <- basename(dirname(exp_dirs))

seurat_data <-  Read10X(exp_dirs)

# Update the names to official symbols where necessary ----

features <- rownames(seurat_data)
alias2sym <- AnnotationDbi::select(org.Hs.eg.db, keys = features,
                                   keytype = "ALIAS", columns = "SYMBOL")

alias2sym <- alias2sym %>%
    dplyr::filter(! is.na(SYMBOL) & ! ALIAS == SYMBOL) %>%
    # Don't rename ambiguous features
    dplyr::group_by(ALIAS) %>%
    dplyr::filter(n_distinct(SYMBOL) == 1) %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::filter(n_distinct(ALIAS) == 1) %>%
    # Remove symbols already present in rowname
    dplyr::filter(! SYMBOL %in% features)

features[match(alias2sym$ALIAS, features)] <- alias2sym[["SYMBOL"]]
rownames(seurat_data) <- features

# Create combined Seurat object ----

samples <- data.frame(Sample = gsub("_[ACTG]+-1", "", colnames(seurat_data)),
                      row.names = colnames(seurat_data))
seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                 min.cells = 10,
                                 min.features = 100,
                                 meta.data = samples)

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj,
                                                   pattern = "^MT-")

rm(seurat_data, features, alias2sym, exp_dirs)

# Add TCR information ----

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


# QC plots for all three Ri samples ----

pdf("figures/qc/seurat_violin.pdf")
VlnPlot(seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)
dev.off()

pdf("figures/qc/seurat_feature_scatter.pdf", width = 9)
plot1 <- FeatureScatter(seurat_obj,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_obj,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# Filter high gene count or high mitochondrial read percentage ----
seurat_obj <- subset(seurat_obj,
                     subset = nFeature_RNA < 5000 &
                         percent.mt < 10)

# Keep cells that have any CD8 count or a TCR ----
cd8_expr <- subset(seurat_obj, features = c("CD8A", "CD8B"))
cd8_expressed <- colSums(cd8_expr@assays$RNA$counts) > 0
has_tcr <- ! is.na(seurat_obj@meta.data$CTgene)
keep <- cd8_expressed | has_tcr
seurat_obj <- subset(seurat_obj, cells = which(keep))

# -------------------------------------------------------------------
# Clustering of Ri01_dis and Ri01_5m
# -------------------------------------------------------------------
# Patient 1 at disease (Ri01_dis) and remission (Ri01_5m) ----

Idents(seurat_obj) <- seurat_obj@meta.data$Sample
ri01 <- subset(seurat_obj, idents = c("Ri01_dis", "Ri01_5m"))

# Reset Idents to clusters for Ri01 
Idents(ri01) <- ri01$orig.ident

# Filter, normalize, scale, PCA ----

ri01 <- NormalizeData(ri01)
ri01 <- FindVariableFeatures(ri01)
ri01 <- ScaleData(ri01) 
ri01 <- RunPCA(ri01)

# PCA QC figures ----

pdf("figures/qc/seurat_ri01_pca_elbow.pdf")
ElbowPlot(ri01)
dev.off()

pdf("figures/qc/seurat_ri01_dimheatmap.pdf",
    width = 10, height = 12)
DimHeatmap(ri01, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# Run clustering ----

ri01 <- FindNeighbors(ri01, dims = 1:ndim)
ri01 <- FindClusters(ri01, resolution = 0.5) 
# min.dist default 0.3, increasing preserves more global structure
ri01 <- RunUMAP(ri01, dims = 1:ndim, min.dist = 0.5, n.neighbors = 20)

# Change the cluster IDs to start from 1 

ri01$seurat_clusters <- as.factor(as.numeric(ri01$seurat_clusters))
Idents(ri01) <- ri01$seurat_clusters

# Write cluster IDs ----
write_rds(Idents(ri01),
          file = sprintf("results/clusters/seurat_ri01_ids_%s_dim_res_0_5.rds",
                         ndim))

# -------------------------------------------------------------------
# Analyses for Figure 4 
# -------------------------------------------------------------------
# Plot UMAP for both Ri01 samples ----
pdf("figures/4a_supp_Ri01_dis_Ri01_5m_colour_default.pdf", width = 8)
p <- DimPlot(ri01, reduction = "umap") +
    guides(colour = guide_legend(title = "Cluster",
                                 override.aes = list(size = 6))) +
    labs(x = "UMAP dimension 1", y = "UMAP dimension 2")
p
dev.off()

pdf("figures/4a_supp_Ri01_dis_Ri01_5m_colour_rainbow.pdf", width = 8)
p + khroma::scale_color_discreterainbow(reverse = TRUE)
dev.off()

# Top 10 clones by sample for Ri01 samples ----

ri01_top_10_clones <- ri01@meta.data %>%
    dplyr::filter(! is.na(TCR2)) %>%
    dplyr::group_by(Sample, TCR2, CTaa2, beta_aa) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::group_by(Sample) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::ungroup() %>%
    # Make sure Sample is in the right order for using
    dplyr::mutate(Sample = 
                      factor(Sample, levels = c("Ri01_dis", "Ri01_5m"))) %>%
    dplyr::arrange(Sample) %>%
    dplyr::mutate(plot_label = ifelse(Sample == "Ri01_dis",
                                      sprintf("P1 C%s", 1:n()), NA)) %>%
    dplyr::group_by(beta_aa) %>%
    tidyr::fill(plot_label, .direction = "updown") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(temp = sprintf("P1 C%s", cumsum(is.na(plot_label)) + 10),
                  plot_label = coalesce(plot_label, temp)) %>%
    dplyr::select(-temp)

write_csv(ri01_top_10_clones,
          "results/ri01_label_to_sequence.csv")

ri01@meta.data <- ri01@meta.data %>%
    dplyr::left_join(ri01_top_10_clones %>% 
                         dplyr::select(beta_aa, plot_label) %>%
                         unique(),
                     relationship = "many-to-one")
rownames(ri01@meta.data) <- Cells(ri01)
ri01$plot_label <- factor(ri01$plot_label,
                          levels = unique(ri01_top_10_clones$plot_label))    

# Plot UMAP, highlighting cells from each sample ----

ri01$Sample <- factor(ri01$Sample, levels = c("Ri01_dis", "Ri01_5m"))

p1_dis <- DimPlot(ri01, group.by = "Sample",
                  cols = c("#DC050C", "lightgray")) +
    guides(color = "none") +
    labs(title = "P1 (disease)",
         x = "UMAP dimension 1",
         y = "UMAP dimension 2")

p1_rem <- DimPlot(ri01, group.by = "Sample",
                  cols = c("lightgray", "#DC050C"), order = "Ri01_5m") +
    guides(color = "none") +
labs(title = "P1 (remission)",
     x = "UMAP dimension 1",
     y = "UMAP dimension 2")

pdf("figures/4a_supp_umap_left_Ri01_dis_right_Ri01_5m.pdf",
    width = 9.5, height = 5)
p1_dis + p1_rem
dev.off()


# Split combined object by sample ----

Ri01_tmp <- SplitObject(ri01, split.by = "Sample")
Ri01_dis <- Ri01_tmp[["Ri01_dis"]]
Ri01_rem <- Ri01_tmp[["Ri01_5m"]]

rm(Ri01_tmp)

# Fig 4a - plot cluster UMAPs separately for Ri01_dis and Ri01_5m ----

# Ri01_dis 

pdf("figures/4a_Ri01_dis_colour_default.pdf",
    width = 8)
p <- DimPlot(Ri01_dis, reduction = "umap") +
    guides(colour = guide_legend(title = "Cluster",
                                 override.aes = list(size = 6))) +
    labs(x = "UMAP dimension 1", y = "UMAP dimension 2")
p
dev.off()

pdf("figures/4a_Ri01_dis_colour_rainbow.pdf",
    width = 8)
p + khroma::scale_color_discreterainbow(reverse = TRUE)
dev.off()


# Ri01_5m 

pdf("figures/4a_supp_Ri01_5m_colour_default.pdf",
    width = 8)
p <- DimPlot(Ri01_rem, reduction = "umap") +
    guides(colour = guide_legend(title = "Cluster",
                                 override.aes = list(size = 6))) +
    labs(x = "UMAP dimension 1", y = "UMAP dimension 2")
p
dev.off()

pdf("figures/4a_supp_Ri01_5m_colour_rainbow.pdf",
    width = 8)
p + khroma::scale_color_discreterainbow(reverse = TRUE)
dev.off()


# Fig 4b - selected markers on UMAP ----

fig2b_markers <- c("BCL2", "CCR7", "TCF7", "PRF1", "CXCR3",
                   "GZMB", "IL2RB", "CD27", "SELL")


pdf("figures/fig4b.pdf", width = 8)
FeaturePlot(Ri01_dis, features = fig2b_markers, keep.scale = "all",
            cols = c("lightgray", "gold", "firebrick")) +
    plot_layout(guides = "collect") &
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) &
    labs(color = "log2 Exp")
dev.off()

# Fig 4c ----

fig4c_markers <- c("BCL2", "CCR7", "SELL", "CD27", "TCF7", "PRF1", "CXCR3",
                   "GZMB", "IL2RB", "FAS", "KLRG1", "IL7R")

# average gene expression by beta_aa combination 
av_expr <- AggregateExpression(Ri01_dis, features = fig4c_markers,
                               return.seurat = TRUE)
av_expr <- FetchData(av_expr, vars = fig4c_markers)
av_expr <- t(av_expr)

# Heatmap of the average normalized expression of the markers,
# rows (genes) scaled 

# Heatmap, viridis palette ----
pdf("figures/fig4c_scaled_by_gene_viridis.pdf", width = 8.6)
p <- pheatmap(av_expr,
         labels_col = gsub("g", "", colnames(av_expr)),
         cluster_cols = FALSE,
         angle_col = 0,
         fontsize = 12,
         scale = "row",
         color = viridis::viridis(100))
print(p)
dev.off()

# Heatmap, Seurat palette ----
pdf("figures/fig4c_scaled_by_gene_purple_yellow.pdf", width = 8.6)
p <- pheatmap(av_expr,
              labels_col = gsub("g", "", colnames(av_expr)),
              cluster_cols = FALSE,
              angle_col = 0,
              fontsize = 12,
              scale = "row",
              color = Seurat::PurpleAndYellow())
print(p)
dev.off()

# Fig 4d ----

ri01_top10_nms <- sprintf("P1 C%s", 1:10) 

ri01_dis_top_10 <- ri01_top_10_clones %>%
    dplyr::filter(Sample == "Ri01_dis") %>%
    dplyr::pull(plot_label)
is_top10 <- Ri01_dis$plot_label %in% ri01_dis_top_10

ri01_dis_clones <- subset(Ri01_dis, 
                          cells = Cells(Ri01_dis)[is_top10])
Idents(ri01_dis_clones) <- ri01_dis_clones$Sample
       
Ri01_dis_by_cl <- SplitObject(ri01_dis_clones, split.by = "plot_label")

full_umap <- DimPlot(Ri01_dis,
                     cols = "lightgray",
                     group.by = "Sample",
                     pt.size = 0.25) +
    labs(title = NULL) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none") 

# Adjust the x and y limits to match the full version
# as the reduced version gets zoomed in
xlims <- layer_scales(full_umap)$x$get_limits()
ylims <- layer_scales(full_umap)$y$get_limits()

pdf("figures/fig4d.pdf", width = 8)
dim_plots <- lapply(names(Ri01_dis_by_cl), function(x){
    colour <- ifelse(x == "P1 C8", "#DC050C", "black")
    DimPlot(Ri01_dis_by_cl[[x]], cols = colour, pt.size = 0.25) +
        theme(#axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none") +
        labs(title = x) +
        coord_cartesian(xlim = xlims, ylim = ylims)
})

wrap_plots(c(dim_plots, list(full_umap)))
dev.off()

# The default colours
cl_names <- levels(Idents(Ri01_dis))
n_clusters <- length(cl_names)
default_pal <- structure(scales::hue_pal()(n_clusters),
                         names = cl_names)
rainbow_pal <- structure(rev(color("discrete rainbow")(n_clusters)),
                         names = cl_names)

# Default colour version -----

pdf("figures/fig4d_colour_by_cluster.pdf", width = 8)
dim_plots <- lapply(names(Ri01_dis_by_cl), function(x){
    DimPlot(Ri01_dis_by_cl[[x]], group.by = "seurat_clusters", pt.size = 0.25) +
        theme(#axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none") +
        labs(title = x) +
        scale_color_manual(values = default_pal, labels = names(default_pal)) +
        coord_cartesian(xlim = xlims, ylim = ylims)
})

full_umap <- DimPlot(Ri01_dis, group.by = "seurat_clusters", pt.size = 0.25) +
    labs(title = NULL) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none") 
wrap_plots(c(dim_plots, list(full_umap)))
dev.off()

# Rainbow colour version -----

pdf("figures/fig4d_colour_by_cluster_rainbow.pdf", width = 8)
dim_plots <- lapply(names(Ri01_dis_by_cl), function(x){
    DimPlot(Ri01_dis_by_cl[[x]], group.by = "seurat_clusters", pt.size = 0.25) +
        theme(#axis.line = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none") +
        labs(title = x) +
        scale_color_manual(values = rainbow_pal, labels = names(rainbow_pal)) +
        coord_cartesian(xlim = xlims, ylim = ylims)
})

full_umap <- full_umap +
    scale_color_discreterainbow(reverse = TRUE)
    
wrap_plots(c(dim_plots, list(full_umap)))
dev.off()

rm(Ri01_dis_by_cl)

# Fig 4e ----

df <- tibble(clone = ri01_dis_clones$plot_label,
             Cluster = droplevels(ri01_dis_clones$seurat_clusters)) %>%
    dplyr::mutate(clone = 
                      factor(clone,
                            levels = rev(levels(ri01_dis_clones$plot_label))))

default_subs <- default_pal[as.character(unique(df$Cluster))]
rainbow_subs <- rainbow_pal[as.character(unique(df$Cluster))]


pdf("figures/fig4e_default.pdf", width = 9)
ggplot(df, aes(x = clone, fill = Cluster)) +
    geom_bar(position = "fill", color = "black") +
    scale_y_continuous(expand = expansion(c(0,0))) +
    scale_x_discrete(expand = expansion(c(0,0))) +
    coord_flip() + 
    scale_fill_manual(values = default_subs[levels(df$Cluster)],
                      labels = names(default_subs[levels(df$Cluster)])) +
    theme_bw() + 
    labs(x = NULL, y = "Proportion of cells")
dev.off()


pdf("figures/fig4e_rainbow.pdf", width = 9)
ggplot(df, aes(x = clone, fill = Cluster)) +
    geom_bar(position = "fill", color = "black") +
    scale_y_continuous(expand = expansion(c(0,0))) +
    scale_x_discrete(expand = expansion(c(0,0))) +
    coord_flip() + 
    scale_fill_manual(values = rainbow_subs[levels(df$Cluster)],
                      labels = names(rainbow_subs[levels(df$Cluster)])) +
    theme_bw() + 
    labs(x = NULL, y = "Proportion of cells")
dev.off()

# -------------------------------------------------------------------
# Analyses for Figure 5 
# -------------------------------------------------------------------
# Fig 5a ----
# Diff expr between top 10 clones for Ri01_dis and Ri01_rem combined----

ri01_clones <- subset(ri01,
                      cells = Cells(ri01)[ri01$plot_label %in% ri01_top10_nms])

# Set Idents to clonotypes 
Idents(ri01_clones) <- ri01_clones$plot_label

# Find marker genes ----
ri01_markers <- FindAllMarkers(ri01_clones,
                               min.pct = -Inf,
                               logfc.threshold = -Inf,
                               min.cells.feature = 1,
                               min.cells.group = 1) %>%
    dplyr::arrange(abs(avg_log2FC))

# Write tables of marker genes ----         
ri01_sig <- ri01_markers %>%
    dplyr::filter(#! grepl("^TR[AB]V", gene), # Remove the variable TCR chains
                  # Remove the novel transcripts, lincRNAs
                  ! grepl("^A[CLP][0-9]{4,}|^LOC[0-9]{4,}|^LINC", gene), 
        abs(avg_log2FC) > 1,
        p_val_adj <= 10e-10,
        # Note - filtering by pct here because all values are needed above
        pct.1 >= 0.25 | pct.2 >= 0.25) 

write_csv(ri01_sig, "results/ri01_dis_ri01_5m/markers_top_10_clones_fig5a.csv")

ri01_sig_c8 <- ri01_sig %>%
    dplyr::filter(cluster == "P1 C8") %>%
    dplyr::arrange(p_val_adj)

write_csv(ri01_sig_c8, "results/ri01_dis_ri01_5m/markers_top_10_clones_c8_fig5a.csv")

# Format for heatmap ----
keep_genes <-  ri01_sig %>% 
    dplyr::pull(gene) %>%
    unique()

ri01_clones_agg <- AggregateExpression(ri01_clones,
                                       features = keep_genes,
                                       return.seurat = TRUE)

ri01_clones_expr <- FetchData(ri01_clones_agg, vars = keep_genes)
ri01_clones_expr <- t(ri01_clones_expr)

# Heatmap of markers distinguishing Ri01 clones, viridis palette ----
pdf("figures/fig5a_av_expr_scaled_by_row_viridis.pdf",
    height = 9)
p <- pheatmap(ri01_clones_expr,
         fontsize_row = 1,
         scale = "row",
         color = viridis::viridis(100))
print(p)
dev.off()

# Heatmap of markers distinguishing Ri01 clones, Seurat palette ----
pdf("figures/fig5a_av_expr_scaled_by_row_purple_yellow.pdf",
    height = 9)
p <- pheatmap(ri01_clones_expr,
              fontsize_row = 1,
              scale = "row",
              color = Seurat::PurpleAndYellow())
print(p)
dev.off()

# Most significant genes in P1 C8 ----
# 
# by_pval <- ri01_markers %>%
#     dplyr::filter(cluster == "P1 C8",
#                   abs(avg_log2FC) > 1,
#                   p_val_adj <= 10e-10,
#                   ! grepl("^A[CLP][0-9]{4,}|^LOC[0-9]{4,}|^LINC", gene)) %>%
#     dplyr::arrange(p_val_adj) %>%
#     head(30) 
# 
# by_logfc <- ri01_markers %>%
#     dplyr::filter(cluster == "P1 C8", 
#                   p_val_adj <= 10e-10,
#                   ! grepl("^A[CLP][0-9]{4,}|^LOC[0-9]{4,}|^LINC", gene)) %>%
#     arrange(desc(abs(avg_log2FC))) %>%
#     head(30)
# 
# by_pval_stringent <- ri01_markers %>%
#     dplyr::filter(cluster == "P1 C8",
#                   abs(avg_log2FC) > 2,
#                   p_val_adj <= 10e-10,
#                   ! grepl("^A[CLP][0-9]{4,}|^LOC[0-9]{4,}|^LINC", gene)) %>%
#     head(30) 
# 
# by_pval_pct <- ri01_markers %>%
#     dplyr::filter(cluster == "P1 C8",
#                   abs(avg_log2FC) > 1,
#                   p_val_adj <= 10e-5,
#                   `pct.1` >= 0.1 | `pct.2` >= 0.1,  
#                   ! grepl("^A[CLP][0-9]{4,}|^LOC[0-9]{4,}|^LINC", gene)) 
# 
# by_pval_dot <- DotPlot(ri01_clones, features = by_pval$gene) +
#     coord_flip() + 
#     scale_color_viridis_c() +
#     theme(axis.text.y = element_text(size = 6),
#           axis.text.x = element_text(angle = 90, size = 10,
#                                      hjust = 1, vjust = 0.5))
# print(by_pval_dot)
# 
# 
# by_pval_str_dot <- DotPlot(ri01_clones, features = by_pval_stringent$gene) +
#     coord_flip() + 
#     scale_color_viridis_c() +
#     theme(axis.text.y = element_text(size = 6),
#           axis.text.x = element_text(angle = 90, size = 10,
#                                      hjust = 1, vjust = 0.5))
# print(by_pval_str_dot)
# 

# Fig 5b and 5d ----

dot_markers <- c("ANXA1", "KIR3DL1", "IL7R", "KLRG1", "TOX", "KLRB1", "CD74", 
  "KLF2", "IFRD1", "MIF", "IRF1", "RUNX3", "ELF1", "AKNA", "IFNGR1", 
  "KLRF1", "EOMES", "PDCD1", "SLAMF6", "GZMB", "GZMH", "GZMK", 
  "TCF7", "XCL1")


ri01_clones_split <- SplitObject(ri01_clones, split.by = "Sample")

# Dotplots for Ri01_dis ----

labs <- levels(Idents(ri01_clones_split[["Ri01_dis"]]))
is_bold_dis <- structure(ifelse(labs == "P1 C8", "bold", "plain"),
                         names = labs)

pdf("figures/fig5b_Ri01_dis_genes_on_y.pdf", height = 5)
p <- DotPlot(ri01_clones_split[["Ri01_dis"]], features = rev(dot_markers)) +
    coord_flip() + 
    scale_color_viridis_c() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 90, size = 14,
                                     hjust = 1, vjust = 0.5,
                                     face = is_bold_dis)) +
    labs(x = NULL, y = NULL)
print(p)
dev.off()

pdf("figures/fig5b_Ri01_dis_genes_on_x.pdf", height = 5)

p <- DotPlot(ri01_clones_split[["Ri01_dis"]], features = rev(dot_markers)) +
    scale_color_viridis_c() +
    theme(axis.text.y = element_text(size = 14, face = is_bold_dis),
          axis.text.x = element_text(angle = 90, size = 12,
                                     hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL)
print(p)
dev.off()

# Dotplots for Ri01_5m ----

labs <- levels(Idents(ri01_clones_split[["Ri01_5m"]]))
is_bold <- structure(ifelse(labs == "P1 C8", "bold", "plain"),
                     names = labs)

pdf("figures/fig5d_Ri01_5m_genes_on_y.pdf", height = 5)
p <- DotPlot(ri01_clones_split[["Ri01_5m"]], features = rev(dot_markers)) +
    coord_flip() + 
    scale_color_viridis_c() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 90, size = 12,
                                     hjust = 1, vjust = 0.5,
                                     face = is_bold)) +
    labs(x = NULL, y = NULL)
print(p)
dev.off()

pdf("figures/fig5d_Ri01_5m_genes_on_x.pdf", height = 5)

p <- DotPlot(ri01_clones_split[["Ri01_5m"]], features = rev(dot_markers)) +
    scale_color_viridis_c() +
    theme(axis.text.y = element_text(size = 14, face = is_bold),
          axis.text.x = element_text(angle = 90, size = 12,
                                     hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL)
print(p)
dev.off()


# Fig 5c ----

fig_5c_markers <- c("KIR3DL1", "KLRG1", "TOX", "ANXA1", "IL7R")

Idents(ri01_dis_clones) <- ri01_dis_clones$plot_label
x <- GetAssayData(ri01_dis_clones)[fig_5c_markers, ]
x <- tibble::as_tibble(t(x)) %>%
    dplyr::bind_cols(cluster = Idents(ri01_dis_clones))%>%
    tidyr::pivot_longer(-cluster) %>%
    dplyr::mutate(name = factor(name, levels = fig_5c_markers))

# Violins for selected markers, scaled by area ----
pdf("figures/fig5c_scale_by_area.pdf", height = 4.5, width = 8)    
p <- ggplot(x, aes(x = cluster, y = value, fill = cluster)) +
    geom_violin(width = 1.5, scale = "area") + # This is ggplot default
    theme_bw() +
    facet_wrap(~name, scales = "free_y", ncol = 2) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10), 
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 10,
                                     face = is_bold_dis),
          legend.text = element_text(size = 10),
          legend.position = c(0.75, 0.0),
          legend.key.size = unit(0.75, "line"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 12),
          axis.line = element_line()) +
    guides(fill = guide_legend(ncol = 5)) +
    labs(x = NULL, y = NULL, fill = NULL)
print(p)
dev.off()

# Violins for selected markers, scaled by width ----
pdf("figures/fig5c_scale_by_width.pdf", height = 4.5, width = 8)    
p <- ggplot(x, aes(x = cluster, y = value, fill = cluster)) +
    geom_rect(xmin = 7.5,
              xmax = 8.5,
              ymin = -Inf,
              ymax = Inf,
              fill = "lightgray",
              alpha = 0.25) +
    geom_violin(scale = "width") +
    theme_bw() +
    facet_wrap(~name, scales = "free_y", ncol = 2) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10), 
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 10,
                                     face = is_bold_dis),
          legend.text = element_text(size = 10),
          legend.position = c(0.75, 0.0),
          legend.key.size = unit(0.75, "line"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 12),
          axis.line = element_line()) +
    guides(fill = guide_legend(ncol = 5)) +
    labs(x = NULL, y = NULL, fill = NULL)
print(p)
dev.off()


# -------------------------------------------------------------------
# Analyses for Figure 6
# -------------------------------------------------------------------

# Setup Seurat object for Ri02 ----
Idents(seurat_obj) <- seurat_obj$Sample
ri02 <- subset(seurat_obj, idents = c("Ri02"))

# Reset Idents to clusters for Ri02
Idents(ri02) <- ri02$orig.ident

# Filter, normalize, scale, PCA, UMAP ----

ri02 <- NormalizeData(ri02)
ri02 <- FindVariableFeatures(ri02)
ri02 <- ScaleData(ri02) 
ri02 <- RunPCA(ri02)

ri02 <- FindNeighbors(ri02, dims = 1:ndim)
ri02 <- FindClusters(ri02, resolution = 0.5) 
# min.dist default 0.3, increasing preserves more global structure
ri02 <- RunUMAP(ri02, dims = 1:ndim, min.dist = 0.5, n.neighbors = 20)

# Change the cluster IDs to start from 1 

ri02$seurat_clusters <- as.factor(as.numeric(ri02$seurat_clusters))
Idents(ri02) <- ri02$seurat_clusters

# DimHeatmap for Ri02 ----

DoHeatmap(ri02, features = VariableFeatures(ri02))


# Fig 6a UMAP ----

pdf("figures/fig6a_Ri02_colour_default.pdf", width = 8)
p <- DimPlot(ri02, reduction = "umap") +
    guides(colour = guide_legend(title = "Cluster",
                                 override.aes = list(size = 6))) +
    labs(x = "UMAP dimension 1", y = "UMAP dimension 2")
p
dev.off()

pdf("figures/fig6a_Ri02_colour_rainbow.pdf", width = 8)
p + khroma::scale_color_discreterainbow(reverse = TRUE)
dev.off()

# Fig 6b UMAP with selected markers ----

fig6b_markers <- fig2b_markers 

pdf("figures/fig6b.pdf", width = 8)
FeaturePlot(ri02, features = fig6b_markers, keep.scale = "all",
            cols = c("lightgray", "gold", "firebrick")) +
    plot_layout(guides = "collect") &
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) &
    labs(color = "log2 Exp")
dev.off()



# Get top 10 clones ----

ri02_top_10_clones <- ri02@meta.data %>%
    dplyr::filter(! is.na(TCR2)) %>%
    dplyr::group_by(Sample, TCR2, CTaa2, beta_aa) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::arrange(desc(n)) %>%
    dplyr::ungroup() %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::mutate(plot_label = sprintf("P2 C%s", 1:n()),
                  # This is the same clone as the clone of interest for Ri01
                  is_bold = ifelse(plot_label == "P2 C5", "bold", "plain"), 
                  plot_label = factor(plot_label,
                                      levels = unique(plot_label))) 

write_csv(ri02_top_10_clones,
          "results/ri02_label_to_sequence.csv")

is_bold_ri02 <- structure(ri02_top_10_clones$is_bold,
                          names = as.character(ri02_top_10_clones$plot_label)) 

ri02@meta.data <- ri02@meta.data %>%
    dplyr::left_join(ri02_top_10_clones %>% 
                         dplyr::select(Sample, beta_aa, plot_label))
rownames(ri02@meta.data) <- Cells(ri02)

# Fig 6d ----

ri02_clones <- subset(ri02, cells = which(!is.na(ri02$plot_label)))
Idents(ri02_clones) <- ri02_clones$plot_label

pdf("figures/fig6d_Ri02_genes_on_y.pdf", height = 5)
p <- DotPlot(ri02_clones, features = rev(dot_markers)) +
    coord_flip() + 
    scale_color_viridis_c() +
    theme(axis.text.y = element_text(size = 12),
          axis.text.x = element_text(angle = 90, size = 14,
                                     hjust = 1, vjust = 0.5,
                                     face = is_bold_ri02)) +
    labs(x = NULL, y = NULL)
print(p)
dev.off()

pdf("figures/fig6d_Ri02_genes_on_x.pdf", height = 5)

p <- DotPlot(ri02_clones, features = rev(dot_markers)) +
    scale_color_viridis_c() +
    theme(axis.text.y = element_text(size = 14, face = is_bold_ri02),
          axis.text.x = element_text(angle = 90, size = 12,
                                     hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL)
print(p)
dev.off()

# Fig 6c ----

x <- GetAssayData(ri02_clones)[fig_5c_markers, ]
x <- tibble::as_tibble(t(x)) %>%
    dplyr::bind_cols(cluster = Idents(ri02_clones))%>%
    tidyr::pivot_longer(-cluster)# %>%
    #dplyr::mutate(name = factor(name, levels = fig_3c_markers))


# Violins for selected markers, scaled by area ----
pdf("figures/fig6c_scale_by_area.pdf", height = 4.5, width = 8)    
p <- ggplot(x, aes(x = cluster, y = value, fill = cluster)) +
    geom_violin(scale = "area") + # This is ggplot default
    theme_bw() +
    facet_wrap(~name, scales = "free_y", ncol = 2) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10), 
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 10,
                                     face = is_bold_ri02),
          legend.text = element_text(size = 10),
          legend.position = c(0.75, 0.0),
          legend.key.size = unit(0.75, "line"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 12),
          axis.line = element_line()) +
    guides(fill = guide_legend(ncol = 5)) +
    labs(x = NULL, y = NULL, fill = NULL)
print(p)
dev.off()

# Violins for selected markers, scaled by width ----
pdf("figures/fig6c_scale_by_width.pdf", height = 4.5, width = 8)    
p <- ggplot(x, aes(x = cluster, y = value, fill = cluster)) +
    geom_rect(xmin = 4.5,
              xmax = 5.5,
              ymin = -Inf,
              ymax = Inf,
              fill = "lightgray",
              alpha = 0.25) +
    geom_violin(scale = "width") +
    theme_bw() +
    facet_wrap(~name, scales = "free_y", ncol = 2) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10), 
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     size = 10,
                                     face = is_bold_ri02),
          legend.text = element_text(size = 10),
          legend.position = c(0.75, 0.0),
          legend.key.size = unit(0.75, "line"),
          panel.grid = element_blank(),
          strip.text = element_text(face = "bold", size = 12),
          axis.line = element_line()) +
    guides(fill = guide_legend(ncol = 5)) +
    labs(x = NULL, y = NULL, fill = NULL)
print(p)
dev.off()

# ---------------------------------------------------------------
# Barplot of Ri02 ----

df <- tibble(clone = ri02_clones$plot_label,
             Cluster = droplevels(ri02_clones$seurat_clusters)) %>%
    dplyr::mutate(clone = 
                      factor(clone,
                             levels = rev(levels(ri02_clones$plot_label))))


df %>% group_by(Cluster) %>%
    mutate(n_cluster = n()) %>%
    group_by(clone, Cluster, n_cluster) %>%
    summarise(n = n()) %>%
    group_by(Cluster) %>%
    mutate(pct = n/n_cluster *100) %>%
    write_csv("ri02_clusters.csv")


default_subs <- default_pal[as.character(unique(df$Cluster))]
rainbow_subs <- rainbow_pal[as.character(unique(df$Cluster))]


pdf("figures/Ri02_clone_barplot.pdf", width = 9)
ggplot(df, aes(x = clone, fill = Cluster)) +
    geom_bar(position = "fill", color = "black") +
    scale_y_continuous(expand = expansion(c(0,0))) +
    scale_x_discrete(expand = expansion(c(0,0))) +
    coord_flip() + 
    scale_fill_manual(values = default_subs[levels(df$Cluster)],
                      labels = names(default_subs[levels(df$Cluster)])) +
    theme_bw() + 
    labs(x = NULL, y = "Proportion of cells")
dev.off()


pdf("figures/Ri02_clone_barplot_rainbow.pdf", width = 9)
ggplot(df, aes(x = clone, fill = Cluster)) +
    geom_bar(position = "fill", color = "black") +
    scale_y_continuous(expand = expansion(c(0,0))) +
    scale_x_discrete(expand = expansion(c(0,0))) +
    coord_flip() + 
    scale_fill_manual(values = rainbow_subs[levels(df$Cluster)],
                      labels = names(rainbow_subs[levels(df$Cluster)])) +
    theme_bw() + 
    labs(x = NULL, y = "Proportion of cells")
dev.off()





# Heatmaps of marker genes distinguishing clusters ---------

seurat_by_smp <- list("Ri01" = ri01,
                      "Ri02" = ri02,
                      "Ri01_dis" = Ri01_dis,
                      "Ri01_5m" = Ri01_rem)

for (smp in names(seurat_by_smp)){
    print(smp)
    smp_markers <- FindAllMarkers(seurat_by_smp[[smp]]) %>% 
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1) %>%
        slice_head(n = 10) %>%
        ungroup() 
    
    png(sprintf("figures/seurat_cluster_marker_heatmap_%s.png",
                smp),
        width = 12, height = 10, units = "in", res = 500)
    
    p <- DoHeatmap(seurat_by_smp[[smp]], features = smp_markers$gene) 
        NoLegend() +
        theme(axis.text.y = element_text(size = 6))
    print(p)
    dev.off()
}

# ---------------------------------------------------------------
