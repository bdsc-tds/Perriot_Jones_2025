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

# Depends on 08_clones_of_interest

parser <- ArgumentParser(description = "Tables for clone of interest")
parser$add_argument('--metadata', '-m',
                    help = 'Metadata from seurat object')
parser$add_argument('--figures',  '-o', 
                    help = 'Directory for saving figres')
parser$add_argument('--clones',  '-c', 
                    help = 'R script containing sequences of interest')

args <- parser$parse_args()



# Cluster 2 pseudobulked top clones 
pseudo_cl2 <- read_rds(file = file.path(project_dir,
                                        "data/processed/cl2_pseudobulk.rds")) 

pseudo_cl2[[]] <- pseudo_cl2[[]] %>%
    dplyr::mutate(coi = case_when(beta_aa == gsub("_", "-", ri01_coi) ~ "Ri01",
                                  beta_aa == gsub("_", "-", ri02_coi) ~ "Ri02",
                                  TRUE ~ "no"))

# Markers distinguishing clone of interest from other Ri
coi_v_other_cl2_ri <- readr::read_csv(file = 
        file.path(fig_dir, "de_coi_v_other_cl2_ri01_dis_ri02.csv"))

# Markers distinguishing Ri from healthy 

ri_v_healthy <- read_csv(file = file.path(fig_dir, "de_Ri_v_healthy_cl_2.csv"))

# Manually curated marker sets 

markers <- c("KIR3DL1", "ANXA1", "IL7R", "CD74", "TOX")
rm_markers <- c("TRAV13-1",
                "TRBV28",
                "TRAV38-1",
                "CREB3L3",
                "TRAV8-1",
                "TRBV13",
                "TRAV20",
                "AC139099.2",
                "TIMP2")

# Top clone heatmap for genes that differ between top clones and others ----

heatmap_genes <- coi_v_other_cl2_ri$gene
heatmap_genes <- heatmap_genes[! (grepl("^TR[AB]", heatmap_genes) |
                                    heatmap_genes %in% rm_markers)]

agg_dat <- GetAssayData(pseudo_cl2,
                        layer = "scale.data")[heatmap_genes, ]

column_ha = HeatmapAnnotation(Clone = pseudo_cl2$coi,
                              Sample = pseudo_cl2$Sample)

pdf(file.path(fig_dir, "heatmap_cl2_clones_scale_data_no_tcr.pdf"),
    width = 12, height = 10)
Heatmap(agg_dat,
        top_annotation = column_ha,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 4),
        show_row_dend = FALSE)

dev.off()


# Heatmap for cluster 2 top clones, with and without TOX ---- 

agg_dat <- GetAssayData(pseudo_cl2,
                        layer = "scale.data")[markers, ]

column_ha = HeatmapAnnotation(Clone = pseudo_cl2$coi,
                              Sample = pseudo_cl2$Sample)

pdf(file.path(fig_dir, "heatmap_cl2_clones_scale_data_markers.pdf"),
    width = 12, height = 5)
Heatmap(agg_dat,
        top_annotation = column_ha,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 4),
        show_row_dend = FALSE)

dev.off()


# Remove TOX

agg_dat <- agg_dat[setdiff(markers, "TOX"), ]

pdf(file.path(fig_dir, "heatmap_cl2_clones_scale_data_markers_no_TOX.pdf"),
    width = 12, height = 5)
Heatmap(agg_dat,
        top_annotation = column_ha,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 4),
        show_row_dend = FALSE)

dev.off()

# Heatmap for top clones, Ri versus healthy ----

ri_v_healthy <- ri_v_healthy %>%
    dplyr::rowwise() %>%
    dplyr::mutate(max_pct = max(`pct.1`, `pct.2`)) %>%
    dplyr::filter(abs(avg_log2FC) >= 1,
                  p_val_adj <= 0.001,
                  max_pct >= 0.25)

agg_dat <- GetAssayData(pseudo_cl2,
                        layer = "scale.data")[ri_v_healthy$gene, ]

pseudo_cl2$Condition <- gsub("[0-9]+(-dis|-5m)?$", "", pseudo_cl2$Sample)

column_ha = HeatmapAnnotation(Clone = pseudo_cl2$coi,
                              Sample = pseudo_cl2$Sample,
                              Condition = pseudo_cl2$Condition)

pdf(file.path(fig_dir, "heatmap_cl2_clones_scale_data_Ri_v_healthy.pdf"),
    width = 12, height = 10)
Heatmap(agg_dat,
        top_annotation = column_ha,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 4),
        show_row_dend = FALSE)

dev.off()