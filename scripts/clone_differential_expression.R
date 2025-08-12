# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ComplexHeatmap")
library("DESeq2")
library("ggplot2")
library("tidyverse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', '-s',
                    help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
parser$add_argument('--clones', help = 'Clone table')
parser$add_argument('--workdir',
                    help = "Working directory, for loading scripts")
parser$add_argument('--min_cells',  '-m', default = 5, 
                    help = 'Minimum cells per clone for pseudobulk')

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/funcs_EnhancedVolcano_mod.R"))
source(file.path(args$workdir, "scripts/funcs_differential_expression.R"))
source(file.path(args$workdir, "scripts/clones_of_interest.R"))
# ----------------------------------------------------------------------------

# filter_tcr_genes ----
filter_tcr_genes <- function(markers){
    markers %>%
        dplyr::filter(! grepl("^TR[ABGD][VDJ]", gene))
} 


# markers ----
markers <- function(){
    return(c("IKZF2", "KLRC2", "KLRC3", "KLRK1", "KIR2DL1", "KIR2DL3", "IL7R",
             "FOS", "JUN", "RELB", "EGR1", "EGR2", #"EGR3", "ZAP70", 
             "JUND", "FOSB", "INPP5D", "CBLB",
             "HLA-DQA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1", "HLA-DRB5",
             "CD69", "CD74", "CD44", "VCAM1", "ICAM1",
             "TNF", "IFNG", "GZMA", "GZMH", "CCL4", "CCL3", "CSF1",
             "ZFP36L1", "NFKBIA",
             "NKG7", "ID2"))
}

# Volcano with markers labelled ----
marker_volcano <- function(de, out_fname,
                           labels = markers(),
                           FCcutoff = 0.5, ...){
    pdf(out_fname)
    p <- volc_mod(de,
                  labels = labels,
                  labelsFrom = "gene",
                  labSize = 2,
                  force = 5,
                  FCcutoff = FCcutoff,
                  col = c("grey50", "darkorange", "royalblue")) +
        theme(axis.title = element_text(size = 30),
              axis.text = element_text(size = 24))
    print(p)
    dev.off()
}


# Run differential expression analyses ----
run_diff_expr <- function(obj, results_dir, idents, ident_1, ident_2,
                          name, max_n = 100,
                          pval_min = 0.05, wd = 7, ht = 7, ...){
    res_dir <- file.path(results_dir, name)
    if (! file.exists(res_dir)) { dir.create(res_dir, recursive = TRUE) }
    
    obj <- ScaleData(obj, features = Features(obj))
    
    # Run diff expr
    Idents(obj) <- obj[[idents]][, 1]
    de <- FindMarkers(obj, ident.1 = ident_1, ident.2 = ident_2, ...) %>%
        dplyr::as_tibble(rownames = "gene") %>%
        dplyr::mutate("comparison" = name) %>%
        relocate(gene, comparison)
    
    # CSV of results, tcr_genes removed
    de <- filter_tcr_genes(de)
    write_csv(de, file.path(res_dir, "supp_table_7_diff_expr_rm_tcr.csv"))
    
    # Volcano of results
    make_volcano(de, file.path(res_dir, "fig_5a.pdf"))

    return(de)
}


# main ----
main <- function(args, min_cells = 5, ...){
    
    # Setup ---- 
    if (! file.exists(args$results)) { 
        dir.create(args$results, recursive = TRUE)
        dir.create(file.path(args$results, "condition_by_reactivity"))
        dir.create(file.path(args$results, "tables"))
    }
    
    seurat_obj <- read_rds(args$seurat)
    selected_clones <- read_csv(args$clones) %>% unique()
    
    # Write table of counts per cluster ---- 
    cluster_counts_table(seurat_obj[[]] %>%
                             dplyr::filter(seurat_clusters == coi_cluster_id),
                      file.path(args$results,
                                sprintf("tables/clone_counts_cl_%s.csv",
                                        coi_cluster_id)))
    
    run_de <- purrr::partial(run_diff_expr, results = args$results, ...)
    run_pseudo <- purrr::partial(run_diff_expr_pb, results = args$results, ...)
    
    # Subset to clones where reactivity information is known ----
    clones <- subset(seurat_obj, reactive %in% c(TRUE, FALSE))
    write_rds(clones, file.path(dirname(args$seurat), "reactive_clones.rds"))
    
    # Write table of reactive counts ----
    clones_counts_table(clones[[]],
                        file.path(args$results,
                                  "tables/reactive_clone_counts.csv"))
    
    # Write a table of cluster counts for clones from cluster of interest ----
    cl_clusters <- seurat_obj[[]] %>%
        dplyr::inner_join(cl_clones, relationship = "many-to-one") %>%
        dplyr::group_by(Sample, beta_aa, tcr_name) %>%
        dplyr::mutate(n_clone = n()) %>%
        dplyr::group_by(Sample, beta_aa, tcr_name, seurat_clusters, n_clone) %>%
        dplyr::summarise(n_clone_cluster = n()) %>%
        dplyr::arrange(desc(n_clone), Sample, beta_aa)
    
    write_csv(cl_clusters,
              file.path(args$results,
                        sprintf("tables/clusters_with_cells_in_cl_%s.csv",
                                coi_cluster_id)))
    
    write_csv(cl_clusters %>%
                   tidyr::pivot_wider(names_from = seurat_clusters,
                                      values_from = n_clone_cluster),
              file.path(args$results,
                        sprintf("tables/clusters_with_cells_in_cl_%s_wide.csv",
                                coi_cluster_id)))
    
    
    # Remove followup sample, adjust condition ---- 
    clones_wo_5m <- subset(clones, Sample != "Ri01_5m")
    
    clones_wo_5m[[]] <- clones_wo_5m[[]] %>%
        dplyr::mutate(condition = case_when(condition == "Ri01_dis" ~ "Ri",
                                            TRUE ~ condition),
                      reactive_lab = case_when(reactive == "TRUE" ~ "Reactive",
                                               reactive == "FALSE" ~ "Non-reactive"))
    
    clones_wo_5m[[]] <- clones_wo_5m[[]] %>%
        dplyr::mutate(rx_by_cnd = paste(condition, reactive, sep = "_"))
    Idents(clones_wo_5m) <- "rx_by_cnd"
    
    #_______________________________________________________
    # Violin plots and heatmaps for cluster 1 reactive 
    
    clones_cl1 <- subset(clones_wo_5m,
                         seurat_clusters == "1" & reactive == "TRUE")
    
    markers <- fig_5f_all()
    
    rx_cl1_markers <- subset(clones_cl1, features = markers)
    rx_cl1_markers <- ScaleData(rx_cl1_markers,
                                features = Features(rx_cl1_markers))
    
    pdf("rx_cl1.pdf")
    h <- DoHeatmap(rx_cl1_markers, features = Features(rx_cl1_markers))
    print(h)
    dev.off()
    
    
    pdf("rx_cl1_by_donor.pdf", width = 12)
    h <- heatmap_w_labs(rx_cl1_markers, col_group = "Sample")
    print(h)
    dev.off()
    
    # Make violins ----
    
    # For just cluster 1
    source(file.path(args$workdir, "scripts/markers_sp1.R"))
    fig5_markers <- fig_5_markers()
    
    make_violins(rx_cl1_markers, args, outdir = ".", markers = fig5_markers)
    
    # For all reactive 
    all_rx <- subset(clones_wo_5m,
                     reactive == "TRUE",
                     features = unlist(markers))
    all_rx <- ScaleData(all_rx,
                        features = Features(all_rx))
    
    make_violins(all_rx, args,
                 out_dir = "results/violins/all_reactive/",
                 markers = fig5_markers)
   
    
    # Subset to just reactive clones, test HD v Ri ----
    rx <- subset(clones_cl1, reactive == "TRUE")
    cnd_de <- run_de(rx, "condition", "Ri", "HD", "reactive_Ri_v_HD")
    
    # Pseudobulk at the sample level ----
    cnd_de_pb_sample <- run_pseudo(rx, "condition", "Ri", "HD", 
                                   group_by = "Sample",
                                   name = "reactive_Ri_v_HD_pseudo_sample")
 }

# ----------------------------------------------------------------------------
main(args)
