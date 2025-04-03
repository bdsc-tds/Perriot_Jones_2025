# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ComplexHeatmap")
library("ggplot2")
library("tidyverse")
library("Seurat")
library("PCAtools")

# Command line arguments ----
parser <- ArgumentParser(description = "Volcano plot for figure 5a")

parser$add_argument('--de',  
                    help = 'csv with differential expression results')
parser$add_argument('--workdir',
                    help = "Working directory, for loading scripts")

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/funcs_EnhancedVolcano_mod.R"))

# ----------------------------------------------------------------------------
# main ----
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

selected_markers <- function(){
    return(list(A = c("EGR2","EGR4","TNF","EGR3","JUND",
                      "CCL4","FOSB","CCL3","ICAM1","HLA-DRA",
                      "FOS","IFNG","CSF1","IL7R","KIR2DL3"),
                
                B = c("JUND","CD74","FOSB","HLA-DRA","ZFP36L1",
                    "KIR2DL3","HLA-DRB1","NFKBIA","TNF","CCL4",
                    "FOS","EGR2","IL7R","JUN","RELB")))
}


marker_volcano <- function(de, out_fname,
                           labels = markers(),
                           FCcutoff = 0.5, ...){
    pdf(out_fname)
    p <- volc_mod(de,
                  labels = labels,
                  labelsFrom = "gene",
                  labSize = 5.5,
                  force = 50,
                  FCcutoff = FCcutoff,
                  col = c("grey50", "#E67E22", "#0C7BDC"),
                  ...) +
        theme(axis.title = element_text(size = 30),
              axis.text = element_text(size = 24))
    print(p)
    dev.off()
}

main <- function(args){
    res_dir <- "results/one_tcr_beta/figures_in_paper/volcano_fig_5/" 
    coi_res <- "results/one_tcr_beta/figures_in_paper/clones_of_interest" 
    de_template <- "results/one_tcr_beta/clone_analyses/%s/diff_expr_full.csv"
    
    # Volcano fig 5a ----
    
    de_cond <- read_csv(sprintf(de_template, "reactive_Ri_v_HD"))
    marker_volcano(de_cond,
                   file.path(res_dir,
                             "volcano_reactive_Ri_v_HD_blue_orange.pdf"))
    
    # Volcano, clones of interest ----
    de_fname <- sprintf(de_template, "Ri01_dis_v_5m_clones_of_interest")
    de_coi <- read_csv(de_fname)
    marker_volcano(de_coi,
                   file.path(coi_res,
                             "volcano_coi_Ri01_dis_v_Ri01_5m_blue_orange.pdf"))
    
    # PCA, clones of interest ----
    ri01_data <- read_rds("data/processed/one_tcr_beta/ri01_cois.rds")
    ri01_data <- ri01_data[markers(), ]
    
    ri01_md <- read_csv("data/processed/one_tcr_beta/ri01_cois_metadata.csv")
    ri01_md <- as.data.frame(ri01_md)
    rownames(ri01_md) <- ri01_md$Cell
    
    ri01_pcs <- pca(ri01_data, metadata = ri01_md)
    
    pdf(file.path(coi_res, "biplot_both_clones.pdf"))
    bp <- biplot(ri01_pcs,
                 showLoadings = TRUE,
                 ntopLoadings = length(markers()),
                 lab = NULL,
                 colby = "Sample",
                 colkey = c("Ri01_dis" = "orange", "Ri01_5m" = "#0C7BDC")) 
    print(bp)
    dev.off()

    pdf(file.path(coi_res, "biplot_both_clones_scaled.pdf"))
    ri01_pcs_scaled <- pca(ri01_data, metadata = ri01_md, scale = TRUE)
    bp <- biplot(ri01_pcs_scaled,
                 showLoadings = TRUE,
                 ntopLoadings = length(markers()),
                 lab = NULL,
                 colby = "Sample",
                 colkey = c("Ri01_dis" = "orange", "Ri01_5m" = "#0C7BDC")) 
    print(bp)
    dev.off()
    
    # PCA, selected markers, one plot per PCA ----
    selected_markers <- selected_markers()
    
    ri01_data <- read_rds("data/processed/one_tcr_beta/ri01_cois.rds")
    ri01_a <- ri01_data[selected_markers[["A"]], ]
    ri01_b <- ri01_data[selected_markers[["B"]], ]
    
    ri01_md <- read_csv("data/processed/one_tcr_beta/ri01_cois_metadata.csv")
    ri01_md <- as.data.frame(ri01_md)
    rownames(ri01_md) <- ri01_md$Cell
    
    ri01_a_pcs <- pca(ri01_a, metadata = ri01_md)
    ri01_b_pcs <- pca(ri01_b, metadata = ri01_md)
    ri01_a_pcs_scaled <- pca(ri01_a, metadata = ri01_md, scale = TRUE)
    ri01_b_pcs_scaled <- pca(ri01_b, metadata = ri01_md, scale = TRUE)
    
    # Unscaled ----
    pdf(file.path(coi_res, "pca_both_clones_set_A.pdf"))
    bp <- biplot(ri01_a_pcs,
                 showLoadings = FALSE,
                 lab = NULL,
                 colby = "Sample",
                 colkey = c("Ri01_dis" = "orange", "Ri01_5m" = "#0C7BDC")) 
    print(bp)
    dev.off()
    
    pdf(file.path(coi_res, "pca_both_clones_set_B.pdf"))
    bp <- biplot(ri01_b_pcs,
                 showLoadings = FALSE,
                 lab = NULL,
                 colby = "Sample",
                 colkey = c("Ri01_dis" = "orange", "Ri01_5m" = "#0C7BDC")) 
    print(bp)
    dev.off()
    
    # Scaled ----
    pdf(file.path(coi_res, "pca_both_clones_set_A_scaled.pdf"))
    bp <- biplot(ri01_a_pcs_scaled,
                 showLoadings = FALSE,
                 lab = NULL,
                 colby = "Sample",
                 colkey = c("Ri01_dis" = "orange", "Ri01_5m" = "#0C7BDC")) 
    print(bp)
    dev.off()
    
    pdf(file.path(coi_res, "pca_both_clones_set_B_scaled.pdf"))
    bp <- biplot(ri01_b_pcs_scaled,
                 showLoadings = FALSE,
                 lab = NULL,
                 colby = "Sample",
                 colkey = c("Ri01_dis" = "orange", "Ri01_5m" = "#0C7BDC")) 
    print(bp)
    dev.off()
    
    
    
    
    
    # ------
    
    pdf(file.path(coi_res, "biplot_both_clones.pdf"))
    bp <- biplot(ri01_pcs,
                 showLoadings = TRUE,
                 ntopLoadings = length(markers()),
                 lab = NULL,
                 colby = "Sample",
                 colkey = c("Ri01_dis" = "orange", "Ri01_5m" = "#0C7BDC")) 
    print(bp)
    dev.off()
    
    
    
    
    # -----
    single_clone <- ## REMOVED SELECTION, CHECK IF ALLOWED  
    single_clone_dat <- ri01_data[, single_clone]
    single_clone_md <- ri01_md[single_clone, ]
    
    clone_pcs <- pca(single_clone_dat, metadata = single_clone_md)
    
    pdf(file.path(coi_res, "biplot_first_clone.pdf"))
    bp <- biplot(clone_pcs,
                 showLoadings = TRUE,
                 ntopLoadings = length(markers()),
                 lab = NULL,
                 colby = "Sample",
                 colkey = c("Ri01_dis" = "orange", "Ri01_5m" = "#0C7BDC")) 
    print(bp)
    dev.off()
    
    pdf(file.path(coi_res, "biplot_first_clone_scaled.pdf"))
    clone_pcs_scaled <- pca(single_clone_dat, metadata = single_clone_md, scale = TRUE)
    bp <- biplot(clone_pcs_scaled,
                 showLoadings = TRUE,
                 ntopLoadings = length(markers()),
                 lab = NULL,
                 colby = "Sample",
                 colkey = c("Ri01_dis" = "orange", "Ri01_5m" = "#0C7BDC")) 
    print(bp)
    dev.off()
    

    
    
    
}

# ----------------------------------------------------------------------------
main(args)

# Fig 5a = reactive Ri v HD