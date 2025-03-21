# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("janitor")
library("tidyverse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--clones', help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
parser$add_argument('--workdir',
                    help = "Working directory, for loading scripts")

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/markers_sp1.R"))

# ----------------------------------------------------------------------------

# Marker lists ----
get_markers <- function(){
    list("A_bis" = c("CCR7", "TCF7", "SELL", "BCL2", "CD27",
                     "CXCR3", "IL2RB", "PRF1", "GZMB"),
         "NK-like" = 
             c("KLRC1", "KLRC2", "KLRG1", "KLRK1", "KIR2DL1", "KIR2DL2", "KIR2DL3", 
               "KIR2DL4", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR3DL1", 
               "KIR3DL2", "KIR3DS1", "KIR3DS2"),
         "Activation" = 
             c("CX3CR1", "TFRC", "HLA-DRA", "HLA-DRB1", "IL2RA", "IL2RB", 
               "CD44", "ITK", "NFATC2", "NFKB1", "NME1", "CD27", "CD28",
               "ICOS", "CD40LG", "TNFRSF4", "TNFRSF9"),
         "Effector function/Cytoxicity" = 
             c("FASLG", "IFNG", "TNF", "PRF1", "GZMA", "GZMB", "GZMK", "GZMM"),
         "Exhaustion" = 
             c("PDCD1", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "CD244", "CD160", 
               "TOX", "EOMES", "BATF", "TBX21", "NR4A1", "PRDM1", "ENTPD1"),
         "NK-like" = 
             c("KLRC1", "KLRC2", "KLRG1", "KLRK1", "KIR2DL1", "KIR2DL2", "KIR2DL3", 
               "KIR2DL4", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS4", "KIR3DL1", 
               "KIR3DL2", "KIR3DS1", "KIR3DS2"),
         "Proliferation/survival" = 
             c("MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MYBL2", "MYC", "CDK1", 
               "CDK2", "CDK4", "CDK6", "CCNA1", "CCNB1", "CCND1", "CCND2", "CCND3", 
               "CCNE1", "E2F1", "FOXM1", "BUB1", "TOP2A", "MKI67", "PCNA", "PLK1", 
               "BCL2", "BCL6", "TCF7"),
         "Tissue residency/migration" =
             c("CD69", "ITGA1", "ITGA4", "ITGAE", "CXCR3", "CXCR4", "CCR5", 
               "CXCR6", "CCR7", "S1PR1", "SELL")
    )
}

# Dotplot of cluster markers ----
make_dotplot <- function(seurat_obj, out_fname, ht = 5, scale = TRUE){
    
    pdf(out_fname, height = ht)
    p <- DotPlot(seurat_obj,
                 scale = scale,
                 features = Features(seurat_obj)) + 
        coord_flip() +
        scale_color_viridis_c() +
        labs(y = NULL, x = NULL) +
        theme(panel.grid = element_line(color = "gray")) +
        guides(size = guide_legend(title = "Percent\nExpressed"),
               colour = guide_legend(title = "Average\nExpression"))
    print(p)
    dev.off()
}

# main ----
main <- function(args){
    
    # Setup ----
    results <- file.path(args$results, "dotplots")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    
    seurat_obj <- read_rds(args$clones)
    seurat_obj <- subset(seurat_obj, Sample != "Ri01_5m")
    
    marker_sets <- get_markers()
    
    seurat_obj[[]] <- seurat_obj[[]] %>%
        dplyr::mutate(condition = case_when(condition == "Ri01_dis" ~ "Ri",
                                            TRUE ~ condition),
                      rx_by_cnd =
                          case_when(condition == "HD" & reactive == "TRUE" ~
                                        "HD reactive",
                                    condition == "HD" & reactive == "FALSE" ~
                                        "HD non-reactive",
                                    condition == "Ri" & reactive == "TRUE" ~
                                        "Ri reactive",
                                    condition == "Ri" & reactive == "FALSE" ~
                                        "Ri non-reactive"))
    # Dotplots for marker lists ----
    template <- file.path(results, "%s_%s.pdf")
    
    for (mk_nm in names(marker_sets)){
        print(mk_nm)
        markers <- marker_sets[[mk_nm]]
        if (length(intersect(markers, Features(seurat_obj))) == 0){
            print(sprintf("no markers found: %s", mk_nm))
            next()
        }
        
        out_nm <- gsub("[-\\/ ]", "_", mk_nm)
        marker_subset <- subset(seurat_obj, features = markers)
        DefaultAssay(marker_subset) <- "RNA"
        marker_subset <- ScaleData(marker_subset, features = markers)
        Idents(marker_subset) <- "rx_by_cnd"
    
        make_dotplot(marker_subset, 
                     sprintf(template, "dotplot_unscaled_no_5m", out_nm),
                     scale = FALSE)
        make_dotplot(marker_subset, 
                     sprintf(template, "dotplot_scaled_no_5m", out_nm),
                     scale = TRUE)
    }
    
    # Average expression and percentage of cells expressing list 1_SP ----
    sp <- markers_1_sp()
    sp_obj <- subset(seurat_obj, features = sp)

    # Get percentage and average expression from dotplot data
    Idents(sp_obj) <- "rx_by_cnd"
    sp_obj <- ScaleData(sp_obj, features = sp)
    
    dp_dat <- DotPlot(object = sp_obj,
                      features = sp,
                      group.by = "rx_by_cnd")$data %>%
        dplyr::select(`features.plot`, `avg.exp`, `pct.exp`, id) %>%
        janitor::clean_names() %>%
        tidyr::pivot_wider(id_cols = features_plot,
                           names_from = id,
                           values_from = c(avg_exp, pct_exp)) %>%
        janitor::clean_names() %>%
        dplyr::rename(gene = features_plot) %>%
        readr::write_csv(file.path(args$results, "tables/markers_1_sp_expr.csv"))
}


