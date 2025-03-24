# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ComplexHeatmap")
library("ggplot2")
library("tidyverse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', '-s',
                    help = 'Seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
args <- parser$parse_args()

# ----------------------------------------------------------------------------
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

# purple and yellow palette ----
purple_and_yellow <- function(disp.min = -2.5, disp.max = 2.5){
    pp_yl <- (disp.max - disp.min)/49
    pp_yl <- seq(disp.min, disp.max, by = pp_yl) 
    pp_yl_pal <- circlize::colorRamp2(pp_yl, PurpleAndYellow())
    return(pp_yl_pal)
}

# heatmap w labs ----
heatmap_w_labs <- function(obj,
                           col_group,
                           col_labs = NULL,
                           disp_min = -2.5,
                           disp_max = 2.5,
                           row_group = FALSE,
                           column_title_gp = gpar(fontsize = 10),
                           column_names_gp = gpar(fontsize = 10),
                           row_names_gp = gpar(fontsize = 3),
                           ...){
    print("in heatmap w labs")
    pp_yl_pal <- purple_and_yellow(disp.min = disp_min, disp.max = disp_max) 
    dat <- LayerData(obj, "scale.data")
    
    if (is.null(col_labs)) {
        temp <- unique(as.character(obj[[col_group]][,1]))
        col_labs <- structure(temp, names = temp)
    }
    
    print(col_labs)
    print(head(obj[[col_group]][,1], n = 20))
    
    col_group <- col_labs[as.character(obj[[col_group]][,1])]
    print(head(col_group))
    print(unique(col_group))
    
    if (! isTRUE(row_group) & ! isFALSE(row_group)) { 
        row_group <- cluster_within_group(t(dat), row_group)
    }
    
    Heatmap(dat,
            cluster_columns = cluster_within_group(dat, col_group),
            cluster_rows = row_group,
            show_column_dend = FALSE,
            show_row_dend = FALSE, 
            column_split = length(unique(col_group)),
            column_title_gp = column_title_gp,
            column_names_gp = column_names_gp, 
            row_names_gp = row_names_gp,
            column_labels = gsub("_.*", "", colnames(dat)),
            col = pp_yl_pal,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            ...)
}



# Heatmap of cluster markers ----
make_heatmap <- function(seurat_obj, markers, out_fname,
                         wd = 7, ht = 7, size = 4){
    pdf(out_fname, width = wd, height = ht)
    Idents(seurat_obj) <- seurat_obj$seurat_clusters 
    p <- DoHeatmap(seurat_obj,
                   size = size,
                   features = markers)
    print(p)
    dev.off()
    
}

# Dotplot of cluster markers ----
make_dotplot <- function(seurat_obj, markers, out_fname,
                         ht = 5, scale = TRUE){
    
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
    results <- file.path(args$results, "dotplots_and_heatmaps")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }
    
    seurat_obj <- read_rds(args$seurat)
    marker_sets <- get_markers()
    
    #cl_levs <- as.character(sort(as.numeric(unique(seurat_obj$seurat_clusters))))
    #seurat_obj$seurat_clusters  <- factor(seurat_obj$seurat_clusters,
    #                                      levels = cl_levs)
    
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
        Idents(marker_subset) <- as.numeric(marker_subset$seurat_clusters)

        # All samples together ----
        #make_dotplot(marker_subset, markers,
        #             sprintf(template, "dotplot_unscaled", out_nm),
        #             scale = FALSE)
        #make_dotplot(marker_subset, markers,
        #             sprintf(template, "dotplot_scaled", out_nm),
        #             scale = TRUE)
        #make_heatmap(marker_subset, markers,
        #             sprintf(template, "heatmap", out_nm))
        
        # Only Ri01_dis ----
        ri01_dis <- subset(marker_subset, Sample == "Ri01_dis")
        make_dotplot(ri01_dis, markers,
                     sprintf(template, "Ri01_dis_dotplot_unscaled", out_nm),
                     scale = FALSE)
        make_dotplot(ri01_dis, markers,
                     sprintf(template, "Ri01_dis_dotplot_scaled", out_nm),
                     scale = TRUE)

        # All samples except Ri01_5m -----
        not_ri01_5m <- subset(marker_subset, Sample != "Ri01_5m")
        make_dotplot(not_ri01_5m, markers,
                     sprintf(template,
                             "not_Ri01_5m_dotplot_unscaled", out_nm),
                     scale = FALSE)
        make_dotplot(not_ri01_5m, markers,
                     sprintf(template,
                             "not_Ri01_5m_dotplot_scaled", out_nm),
                     scale = TRUE)
        
        # Clustered heatmaps, only Ri01_5m
        tryCatch({
            ri01_5m <- subset(marker_subset, Sample == "Ri01_5m")
            print(dim(ri01_5m))
            pdf(sprintf(template, "heatmap_ri01_5m", out_nm),
                width = 10)
            h <- heatmap_w_labs(ri01_5m,
                                col_group = "seurat_clusters",
                                row_names_gp = gpar(fontsize = 8),
                                show_column_names = FALSE)
            print(h)
            dev.off()
            
        }, error = function(e){return()})
        
        # Clustered heatmaps, all samples except Ri01_5m
        tryCatch({
            dis <- subset(marker_subset, Sample != "Ri01_5m")
            h <- heatmap_w_labs(dis, col_group = "seurat_clusters",
                                row_names_gp = gpar(fontsize = 8),
                                show_column_names = FALSE)
            pdf(sprintf(template, "heatmap_excluding_ri01_5m", out_nm),
                width = 10)
            print(h)
            dev.off()
        }, error = function(e){return()})
    }
}

# ----------------------------------------------------------------------------
main(args)
