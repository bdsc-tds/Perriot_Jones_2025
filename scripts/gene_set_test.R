library("DESeq2")
library("edgeR")
library("limma")
library("tidyverse")
library("PCAtools")

source(file.path(".", "scripts/markers_sp1.R"))
source(file.path(".", "scripts/funcs_custom_heatmaps.R"))

fig5_markers <- fig_5_markers()
all_markers <- unname(unlist(fig5_markers))
marker_df <- tibble(gene = unlist(fig5_markers),
                    category = rep(names(fig5_markers), lengths(fig5_markers)))

# make heatmap ----
make_heatmap <- function(dat,
                         col_group,
                         col_labs = NULL,
                         disp_min = -2.5,
                         disp_max = 2.5,
                         row_group = FALSE,
                         palette = viridis::viridis(100),
                         column_title_gp = gpar(fontsize = 10),
                         column_names_gp = gpar(fontsize = 10),
                         row_names_gp = gpar(fontsize = 3),
                         ...){
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
            col = palette,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            ...)
    
}


make_deseq_heatmap <- function(pb, samples, all_markers = all_markers){
    
    dds1 <- DESeq2::DESeqDataSetFromMatrix(countData = pb[keep_genes, ], 
                                           colData = samples,
                                           design = ~condition)
    
    dds1 <- DESeq2::estimateSizeFactors(object = dds1)
    dds1 <- DESeq2::estimateDispersions(object = dds1, fitType = "local")
    
    fpms <- fpm(dds1)
    fpms <- fpms[intersect(rownames(fpms), all_markers), ]
    cats <- marker_df[match(rownames(fpms), marker_df$gene), ]$category

    
    pdf(file.path(out_dir, "deseq2_log_fpm_plus_1.pdf"), height = 12)
    h <- make_heatmap(log2(fpms + 1),
                      col_group = samples$condition,
                      row_names_gp = gpar(fontsize = 6),
                      row_group = cats,
                      row_split = length(unique(cats)),
                      row_title_gp = gpar(fontsize = 7.4))
    print(h)
    dev.off()
    
    #dds1 <- DESeq2::nbinomWaldTest(object = dds1)
    #res <- DESeq2::results(object = dds1,
    #                       contrast = c("condition", "Ri", "HD"),
    #                       alpha = 0.05)
    
    #to.return <- data.frame(p_val = res$pvalue,
    #                        row.names = rownames(res))
    #return(to.return)
    
}

# Analysis ----

base_dir <- "reactive_cl1/"
cd3 <- file.path(base_dir, "reactive_Ri_v_HD_pseudo_cdr3")
clone <- file.path(base_dir, "reactive_Ri_v_HD_pseudo_clone")
smp <- file.path(base_dir, "reactive_Ri_v_HD_pseudo_sample")

