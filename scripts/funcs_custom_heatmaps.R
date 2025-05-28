library("ComplexHeatmap")

# pseudobulk heatmap (for figure 4f and 5c-h) ----
pb_heatmap <- function(pseudo, markers, palette, ...) {
    anno <- markers$cat_label[match(Features(pseudo), markers$gene)]
    cats <- unique(anno) 
    dat <- LayerData(pseudo, "scale.data")
    
    if (length(cats) == 1){
        cluster_rows <- TRUE
        row_split <- NULL
        row_ha <- NULL
    } else {
        cluster_rows <- cluster_within_group(t(dat), anno)
        row_split <- length(cats)
        row_ha <- rowAnnotation(Category = anno,
                                col = list(Category = 
                                               structure(color("light")(length(cats)),
                                                         names = cats)),
                                show_legend = c(FALSE),
                                show_annotation_name = FALSE)
    }
    
    heatmap_args <- list(matrix = dat,
                         col = palette,
                         cluster_columns = TRUE,
                         show_column_dend = FALSE,
                         show_row_dend = FALSE, 
                         column_title_gp = gpar(fontsize = 10),
                         row_names_gp = gpar(fontsize = 7),
                         column_names_gp = gpar(fontsize = 7),
                         row_split = row_split,
                         row_title_gp = gpar(fontsize = 7.4),
                         column_labels = unique(pseudo$orig.ident),
                         heatmap_legend_param = list(title = "Scaled \nexpression"),
                         left_annotation = row_ha,
                         cluster_rows = cluster_rows,
                         column_names_rot = 45,
                         row_gap = unit(2, "mm"))
    
    heatmap_args <- modifyList(heatmap_args, list(...))
    
    return(do.call(Heatmap, heatmap_args))
}

# reactivity analyses for single marker set ----
pb_marker_set <- function(all_clones,
                          markers,
                          palettes,
                          group_by = "rx_by_cnd",
                          name = "category_by_reactivity.pdf",
                          width = 3,
                          height = 8,
                          agg_method = "pseudobulk",
                          ...){
    
    # Subset to genes of interest, aggregate expression
    clones <- subset(all_clones,
                     Sample != "Ri01_5m",
                     features = markers$gene)
    
    if (agg_method == "pseudobulk"){
        pseudo_cat <- AggregateExpression(clones,
                                          group.by = group_by,
                                          return.seurat = TRUE)
    } else {
        pseudo_cat <- AverageExpression(clones,
                                        group.by = group_by,
                                        return.seurat = TRUE)
    }
    
    
    if (group_by == "seurat_clusters"){
        pseudo_cat$orig.ident <- gsub("g", "", pseudo_cat$orig.ident)
    }
    
    # For each palette, make a heatmap of aggregated values
    dummy <- lapply(names(palettes), function(pal_dir){
        print(file.path(pal_dir, name))
        pdf(file.path(pal_dir, name), width = width, height = height)
        h <- pb_heatmap(pseudo_cat, markers, palette = palettes[[pal_dir]], ...)
        print(h)
        dev.off()
    })
}


palette_steps <- function(disp.min, disp.max, n_steps = 49){
    steps <- (disp.max - disp.min)/n_steps
    steps <- seq(disp.min, disp.max, by = steps) 
    return(steps)
}


# purple_and_yellow ----
pp_and_yl <- function(disp.min = -2.5, disp.max = 2.5){
    steps <- palette_steps(disp.min = disp.min, disp.max = disp.min)
    pp_yl_pal <- circlize::colorRamp2(steps, PurpleAndYellow())
    return(pp_yl_pal)
}

# heatmap_w_labs ----
heatmap_w_labs <- function(obj,
                           col_group,
                           col_labs = NULL,
                           disp_min = -2.5,
                           disp_max = 2.5,
                           row_group = FALSE,
                           palette = NULL,
                           column_title_gp = gpar(fontsize = 10),
                           column_names_gp = gpar(fontsize = 10),
                           row_names_gp = gpar(fontsize = 3),
                           ...){
    if (is.null(palette)){
        #palette <- pp_and_yl(disp.min = disp_min, disp.max = disp_max) 
        palette <- viridis::viridis(100)
    }
    
    dat <- LayerData(obj, "scale.data")
    
    if (is.null(col_labs)) {
        temp <- unique(as.character(obj[[col_group]][,1]))
        col_labs <- structure(temp, names = temp)
    }
    
    col_group <- col_labs[as.character(obj[[col_group]][,1])]
    
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