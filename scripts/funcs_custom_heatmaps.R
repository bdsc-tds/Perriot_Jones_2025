library("ComplexHeatmap")

# purple_and_yellow ----
purple_and_yellow <- function(disp.min = -2.5, disp.max = 2.5){
    pp_yl <- (disp.max - disp.min)/49
    pp_yl <- seq(disp.min, disp.max, by = pp_yl) 
    pp_yl_pal <- circlize::colorRamp2(pp_yl, PurpleAndYellow())
    return(pp_yl_pal)
}


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
            col = pp_yl_pal,
            heatmap_legend_param = list(title = "Scaled \nexpression"),
            ...)
}