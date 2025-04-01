library("ComplexHeatmap")

palette_steps <- function(disp.min, disp.max, n_steps = 49){
    steps <- (disp.max - disp.min)/n_steps
    steps <- seq(disp.min, disp.max, by = steps) 
    return(steps)
}

# purple_and_yellow ----
purple_and_yellow <- function(disp.min = -2.5, disp.max = 2.5){
    steps <- palette_steps(disp.min = disp.min, disp.max = disp.min)
    pp_yl_pal <- circlize::colorRamp2(steps, PurpleAndYellow())
    return(pp_yl_pal)
}

blue_and_yellow <- function(disp.min = -2.5, disp.max = 2.5){
    bu_yl <- c("#0C7BDC", "#2E73CA", "#7E6C6E", "#B28808", "#FFC20A")
    steps <- palette_steps(disp.min, disp.max, n_steps = 4)
    bu_yl_pal <- circlize::colorRamp2(steps, bu_yl)
    return(bu_yl_pal)
}

blue_white_yellow <- function(disp.min = -2.5, disp.max = 2.5){
    bu_yl <- c("#0C7BDC", "#2E73CA", "#FFFFFF", "#B28808", "#FFC20A")
    steps <- palette_steps(disp.min, disp.max, n_steps = 4)
    bu_yl_pal <- circlize::colorRamp2(steps, bu_yl)
    return(bu_yl_pal)
}


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
        palette <- purple_and_yellow(disp.min = disp_min, disp.max = disp_max) 
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