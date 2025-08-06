# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
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
    return(c("CCR7", "TCF7", "SELL", "BCL2", "CD27",
             "CXCR3", "IL2RB", "PRF1", "GZMB"))
}


# Dotplot of cluster markers ----
make_dotplot <- function(seurat_obj, markers, out_fname,
                         ht = 5, scale = TRUE){
    
    pdf(out_fname, height = ht)
    p <- DotPlot(seurat_obj,
                 scale = scale,
                 features = Features(seurat_obj)) + 
        coord_flip() +
        labs(y = NULL, x = NULL) +
        scale_color_gradient2(low = "lightgray",
                              mid = "#E1C7C2",
                              high = "#e60000") +
        theme(panel.grid = element_line(color = "gray")) +
        guides(size = guide_legend(title = "Percent\nExpressed"),
               colour = guide_legend(title = "Average\nExpression"))
    print(p)
    dev.off()
    return(p)
}


# main ----
main <- function(args){
    if (! file.exists(args$results)) { dir.create(args$results, recursive = TRUE) }
    
    seurat_obj <- read_rds(args$seurat)
    markers <- get_markers()
        
    marker_subset <- subset(seurat_obj, features = markers)
    DefaultAssay(marker_subset) <- "RNA"
    marker_subset <- ScaleData(marker_subset, features = markers)
    Idents(marker_subset) <- as.numeric(marker_subset$seurat_clusters)
        
    p <- make_dotplot(marker_subset, markers,
                      file.path(args$results, "supp_fig_1b.pdf"), 
                      scale = FALSE)
    plot_dat <- p$data 
    write_csv(plot_dat, file.path(args$results, "supp_fig_1b_data.csv"))
    
}

# ----------------------------------------------------------------------------
main(args)
