# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("tidyverse")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--metadata', '-m',
                    help = 'Metadata from seurat object')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
args <- parser$parse_args()

# ----------------------------------------------------------------------------
# Functions ----


# barplot_data ----
# recreate the data used in a barplot
barplot_data <- function(df, x, fill){
    x_str <- deparse(substitute(x))
    fill_str <- deparse(substitute(fill))
    n_x_str <- paste0("n_", x_str)
    n_fill_str <- paste0("n_", fill_str)
    n_both <- paste(n_x_str, fill_str, sep = "_")
    gp_by <- c(x_str, fill_str, n_x_str, n_fill_str)
    print(gp_by)
    
    df %>%
        dplyr::group_by({{x}}) %>%
        dplyr::mutate("n_{{x}}" := n()) %>%
        dplyr::group_by({{fill}}) %>%
        dplyr::mutate("n_{{fill}}" := n()) %>%
        dplyr::group_by(across(all_of(gp_by))) %>%
        dplyr::summarise("n_{{x}}_{{fill}}" := n()) %>%
        dplyr::mutate("pct_{{x}}" := .data[[n_both]]/.data[[n_x_str]] * 100,
                      "pct_{{fill}}" := 
                          .data[[n_both]]/.data[[n_fill_str]] * 100) %>%
        dplyr::ungroup()
}


# Bar cluster per sample----
make_barplot <- function(md, x, fill, out_fname,
                         ylab = "Proportion of cells", ht = 7, wd = 7){
    
    pdf(out_fname, height = ht, width = wd)
    p <- ggplot(md, aes(x = {{x}}, fill = {{fill}})) +
        geom_bar(position = "fill", color = "black") +
        scale_y_continuous(expand = expansion(c(0,0))) +
        scale_x_discrete(expand = expansion(c(0,0))) +
        coord_flip() + 
        theme_bw() + 
        labs(x = NULL, y = ylab)
    print(p)
    dev.off()
}


# make_barplots ----
make_barplots <- function(md, res_dir){
    # Clusters per sample
    make_barplot(md, x = seurat_clusters, fill = Sample,
                 out_fname = file.path(out_dir, "supp_fig_2b.pdf"))
    barplot_data(md, x = seurat_clusters, fill = Sample) %>%
        write_csv(file.path(res_dir, "supp_fig_2b_data.csv"))
}


# main ----
main <- function(args){
    if (! file.exists(args$results)) { 
        dir.create(args$results, recursive = TRUE)
    }
    
    # Metadata with UMAP coordinates
    md <- read_csv(args$metadata) %>%
        dplyr::mutate(seurat_clusters = 
                          factor(seurat_clusters,
                                 levels = rev(sort(unique(seurat_clusters)))))
    
    make_barplots(md, args$results)
    
    barplot_data(md, x = Sample, fill = seurat_clusters) %>%
        dplyr::mutate(condition = gsub("(HD|Ri).*", "\\1", Sample)) %>%
        write_csv(file.path(args$results, "supp_fig_2c_data.csv"))
    
}

# ----------------------------------------------------------------------------
main(args)
