# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("tidyverse")
library("scales")
library("Seurat")
library("ComplexHeatmap")
library("khroma")
library("viridis")
library("RColorBrewer")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--clones', '-c',
                    help = 'Tested clones, rds')
parser$add_argument('--results',  '-f', 
                    help = 'Directory for saving results')
parser$add_argument('--workdir',  '-w', 
                    help = 'working directory')

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/funcs_custom_heatmaps.R"))
source(file.path(args$workdir, "scripts/markers_sp1.R"))

# ----------------------------------------------------------------------------
# Functions ----

# palettes ----
palettes <- function(results){
    palettes <- list("viridis" = viridis(100))
    names(palettes) <- file.path(results, names(palettes))
    
    dummy <- lapply(names(palettes), function(nm){
        if (! file.exists(nm)) { dir.create(nm) }
    })

    return(palettes)
}


# make volcano ----
make_volcano <- function(de, res_dir){
    pdf(file.path(res_dir, "volcano.pdf"))
    p <- volc_mod(de,
                  xlab = expression("Average " * log[2] * "FC"),
                  ylab = expression(-log[10]*"(adjusted p-value)"),
                  x = "avg_log2FC",
                  y = "p_val_adj")
    print(p)
    dev.off()
}


de_reactive_Ri_v_HD <- function(seurat_obj){
    # Subset to just reactive clones, test HD v Ri ----
    rx <- subset(clones_wo_5m, reactive == "TRUE")
    cnd_de <- run_de(rx, "condition", "Ri", "HD", "reactive_Ri_v_HD")
    cnd_de_pb <- run_pseudo(rx, "condition", "Ri", "HD", "reactive_Ri_v_HD_pseudo")
    
}


# main ----
main <- function(args){
    
    results <- file.path(args$results, "heatmaps_fig_4_5")
    if (! file.exists(results)) { dir.create(results, recursive = TRUE) }

    palettes <- palettes(results) 
    
    # Set up for cells with tested reactivity ----
    
    clones <- read_rds(args$clones)
    
    clones[[]] <- clones[[]] %>%
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
                                        "Ri non-reactive"),
                      rx_by_cnd = factor(rx_by_cnd,
                                         levels = c("Ri reactive",
                                                    "HD reactive",
                                                    "Ri non-reactive",
                                                    "HD non-reactive")))
    
    # Heatmaps for curated marker sets ----
    markers <- fig_5_markers()
    rx_partial_small <- purrr::partial(rx_marker_set,
                                       all_clones = clones,
                                       palettes = palettes,
                                       width = 7, height = 7,
                                       row_names_gp = gpar(fontsize = 14),
                                       row_names_side = "left")
    
    dummy <- lapply(names(markers), function(name){
        print_nm <- gsub("[[:punct:][:space:]]", "_", name)
        rx_partial_small(markers = data.frame(gene = markers[[name]],
                                              cat_label = name),
                         name = sprintf("%s.pdf", print_nm))
    })
    
    dummy <- lapply(names(markers), function(name){
        print_nm <- gsub("[[:punct:][:space:]]", "_", name)
        rx_partial_small(markers = data.frame(gene = markers[[name]],
                                              cat_label = name),
                         name = sprintf("%s_average.pdf", print_nm),
                         agg_method = "average")
    })   
    
    # ____________
    # TO DO -----
    # Make pseudobulk by clone, average across
    
    # ___________
}

# ----------------------------------------------------------------------------
main(args)
