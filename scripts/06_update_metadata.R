# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', help = 'Seurat object')
parser$add_argument('--clones', help = 'csv of tested clones')

args <- parser$parse_args()

source(file.path(args$workdir, "scripts/clones_of_interest.R"))

# ----------------------------------------------------------------------------


# Get the cluster with the most clones of interest ----
get_cluster_coi <- function(md){
    md %>%
        dplyr::filter(is_coi == "coi") %>%
        dplyr::group_by(seurat_clusters) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::arrange(desc(n)) %>%
        dplyr::slice_head(n = 1) %>%
        dplyr::pull(seurat_clusters)
}

# Subset to reactive clones from the cluster of interest ----
subset_rx_coi_dis <- function(obj){
    
    
}

# Subset to the cluster of interest ----
subset_coi_dis <- function(obj){
    
    
}

# main ----
main <- function(args){
    seurat_obj <- read_rds(args$seurat)
    selected_clones <- read_csv(args$clones) %>% unique()
    
    # Add tested clone names into metadata 
    # (This works as data was filtered for exactly one tcr beta,
    # otherwise it would be necessary to first separate multiple chains)
    seurat_obj[[]] <- seurat_obj[[]] %>%
        dplyr::mutate(cdr3_aa_beta = CTaa2) %>%
        tidyr::separate(TCR2,
                        into = c("trbv", "trbd", "trbj", "trbc"),
                        sep = "\\.",
                        remove = FALSE) %>%
        dplyr::left_join(selected_clones,
                         by = c("Sample","cdr3_aa_beta", "trbv", "trbj"),
                         relationship = "many-to-one")
    
    
    # Add reactivity x condition ----
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
                                        "Ri non-reactive"),
                      reactive_lab = case_when(reactive == 
                                                   "TRUE" ~ "Reactive",
                                               reactive ==
                                                   "FALSE" ~ "Non-reactive"))
    
    # Subset to cluster of interest, save
    coi_cluster_id <- get_cluster_coi(seurat_obj[[]])
    coi <- subset(seurat_obj, seurat_clusters == coi_cluster_id)
    write_rds(coi, file.path(dirname(args$seurat), "cluster_of_interest.rds"))
    
    # Subset to clones where reactivity information is known ----
    clones <- subset(seurat_obj, reactive %in% c(TRUE, FALSE))
    write_rds(clones, file.path(dirname(args$seurat), "reactivity_tested.rds"))
    
    # Remove followup sample, adjust condition ---- 
    clones_wo_5m <- subset(clones, Sample != "Ri01_5m")
    
    clones_wo_5m[[]] <- clones_wo_5m[[]] %>%
        dplyr::mutate(condition = case_when(condition == "Ri01_dis" ~ "Ri",
                                            TRUE ~ condition)) %>%
        dplyr::mutate(rx_by_cnd = paste(condition, reactive, sep = "_"))
    Idents(clones_wo_5m) <- "rx_by_cnd"
    
    # Subset to just reactive clones ----
    rx <- subset(clones_wo_5m, reactive == "TRUE")
    write_rds(clones, file.path(dirname(args$seurat), "reactive_clones.rds"))
    
    # Subset to just reactive clones from cluster of interest ----
    rx <- subset(clones_wo_5m,
                 reactive == "TRUE",
                 seurat_clusters == coi_cluster_id)
    write_rds(clones, file.path(dirname(args$seurat),
                                "reactive_cluster_of_interest.rds"))
    
    # Re-save seurat object and metadata
    write_rds(seurat_obj, args$seurat)
    write_csv(seurat_obj[[]], gsub("rds$", "csv.gz", args$seurat))

}


