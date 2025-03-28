# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("Seurat")

# Command line arguments ----
parser <- ArgumentParser(description = "Differential expression analyses")

parser$add_argument('--seurat', help = 'Seurat object')
parser$add_argument('--clones', help = 'csv of tested clones')

args <- parser$parse_args()


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
    
    # Add reactivity x condition
    # TO DO - check if Ri01_5m is a problem here
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
    
    # Re-save seurat object and metadata
    write_rds(seurat_obj, args$seurat)
    write_csv(seurat_obj[[]], gsub("rds$", "csv.gz", args$seurat))

}


