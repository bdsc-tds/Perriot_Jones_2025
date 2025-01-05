# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("tidyverse")

parser <- ArgumentParser(description = "Tables for clone of interest")
parser$add_argument('--metadata', '-m',
                    help = 'Metadata from seurat object')
parser$add_argument('--output',  '-o', 
                    help = 'Directory for saving tables')
parser$add_argument('--clones',  '-c', 
                    help = 'R script containing sequences of interest')

args <- parser$parse_args()

# ----------------------------------------------------------------------------
# Functions ----
# Clone composition of clusters containing clones of interest ----
clones_in_coi <- function(meta, out_dir){
    # Select clusters with at least 5 clones of interest,
    # regardless of sample and sequence
    clusters_w_roi = meta %>%
        dplyr::filter(coi == TRUE) %>%
        dplyr::group_by(seurat_clusters) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::filter(n >= 5) %>%
        dplyr::pull(seurat_clusters) %>%
        unique()
    
    clones_in_coi <- meta %>%
        dplyr::filter(seurat_clusters %in% clusters_w_roi) %>%
        dplyr::group_by(seurat_clusters, Sample, coi) %>%
        dplyr::mutate(n_sample_cluster = n()) %>%
        dplyr::group_by(seurat_clusters,
                        Sample,
                        coi,
                        CTaa2,
                        beta_aa,
                        n_sample_cluster) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::filter(! beta_aa == "NA_NA") %>%
        dplyr::arrange(seurat_clusters, desc(n)) %>%
        dplyr::mutate(pct_sample_cluster = n / n_sample_cluster * 100)
    
    write_csv(clones_in_coi, 
              file.path(out_dir, "clones_in_clusters_with_coi.csv"))
}


# Top clones per cluster and donor ----
top_clones <- function(meta, out_dir){
    all_clones <- meta %>%
        
        dplyr::group_by(Sample, seurat_clusters) %>%
        dplyr::mutate(n_sample_cluster = n()) %>%
        
        dplyr::group_by(Sample) %>%
        dplyr::mutate(n_sample = n()) %>%
        
        dplyr::group_by(seurat_clusters) %>%
        dplyr::mutate(n_cluster = n()) %>%
        
        dplyr::group_by(Sample,
                        seurat_clusters,
                        beta_aa,
                        n_sample,
                        n_cluster,
                        n_sample_cluster) %>%
        
        dplyr::summarise(n = n()) %>%
        dplyr::mutate(pct_sample_cluster = n / n_sample_cluster * 100,
                      pct_sample = n / n_sample * 100,
                      pct_cluster = n / n_cluster * 100) %>%
        
        dplyr::arrange(Sample, seurat_clusters, desc(n)) %>%
        dplyr::relocate(Sample, seurat_clusters, beta_aa,
                        n_sample_cluster, pct_sample_cluster,
                        n_sample, pct_sample,
                        n_cluster, pct_cluster)
    
    write_csv(all_clones,
              file.path(out_dir, "all_clones_counts_by_sample_cluster.csv"))
    
    # Filter for clones that account for at least 1% of the sample
    top_clones <- all_clones %>%
        dplyr::filter(pct_sample >= 1) %>%
        write_csv(file.path(out_dir,
                            "top_clones_counts_by_sample_cluster.csv"))
    
}


# Count clone of interest ----
count_coi <- function(md, out_dir){
    # Counts of clone of interest per sample
    md %>%
        dplyr::filter(coi == TRUE) %>%
        dplyr::select(Sample, seurat_clusters, beta_aa) %>%
        dplyr::group_by(across(everything())) %>%
        dplyr::summarise(n = n()) %>%
        dplyr::arrange(desc(n)) %>%
        write_csv(file.path(out_dir, "clone_of_interest_counts.csv"))
}

# make_tables ----
make_tables <- function(args){
    if (! dir.exists(args$output)) { dir.create(args$output) }
    
    source(args$clone)
    coi <- get_coi()
    
    md <- read_csv(args$metadata) %>%
        dplyr::mutate(coi = beta_aa %in% coi)
    
    count_coi(md, args$output)
    top_clones(md, args$output)
    clones_in_coi(md, args$output)
} 

# ----------------------------------------------------------------------------
make_tables(args)
    
    