# ----------------------------------------------------------------------------
# Libraries and setup ----

library("argparse")
library("ggplot2")
library("khroma") 
library("patchwork")
library("scales")
library("tidyverse")
library("viridis")

parser <- ArgumentParser(description = "Tables for clone of interest")
parser$add_argument('--metadata', '-m',
                    help = 'Metadata from seurat object')
parser$add_argument('--figures',  '-f', 
                    help = 'Directory for saving figures')
parser$add_argument('--clones',  '-c', 
                    help = 'R script containing sequences of interest')

args <- parser$parse_args()

# ----------------------------------------------------------------------------
# Functions ----
# umap_by_cluster ----
umap_by_cluster <- function(meta, fig_dir){
    
    # Get the default colours used in the seurat generated UMAP
    cl_names <- as.character(sort(unique(meta$seurat_clusters)))
    n_clusters <- length(cl_names)
    default_pal <- structure(scales::hue_pal()(n_clusters),
                             names = cl_names)
    
    cl_umaps <- lapply(cl_names, function(cl){
        meta_cl <- meta %>%
            dplyr::mutate(is_cl = seurat_clusters == cl) %>%
            dplyr::arrange(is_cl)
        
        pl <- ggplot(meta_cl, 
                     aes(x = UMAPfull_1, y = UMAPfull_2)) +
            geom_hex(aes(fill = is_cl)) + 
            theme_minimal(base_size = 15) +
            scale_fill_manual(labels = c(FALSE, TRUE),
                              values = c("lightgray", default_pal[[cl]]),
                              guide = NULL) +
            labs(title = cl,
                 x = NULL,
                 y = NULL) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
        pl
    })
    
    pdf(file.path(fig_dir, "umap_split_by_clusters.pdf"),
        width = 12, height = 13)
    p <- wrap_plots(cl_umaps, ncol = 4)
    print(p)
    dev.off()
}

# umap_coi ----
umap_coi <- function(meta, fig_dir){
    
    # UMAP with clone of interest highlighted, separated by sample ----
    
    meta$is_coi = c(2, 0.5)[as.numeric(meta$coi == FALSE) + 1]
    meta <- meta %>%
        dplyr::mutate(coi = factor(coi, levels = c("no", "Ri01", "Ri02"))) %>%
        dplyr::arrange(coi) %>%
        dplyr::mutate(coi = as.character(coi))
    
    coi_colours <- structure(c("lightgray","#DC050C","steelblue"),
                             names = c("no", "Ri01", "Ri02"))
    
    pdf(file.path(fig_dir, "umap_coi_by_sample_manual.pdf"),
        height = 12, width = 12)
    p <- ggplot(meta, 
           aes(x = UMAPfull_1, y = UMAPfull_2)) +
        geom_point(aes(size = is_coi, colour = coi)) +
        scale_size_identity() +
        facet_wrap(~Sample) +
        theme_minimal(base_size = 15) +
        scale_color_manual(labels = names(coi_colours), values = coi_colours) +
        labs(x = "UMAP dimension 1", y = "UMAP dimension 2") +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    print(p)
    dev.off()
}

# make_custom_umaps ----
make_custom_umaps <- function(args){
    if (! dir.exists(args$figures)) { dir.create(args$figures) }

    source(args$clone)
    coi <- get_coi()
    
    md <- read_csv(args$metadata) %>%
        dplyr::mutate(coi = case_when(beta_aa == coi[["ri01_coi"]] ~ "Ri01",
                                      beta_aa == coi[["ri02_coi"]] ~ "Ri02",
                                      TRUE ~ "no"))
    
    umap_by_cluster(meta, args$figures)
    umap_coi(meta, fig_dir)
}

# ----------------------------------------------------------------------------

make_custom_umaps(args)
