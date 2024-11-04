# Libraries and setup ----

library("ggplot2")
library("khroma") 
library("patchwork")
library("scales")
library("tidyverse")
library("viridis")

project_dir <- tryCatch({Sys.getenv()[["project_dir"]]},
                        error = function(cond){return(".")})
fig_dir <- file.path(project_dir, "figures")

# Metadata with UMAP coordinates
meta <- read_csv(file.path(project_dir,
                           "results/integrated_sketch_rpca_md.csv"))

# Get sequences for clone of interest
source(file.path(project_dir, "scripts/clones_of_interest.R"))

# UMAP by sample (manually) ----

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
        geom_point(aes(colour = is_cl), size = 0.5) +
        theme_minimal(base_size = 15) +
        scale_color_manual(labels = c(FALSE, TRUE),
                           values = c("lightgray", default_pal[[cl]]),
                           guide = NULL) +
        labs(title = cl,
             x = NULL,
             y = NULL) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    pl
})

pdf(file.path(project_dir, "figures/umap_split_by_clusters.pdf"),
    width = 10, height = 15)
wrap_plots(cl_umaps, ncol = 3)
dev.off()


# UMAP with clone of interest highlighted, separated by sample ----

meta$is_coi = c(2, 0.5)[as.numeric(meta$coi == "no") + 1]
meta <- meta %>%
    dplyr::mutate(coi = factor(coi, levels = c("no", "Ri01", "Ri02"))) %>%
    dplyr::arrange(coi) %>%
    dplyr::mutate(coi = as.character(coi))

coi_colours <- structure(c("lightgray","#DC050C","steelblue"),
                         names = c("no", "Ri01", "Ri02"))


pdf(file.path(fig_dir, "umap_coi_by_sample_manual.pdf"),
    height = 12, width = 12)
ggplot(meta, 
       aes(x = UMAPfull_1, y = UMAPfull_2)) +
    geom_point(aes(size = is_coi, colour = coi)) +
    scale_size_identity() +
    facet_wrap(~Sample) +
    theme_minimal(base_size = 15) +
    scale_color_manual(labels = names(coi_colours), values = coi_colours) +
    labs(x = "UMAP dimension 1", y = "UMAP dimension 2") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

dev.off()
