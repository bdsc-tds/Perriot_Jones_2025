# ----------------------------------------------------------------------------
# Libraries and setup ----

library("Seurat")
library("SeuratDisk")
library("tidyverse")
library("argparse")

# Command line arguments ----
parser <- ArgumentParser(description = "Integrate seurat object")

parser$add_argument('--input', '-i',
                    help = 'Directory containing seurat object')
parser$add_argument('--output', '-o',
                    help = 'Directory to write integrated object')
parser$add_argument('--integration', 
                    help = 'Method to use for integration')
parser$add_argument('--figures',  '-f', 
                    help = 'Directory for storing UMAPs')
parser$add_argument('--n-sketch', '-n', dest = "n_sketch", 
                    help = 'Number of cells for sketch', default = 5000)
parser$add_argument('--cluster-res', '-c', dest = "cluster_res", 
                    help = 'Cluster resolution', default = 0.5)
args <- parser$parse_args()


# ----------------------------------------------------------------------------
# Functions ----
# make_umaps ----
# FIGURES NOT PRINTING
make_umaps <- function(seurat_obj, fig_dir){
    
    pdf(file.path(fig_dir, "umap_rpca.pdf"), width = 12, height = 12)
    p <- DimPlot(seurat_obj,
            reduction = "umap.full",
            group.by = "Sample")
    print(p)
    dev.off()
    
    pdf(file.path(fig_dir, "umap_rpca_by_sample.pdf"), width = 12, height = 12)
    p <- DimPlot(seurat_obj,
            reduction = "umap.full",
            group.by = "Sample",
            split.by = "Sample",
            ncol = 4)
    print(p)
    dev.off()
}
    
# run_sketch: Prep and sketch cells ----
run_sketch <- function(obj, n_sketch = 5000){
    # Split and find variable features ----
    obj[["RNA"]] <- split(obj[["RNA"]],
                          f = obj$Sample)
    
    # Splitting means that variable features are found for each sample
    obj <- FindVariableFeatures(obj)
    
    # Sketch data ----
    obj <- SketchData(object = obj,
                      ncells = n_sketch,
                      method = "LeverageScore",
                      sketched.assay = "sketch")
    DefaultAssay(obj) <- "sketch"
    
    # Find variable features and run PCA on sketched data ----
    obj <- FindVariableFeatures(obj)
    obj <- ScaleData(obj)
    obj <- RunPCA(obj)
    obj
}

# .run_integrate ----
.run_integrate <- function(obj,
                           integration_method,
                           nm,
                           cluster_res,
                           dims = 1:30,
                           k_anchor = 20){
    layer_nm <- sprintf("integrated.%s", nm)
    
    integrated_obj <- IntegrateLayers(obj,
                                      method = integration_method,
                                      orig = "pca",
                                      new.reduction = layer_nm,
                                      dims = dims,
                                      k.anchor = k_anchor)
    
    integrated_obj <- FindNeighbors(integrated_obj,
                                    reduction = layer_nm,
                                    dims = dims)
    
    integrated_obj <- FindClusters(integrated_obj, resolution = cluster_res)
    integrated_obj <- RunUMAP(integrated_obj,
                              reduction = layer_nm,
                              dims = 1:30)
    
    # Project to full data
    # At this point, only the sketched cells have cluster annotation 
    # Project the other cells into the integrated space 
    
    # This adds the dimension reduction integrated.rpca.full 
    integrated_obj <- ProjectIntegration(object = integrated_obj,
                                         sketched.assay = "sketch",
                                         assay = "RNA",
                                         reduction = layer_nm)
    
    full_layer_nm <- paste0(layer_nm, ".full")
    integrated_obj <- ProjectData(object = integrated_obj, 
                                  sketched.assay = "sketch",
                                  assay = "RNA",
                                  sketched.reduction = layer_nm,
                                  full.reduction = full_layer_nm,
                                  refdata = "seurat_clusters",
                                  dims = 1:30)
    
    integrated_obj[["sketch"]] <- JoinLayers(integrated_obj[["sketch"]])
    
    # Run UMAP on the full projected data 
    integrated_obj <- RunUMAP(integrated_obj,
                              reduction = full_layer_nm,
                              dims = 1:30,
                              reduction.name = "umap.full",
                              reduction.key = "UMAPfull_")
    integrated_obj
}

get_integration_method <- function(name){
    if (name == "rpca") return(RPCAIntegration)
    if (name == "harmony") return(HarmonyIntegration)
}


# Run integration on sketched cells, project to full ----
run_integration <- function(args){
    # Create directories if they don't exist
    fig_dir <- file.path(args$figures, args$integration)
    if (! file.exists(args$figures)){ dir.create(args$figures) }
    if (! file.exists(fig_dir)){ dir.create(fig_dir) }
    if (! file.exists(args$output)){ dir.create(args$output) }
    
    integration_f <- get_integration_method(args$integration)
    
    # Load seurat object
    seurat_obj <- read_rds(args$input)
    
    # Sketch cells
    seurat_obj <- run_sketch(seurat_obj, n_sketch = args$n_sketch)
    
    # Integrate and project to whole data set
    seurat_obj <- .run_integrate(seurat_obj,
                                 integration_f,
                                 args$integration,
                                 args$cluster_res)
    
    make_umaps(seurat_obj, fig_dir)
    
    # Add UMAP coordinates to metadata
    um <- Embeddings(seurat_obj[["umap.full"]])
    stopifnot(identical(rownames(um), rownames(seurat_obj[[]])))
    seurat_obj[[]] <- dplyr::bind_cols(seurat_obj[[]], um)
    
    # Save integrated seurat object 
    write_rds(seurat_obj, file.path(args$output, "integrated_seurat.rds"))
    
    # Save metadata 
    md <- as_tibble(seurat_obj[[]], rownames = "Cell")
    write_csv(md,
              file.path(args$output, "integrated_seurat.csv.gz"))
}
# ----------------------------------------------------------------------------
run_integration(args)