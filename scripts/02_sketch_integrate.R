# Libraries and setup ----

library("hdf5r")
library("rhdf5")
library("Seurat")
library("SeuratDisk")
library("tidyverse")

set.seed(7620)

min_cells <- 10
min_features <- 100
n_pcs <- 30
n_integration_features <- 3000

project_dir <- tryCatch({Sys.getenv()[["project_dir"]]},
                        error = function(cond){return(".")})
scratch_dir <- tryCatch({Sys.getenv()[["scratch_dir"]]},
                        error = function(cond){return(".")})
data_dir <- file.path(project_dir, "data")

# Functions ----

h5_to_seurat <- function(nm, exp_dir, min_cells, min_features){
    
    h5_data <- hdf5r::H5File$new(exp_dir, mode = 'r')
    
    decont_m <- Matrix::sparseMatrix(
        i = h5_data[['matrix/indices']][],
        p = h5_data[['matrix/indptr']][],
        x = h5_data[['matrix/data']][],
        dimnames = list(
            h5_data[['matrix/features/name']][],
            h5_data[['matrix/barcodes']][]
        ),
        dims = h5_data[['matrix/shape']][],
        index1 = FALSE
    )
    
    samples <- data.frame(Sample = rep(nm, ncol(decont_m)),
                          row.names = colnames(decont_m))
    
    seurat_obj <- CreateSeuratObject(counts = decont_m[keep_rn, ], 
                                     min.cells = min_cells,
                                     min.features = min_features,
                                     meta.data = samples)
}    

# Read data ----

# Assume the combined_tcr file only contains the samples of interest
combined_tcr <- read_rds(file.path(data_dir, "processed/combined_tcr.rds"))
samples <- names(combined_tcr)

# Read cellbender filtered data ----

exp_dirs <- list.files(file.path(data_dir, "processed/cellbender/filtered"),
                       full.names = TRUE,
                       pattern = paste(samples, collapse = "|"))

names(exp_dirs) <- gsub("_cellbender.*", "", basename(exp_dirs))

# Get duplicated row names to drop
# (note that after cellbender, samples do not have the same features)
rn <- H5File$new(exp_dirs[[1]])[['matrix/features/name']][]
rhdf5::h5closeAll()
keep_rn <- ! ( duplicated(rn) | duplicated(rn, fromLast = TRUE) )

# Create Seurat objects from cellbender output ----
seurat_objs <- lapply(names(exp_dirs), function(nm){
    obj <- h5_to_seurat(nm, exp_dirs[[nm]], min_cells, min_features)
    
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj,
                                                        pattern = "^MT-")
    obj
})

# Run Seurat SketchData ----

# Merge and normalize
seurat_obj <- merge(seurat_objs[[1]], seurat_objs[2:length(seurat_objs)])
seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
seurat_obj <- NormalizeData(seurat_obj)

# Split and find variable features
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]],
                             f = seurat_obj$Sample)
seurat_obj <- FindVariableFeatures(seurat_obj)

# Sketch data
seurat_obj <- SketchData(object = seurat_obj,
                     ncells = 5000,
                     method = "LeverageScore",
                     sketched.assay = "sketch")
DefaultAssay(seurat_obj) <- "sketch"

# PCA on sketched data
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# Integrated with RPCA
#seurat_obj <- IntegrateLayers(seurat_obj,
#                          method = RPCAIntegration,
#                          orig = "pca",
#                          new.reduction = "integrated.rpca",
#                          dims = 1:30,
#                          k.anchor = 20)

# Integrate with Harmony ----
seurat_obj <- IntegrateLayers(seurat_obj,
                              method = HarmonyIntegration,
                              orig = "pca",
                              new.reduction = "integrated.harmony")


seurat_obj <- FindNeighbors(seurat_obj,
                            reduction = "integrated.harmony",
                            dims = 1:30)
seurat_obj <- FindClusters(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj,
                      reduction = "integrated.harmony",
                      dims = 1:30)

seurat_obj[["sketch"]] <- JoinLayers(seurat_obj[["sketch"]])

pdf("umap_sketch_harmony.pdf")
DimPlot(seurat_obj, group.by = "Sample", reduction = "umap")
dev.off()


write_rds(seurat_obj, file = "integrated_sketch.rds")



seurat_obj[["sketch"]] <- split(seurat_obj[["sketch"]],
                                f = seurat_obj$Sample)

seurat_obj <- ProjectIntegration(object = seurat_obj,
                                 sketched.assay = "sketch",
                                 assay = "RNA",
                                 reduction = "integrated.harmony")


seurat_obj <- ProjectData(object = seurat_obj, 
                          sketched.assay = "sketch",
                          assay = "RNA",
                          sketched.reduction = "integrated.harmony",
                          full.reduction = "integrated.harmony.full",
                          dims = 1:30)


write_rds(seurat_obj, file = "integrated_sketch.rds")

