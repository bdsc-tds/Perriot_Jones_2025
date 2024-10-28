# Libraries and setup ----

library("hdf5r")
library("rhdf5")
library("scRepertoire")
library("Seurat")
library("SeuratDisk")
library("tidyverse")

set.seed(7620)


# Filtering and integration parameters
min_cells <- 10 # keep features expressed in at least 10 cells
min_features <- 100 # require at least 100 genes expressed
max_percent_mt <- 10 # keep if not more than 10% mitochondrial reads 
n_sketch <- 5000 # Number of cells for sketch


# Check for directory environment variables 
project_dir <- tryCatch({Sys.getenv()[["project_dir"]]},
                        error = function(cond){return(".")})
scratch_dir <- tryCatch({Sys.getenv()[["scratch_dir"]]},
                        error = function(cond){return(".")})
data_dir <- file.path(project_dir, "data")
fig_dir <- file.path(project_dir, "figures")
log_dir <- file.path(project_dir, "logs")


# Select the samples to use
samples <- read_delim(file.path(data_dir, "cellbender_input.txt"),
                      col_names = FALSE)[[1]]

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
    
    # Add sample to cell names for easier integration of TCR
    colnames(decont_m) <- paste(nm, colnames(decont_m), sep = "_")
    
    samples <- data.frame(Sample = rep(nm, ncol(decont_m)),
                          row.names = colnames(decont_m))
    
    # Create Seurat objects, applying filtering for cells and features
    seurat_obj <- CreateSeuratObject(counts = decont_m[keep_rn, ], 
                                     min.cells = min_cells,
                                     min.features = min_features,
                                     meta.data = samples)
}    

# Filter for expression of CD8A or B or presence of TCR in metadata CTgene
filter_CD8 <- function(seurat_obj){
    cd8_expr <- subset(seurat_obj, features = c("CD8A", "CD8B"))
    cd8_expressed <- colSums(cd8_expr@assays$RNA$counts) > 0
    has_tcr <- ! is.na(seurat_obj@meta.data$CTgene)
    keep <- cd8_expressed | has_tcr
    seurat_obj <- subset(seurat_obj, cells = which(keep))
}


# Read TCR data ----

combined_tcr <- read_rds(file.path(data_dir, "processed/combined_tcr.rds"))

# Check that sample names are the same as in the config file
length(intersect(names(combined_tcr), samples)) == length(samples)

# Read cellbender filtered data ----

exp_dirs <- list.files(file.path(data_dir, "processed/cellbender/filtered"),
                       full.names = TRUE,
                       pattern = paste(samples, collapse = "|"))

names(exp_dirs) <- gsub("_cellbender.*", "", basename(exp_dirs))

all(samples %in% names(exp_dirs))
exp_dirs <- exp_dirs[samples]

# Get duplicated row names to drop
# (note that after cellbender, samples do not have the same features)
rn <- H5File$new(exp_dirs[[1]])[['matrix/features/name']][]
rhdf5::h5closeAll()
keep_rn <- ! ( duplicated(rn) | duplicated(rn, fromLast = TRUE) )

# Create Seurat objects from cellbender output ----

seurat_objs <- lapply(names(exp_dirs), function(nm){
    print(nm)
    obj <- h5_to_seurat(nm, exp_dirs[[nm]], min_cells, min_features)
    
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj,
                                                        pattern = "^MT-")
    
    # Filter high mitochondrial read percentage 
    obj <- subset(obj, subset = percent.mt < max_percent_mt)
    obj
    
})

# Prep for Seurat SketchData ----

# Merge and normalize
seurat_obj <- merge(seurat_objs[[1]], seurat_objs[2:length(seurat_objs)])
seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
seurat_obj <- NormalizeData(seurat_obj) # normalise for total count per cell

# Add TCR information 
seurat_obj <- combineExpression(combined_tcr, seurat_obj)

# Add beta chain cdr3
seurat_obj@meta.data <- seurat_obj@meta.data %>%
    tidyr::separate(CTgene,
                    into = c("TCR1", "TCR2"),
                    sep = "_", remove = FALSE) %>%
    tidyr::separate(CTaa,
                    into = c("CTaa1", "CTaa2"),
                    sep = "_", remove = FALSE) %>%
    dplyr::mutate(across(matches("TCR|CTaa"), ~na_if(.x, "NA")),
                  beta_aa = paste(TCR2, CTaa2, sep = "_")) 

sink(file.path(log_dir, "cd8_tcr_filtering.txt"))
print("before filtering")
table(seurat_obj$Sample)

print("after filtering")
seurat_obj <- filter_CD8(seurat_obj)
table(seurat_obj$Sample)
sink()

# Split and find variable features ---
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]],
                             f = seurat_obj$Sample)
seurat_obj <- FindVariableFeatures(seurat_obj)

# Sketch data ----
seurat_obj <- SketchData(object = seurat_obj,
                     ncells = n_sketch,
                     method = "LeverageScore",
                     sketched.assay = "sketch")
DefaultAssay(seurat_obj) <- "sketch"

# Find variable features and run PCA on sketched data ----
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# Integrate with RPCA ----
seurat_rpca <- IntegrateLayers(seurat_obj,
                          method = RPCAIntegration,
                          orig = "pca",
                          new.reduction = "integrated.rpca",
                          dims = 1:30,
                          k.anchor = 20)

seurat_rpca <- FindNeighbors(seurat_rpca,
                            reduction = "integrated.rpca",
                            dims = 1:30)
seurat_rpca <- FindClusters(seurat_rpca)
seurat_rpca <- RunUMAP(seurat_rpca,
                      reduction = "integrated.rpca",
                      dims = 1:30)

# UMAP on the integrated, sketched cells
pdf(file.path(fig_dir, "umap_sketch_rpca.pdf"))
DimPlot(seurat_rpca,
        group.by = "Sample",
        reduction = "umap")
dev.off()

#seurat_rpca[["sketch"]] <- JoinLayers(seurat_rpca[["sketch"]])

# At this point, only the 5000 sketched cells have cluster annotation 
# Project the other cells into the integrated space 

#seurat_rpca[["sketch"]] <- split(seurat_rpca[["sketch"]],
#                                f = seurat_rpca$Sample)

# This adds the dimension reduction integrated.rpca.full 
seurat_rpca <- ProjectIntegration(object = seurat_rpca,
                                 sketched.assay = "sketch",
                                 assay = "RNA",
                                 reduction = "integrated.rpca")

seurat_rpca <- ProjectData(object = seurat_rpca, 
                          sketched.assay = "sketch",
                          assay = "RNA",
                          sketched.reduction = "integrated.rpca",
                          full.reduction = "integrated.rpca.full",
                          refdata = "seurat_clusters",
                          #umap.model = "umap", to use the sketched umap?
                          dims = 1:30)

#DefaultAssay(obj) <- "RNA"? 
# Vignette https://satijalab.org/seurat/articles/seurat5_sketch_analysis
# extends UMAP, other vignette reruns it

seurat_rpca[["sketch"]] <- JoinLayers(seurat_rpca[["sketch"]])

# Run UMAP on the full projected data 
seurat_rpca <- RunUMAP(seurat_rpca,
                  reduction = "integrated.rpca.full",
                  dims = 1:30,
                  reduction.name = "umap.full",
                  reduction.key = "UMAPfull_")

# UMAP of integrated data, all samples
pdf(file.path(fig_dir, "umap_sketch_rpca_full_scratch.pdf"))
DimPlot(seurat_rpca,
        reduction = "umap.full", 
        group.by = "Sample")
dev.off()

# FeaturePlot with CD8A and CD8B 
pdf(file.path(fig_dir, "umap_sketch_rpca_cd8a_cd8b_scratch.pdf"),
    width = 10)
FeaturePlot(seurat_rpca,
            reduction = "umap.full", 
            features = c("CD8A", "CD8B"))
dev.off()


# Save integrated data ----
write_rds(seurat_rpca,
          "data/processed/integrated_sketch_rpca.rds")


seurat_rpca[[]] %>%
    dplyr::select(Sample, seurat_clusters) %>%
    dplyr::group_by(Sample, seurat_clusters) %>%
    dplyr::summarise(n_sample_cluster = n()) %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(n_sample = sum(n_sample_cluster),
                  pct_cluster = n_sample_cluster/n_sample * 100) %>%
    dplyr::arrange(desc(pct_cluster))
