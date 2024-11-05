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
max_features <- 6000 # require at most 6000 genes expressed
max_percent_mt <- 10 # keep if not more than 10% mitochondrial reads 
n_sketch <- 5000 # Number of cells for sketch

# Check for directory environment variables 
project_dir <- tryCatch({Sys.getenv()[["project_dir"]]},
                        error = function(cond){return(".")})
scratch_dir <- tryCatch({Sys.getenv()[["scratch_dir"]]},
                        error = function(cond){return(".")})
data_dir <- file.path(project_dir, "data")
de_results <- file.path(project_dir, "results/differential_expression")
fig_dir <- file.path(project_dir, "figures")
log_dir <- file.path(project_dir, "logs")


# Select the samples to use
samples <- read_delim(file.path(data_dir, "cellbender_input.txt"),
                      col_names = FALSE)[[1]]

source(file.path(project_dir, "scripts/clones_of_interest.R"))

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

# Add TCR information ----
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

# Filtering ----

# Filter to require either CD8 expression or TCR sequence present
sink(file.path(log_dir, "cd8_tcr_filtering.txt"))
print("before CD8 filtering")
table(seurat_obj$Sample)

seurat_obj <- filter_CD8(seurat_obj)
print("after CD8 filtering")
table(seurat_obj$Sample)

# Filter for maximum genes expressed
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA < max_features)
print("after max gene filtering")
table(seurat_obj$Sample)
sink()

# Split and find variable features ----
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]],
                             f = seurat_obj$Sample)
# Splitting means that variable features are found for each sample
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

# Project to full data ----
# At this point, only the 5000 sketched cells have cluster annotation 
# Project the other cells into the integrated space 

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

# Note - vignette https://satijalab.org/seurat/articles/seurat5_sketch_analysis
# extends UMAP, other vignette reruns it

seurat_rpca[["sketch"]] <- JoinLayers(seurat_rpca[["sketch"]])

# Run UMAP on the full projected data ---- 

seurat_rpca <- RunUMAP(seurat_rpca,
                       reduction = "integrated.rpca.full",
                       dims = 1:30,
                       reduction.name = "umap.full",
                       reduction.key = "UMAPfull_")

# Save integrated data ----
write_rds(seurat_rpca,
          "integrated_sketch_rpca.rds")


# Add UMAP coordinates and clone of interest to metadata and save ----

um <- Embeddings(seurat_rpca[["umap.full"]])
identical(rownames(um), rownames(seurat_rpca[[]]))
seurat_rpca[[]] <- dplyr::bind_cols(seurat_rpca[[]], um)

seurat_rpca[[]] <- seurat_rpca[[]] %>%
    dplyr::mutate(coi = case_when(beta_aa == ri01_coi ~ "Ri01",
                                  beta_aa == ri02_coi ~ "Ri02",
                                  TRUE ~ "no"))

write_csv(seurat_rpca[[]],
          "integrated_sketch_rpca_md.csv")

### TO DO - write sketch umap coords

# ---------------
meta <- seurat_rpca[[]] 
# ---------------
# Differential expression ----

# Add condition * cluster to metadata
seurat_rpca$Condition <- as.factor(gsub("[0-9]+(-dis)?$", "", seurat_rpca$Sample))
seurat_rpca$Cond_Cl <- paste(seurat_rpca$Condition,
                             seurat_rpca$seurat_clusters,
                             sep = "_")

DefaultAssay(seurat_rpca) <- "RNA"
seurat_rpca[["RNA"]] <- JoinLayers(seurat_rpca[["RNA"]])

# Pseudobulk
pseudo_rpca <- AggregateExpression(seurat_rpca,
                                   assays = "RNA",
                                   return.seurat = TRUE,
                                   group.by = c("Sample", "seurat_clusters"))

# Differential expression between clusters ----
Idents(pseudo_rpca) <- "seurat_clusters"
bulk_cl_markers <- FindAllMarkers(object = pseudo_rpca,
                                  min.pct = 0,
                                  test.use = "DESeq2")

write_csv(bulk_cl_markers,
          file.path(de_results, "cluster_markers_all_samples.csv"))

features <- bulk_cl_markers %>%
    dplyr::filter(! grepl("^TR[AB]", gene),
                  p_val_adj <= 0.01,
                  abs(avg_log2FC) > 0.5) %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = 5) %>%
    pull(gene) %>%
    unique()

# Cluster heatmap for pseudobulked data
pdf(file.path(fig_dir, "sketch_rpca_cluster_heatmap_pseudobulk.pdf"),
    width = 12, height = 12)
DoHeatmap(pseudo_rpca,
          features = features) +
    theme(axis.text = element_text(size = 8))
dev.off()

#seurat_rpca <- ScaleData(seurat_rpca)

# Differential expression between clusters, without Ri01-5m ----

Idents(pseudo_rpca) <- pseudo_rpca$Sample
pseudo_subs <- subset(pseudo_rpca,
                    idents = setdiff(unique(pseudo_rpca$Sample), c("Ri01-5m")))


# Differential expression between HD and CL, per cluster ----

pseudo_subs$Condition <- as.factor(gsub("[0-9]+(-dis)?$", "", pseudo_subs$Sample))

Idents(pseudo_subs) <- pseudo_subs$Condition
pseudo_rpca_by_cl <- SplitObject(pseudo_subs,
                                 split.by = "seurat_clusters")


hd_vs_ri <- lapply(seq_along(pseudo_rpca_by_cl), function(i){
    mk <- FindAllMarkers(object = pseudo_rpca_by_cl[[i]],
                min.pct = 0,
                test.use = "DESeq2") 
    if (nrow(mk) > 0){
        mk <- mk %>%
            # Filter so logFC always refers to HD vs Ri
            dplyr::filter(cluster == "HD") %>%
            dplyr::mutate(cluster = i)
    }
    mk
})

hd_vs_ri <- bind_rows(hd_vs_ri)

write_csv(hd_vs_ri,
          file.path(de_results, "cluster_markers_hd_v_ri_without_Ri01_5m.csv"))

hd_vs_ri_sig <- hd_vs_ri %>%
    dplyr::filter(p_val_adj <= 0.01,
                  abs(avg_log2FC) > 0.5)

write_csv(hd_vs_ri_sig,
          file.path(de_results,
                    "cluster_markers_hd_v_ri_without_Ri01_5m_sig.csv"))



# Run scRepertoire positional entropy

clusters_w_roi = seurat_rpca[[]] %>%
    dplyr::filter(! coi == "no") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n >= 5) %>%
    dplyr::pull(seurat_clusters) %>%
    unique()

coi_cl_subs <- subset(seurat_rpca,
                      subset = seurat_clusters %in% clusters_w_roi)
coi_cl_subs <- SplitObject(coi_cl_subs, split.by = "seurat_clusters")

dummy <- lapply(names(coi_cl_subs), function(nm){
    pdf(file.path(fig_dir, 
                  sprintf("positional_entropy_cl_%s.pdf", nm)))
    scRepertoire::positionalEntropy(coi_cl_subs)
    dev.off()
})

# Error
#  Please provide rownames to the matrix before converting to a Graph
# WORKING HERE 
    

# 




# TO DO - adjust TCR for cellbender ----

