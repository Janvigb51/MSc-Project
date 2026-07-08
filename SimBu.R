# Load Packages & Libraries
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("SimBu")
BiocManager::install("SummarizedExperiment")
BiocManager::install("BiocParallel")

install.packages("Seurat")
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(Matrix)
  library(harmony)
  library(SimBu)
  library(tibble)
  library(tidyr)
  library(SummarizedExperiment)})

## Introduction
### The purpose of this Rmd is to generate 100 pseudobulk datasets of CAF subpopulations using the `SimBu` R pacakge.
### Thanks to John O'Grady for suggesting this R package.
breast_data <- readRDS("C:/Users/janvi/Desktop/MSc Project/scRNA-seq_dataobjects/scRNA-seq/BREAST_fibro_tumour.rds")
lung_data <- readRDS("C:/Users/janvi/Desktop/MSc Project/scRNA-seq_dataobjects/scRNA-seq/LUNG_fibro_tumour.rds")
# don't read in with header as R does not like duplicate colnames
# load in metadata
breast_meta <- breast_data@meta.data #output: 18984 genes, 16704 cells
lung_meta <- lung_data@meta.data #output: 25765 genes, 1377 cells

### 1) Create the Sparse Count Matrix (genes x cells)
sparse_matrix_breast <- breast_data@assays$RNA@counts
sparse_matrix_lung <- lung_data@assays$RNA@counts

# Make sure metadata order matches count matrix columns
breast_meta <- breast_meta[colnames(sparse_matrix_breast), ]
lung_meta <- lung_meta[colnames(sparse_matrix_lung), ]
# check match
stopifnot(all(rownames(breast_meta) == colnames(sparse_matrix_breast)))
stopifnot(all(rownames(lung_meta) == colnames(sparse_matrix_lung)))

# Create new simple cell IDs
new_breast_cell_ids <- paste0("cell", seq_len(ncol(sparse_matrix_breast)))
new_lung_cell_ids <- paste0("cell", seq_len(ncol(sparse_matrix_lung)))
  
# Save old-to-new mapping in case we need it later
breast_cell_id_map <- data.frame(old_id = colnames(sparse_matrix_breast), new_id = new_breast_cell_ids) 
lung_cell_id_map <- data.frame(old_id = colnames(sparse_matrix_lung), new_id = new_lung_cell_ids)

# Rename columns of count matrix
colnames(sparse_matrix_breast) <- new_breast_cell_ids
colnames(sparse_matrix_lung) <- new_lung_cell_ids

# Rename rows of metadata
rownames(breast_meta) <- new_breast_cell_ids
rownames(lung_meta) <- new_lung_cell_ids

### 2) Create the Annotation Table (cells x info)
breast_annotation <- data.frame(ID = rownames(breast_meta), cell_type = breast_meta$CAFtype)
view(breast_annotation)
lung_annotation <- data.frame(ID = rownames(lung_meta), cell_type = lung_meta$cluster_ft)
view(lung_annotation)
# check match
stopifnot(all(breast_annotation$ID == colnames(sparse_matrix_breast)))
stopifnot(all(lung_annotation$ID == colnames(sparse_matrix_lung)))

# Save cleaned SimBu input objects
simbu_input_dir <- "simbu_inputs"
dir.create(simbu_input_dir, showWarnings = FALSE)
saveRDS(sparse_matrix_breast,
        file.path(simbu_input_dir, "breast_sparse_matrix.rds"))
saveRDS(sparse_matrix_lung,
        file.path(simbu_input_dir, "lung_sparse_matrix.rds"))
saveRDS(breast_annotation,
        file.path(simbu_input_dir, "breast_annotation.rds"))
saveRDS(lung_annotation,
        file.path(simbu_input_dir, "lung_annotation.rds"))
saveRDS(breast_cell_id_map,
        file.path(simbu_input_dir, "breast_cell_id_map.rds"))
saveRDS(lung_cell_id_map,
        file.path(simbu_input_dir, "lung_cell_id_map.rds"))

### 3) Simulation
cords_breast <- SimBu::dataset(
  annotation = breast_annotation, #cell_id belonging to a cell type
  count_matrix = sparse_matrix_breast, #genes per unique cell_ids
  #tpm_matrix = NULL,
  name = "CAFtype")

cords_lung <- SimBu::dataset(
  annotation = lung_annotation, #cell_id belonging to a cell type
  count_matrix = sparse_matrix_lung, #genes per unique cell_ids
  #tpm_matrix = NULL,
  name = "CAFtype")

## SAVE SIMBU DATASET OBEJCTS
simbu_input_dir <- "simbu_inputs"
dir.create(simbu_input_dir, showWarnings = FALSE)
saveRDS(cords_breast,
        file.path(simbu_input_dir, "cords_breast_simbu.rds"))
saveRDS(cords_lung,
        file.path(simbu_input_dir, "cords_lung_simbu.rds"))

## RUN SIMULATION
print("simulating data")
breast_simulation <- SimBu::simulate_bulk(
  data = cords_breast,
  scenario = "mirror_db",
  ncells = 2000,
  nsamples = 100,
  scaling_factor = "read_number", # use number of reads as scaling factor
  BPPARAM = BiocParallel::MulticoreParam(workers = 3), # this will use 4 threads to run the simulation
  run_parallel = TRUE,
  balance_even_mirror_scenario = 0.01,
  seed = 20240618)

print("simulating data")
lung_simulation <- SimBu::simulate_bulk(
  data = cords_lung,
  scenario = "mirror_db",
  ncells = 2000,
  nsamples = 100,
  scaling_factor = "read_number", # use number of reads as scaling factor
  BPPARAM = BiocParallel::MulticoreParam(workers = 3), # this will use 4 threads to run the simulation
  run_parallel = TRUE,
  balance_even_mirror_scenario = 0.01,
  seed = 20240618)

bulk_counts_breast <- as.data.frame(SummarizedExperiment::assays(breast_simulation$bulk)[["bulk_counts"]]) %>% rownames_to_column(var = "Gene")
print("writing to file")
write.table(bulk_counts_breast, file = "C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-breast-100samples-2000cells-20240618/bulk_counts_2000cells_100samps.txt", quote = F, sep = "\t", row.names = F)
write.table(breast_simulation$cell_fractions, file = "C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-breast-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt", quote = F, sep = "\t", row.names = T)
# breast_simulation$bulk is our pseudo bulk samples that we feed in Insta Prism
# breast_simulation$cell_fractions is our ground truth data that we use to compare deconv results (the estimated)

bulk_counts_lung <- as.data.frame(SummarizedExperiment::assays(lung_simulation$bulk)[["bulk_counts"]]) %>% rownames_to_column(var = "Gene")
print("writing to file")
write.table(bulk_counts_lung, file = "C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-lung-100samples-2000cells-20240618/bulk_counts_2000cells_100samps.txt", quote = F, sep = "\t", row.names = F)
write.table(lung_simulation$cell_fractions, file = "C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-lung-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt", quote = F, sep = "\t", row.names = T)

