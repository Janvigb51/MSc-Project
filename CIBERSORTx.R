# Load any required packages
library(omnideconv)
library(SimBu)
library(Matrix)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(harmony)
library(tibble)
library(tidyr)
library(ggplot2)

# Load the datasets to be analysed
breast_data <- readRDS("../scRNA-seq_dataobjects/scRNA-seq/BREAST_fibro_tumour.rds")
lung_data <- readRDS("../scRNA-seq_dataobjects/scRNA-seq/LUNG_fibro_tumour.rds")
# load in metadata
breast_meta <- breast_data@meta.data #output: 18984 genes, 16704 cells
lung_meta <- lung_data@meta.data #output: 25765 genes, 1377 cells

### FOR OMNIDECONV CIBERSORTx
# set_cibersortx_credentials(email,token) # replaced with personal credentials

### 1) Create single_cell_object (genes x cells)
breast_sc_counts <- breast_data@assays$RNA@counts
lung_sc_counts <- lung_data@assays$RNA@counts

### 2) Create cell_type_annotations (CAF subtypes)
breast_cell_types <- as.character(breast_meta$CAFtype)
lung_cell_types <- as.character(lung_meta$cluster_ft)
# Clear any spaces if needed
breast_cell_types <- gsub(" ", "_", breast_cell_types)
lung_cell_types <- gsub(" ", "_", lung_cell_types)
# Check each dataset's CAF types
table(breast_cell_types)
table(lung_cell_types)
# Labels and cells should be in the same order
stopifnot(all(rownames(breast_meta) == colnames(breast_sc_counts)))
stopifnot(all(rownames(lung_meta) == colnames(lung_sc_counts)))

### 3) Create 3rd object: bulk_tpm  (genes x bulk samples, TPM/CPM-normalised)
breast_bulk <- read.delim("../deconv_inputs/cords-breast-100samples-2000cells-20240618/bulk_counts_2000cells_100samps.txt",
header = TRUE, sep = "\t", check.names = FALSE)
lung_bulk <- read.delim("../deconv_inputs/cords-lung-100samples-2000cells-20240618/bulk_counts_2000cells_100samps.txt",
header = TRUE, sep = "\t", check.names = FALSE)
# Convert the gene columns into row names
rownames(breast_bulk) <- breast_bulk$Gene
breast_bulk$Gene <- NULL
rownames(lung_bulk) <- lung_bulk$Gene
lung_bulk$Gene <- NULL
# Convert bulk data to a numeric matrix
breast_bulk <- as.matrix(breast_bulk)
mode(breast_bulk) <- "numeric"
lung_bulk <- as.matrix(lung_bulk)
mode(lung_bulk) <- "numeric"

breast_bulk_tpm <- sweep(breast_bulk, 2, colSums(breast_bulk), "/") * 1e6
lung_bulk_tpm <- sweep(lung_bulk, 2, colSums(lung_bulk), "/") * 1e6

### 4) True CAF proportions from SimBu (for performance evaluation)
breast_truth <- read.delim("../deconv_inputs/cords-breast-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt",
header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
lung_truth <- read.delim("../deconv_inputs/cords-lung-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt",
header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

### 5) Check and Match genes between single cell and bulk
common_genes_breast <- intersect(rownames(breast_sc_counts), rownames(breast_bulk_tpm))
breast_sc_counts2 <- breast_sc_counts[common_genes_breast, , drop = FALSE]
breast_bulk_tpm2  <- breast_bulk_tpm[common_genes_breast, , drop = FALSE]
all(rownames(breast_sc_counts2) == rownames(breast_bulk_tpm2))
common_genes_lung <- intersect(rownames(lung_sc_counts), rownames(lung_bulk_tpm))
lung_sc_counts2 <- lung_sc_counts[common_genes_lung, , drop = FALSE]
lung_bulk_tpm2 <- lung_bulk_tpm[common_genes_lung, , drop = FALSE]
all(rownames(lung_sc_counts2) == rownames(lung_bulk_tpm2))

### FOR OMNIDECONV METHODS:
### Build Signature Matrix for CIBERSORTx

dir.create("cibersortx_inputs/breast", recursive = TRUE, showWarnings = FALSE)
dir.create("cibersortx_results/breast", recursive = TRUE, showWarnings = FALSE)

# Make sure CAF annotations are plain character labels
breast_cell_types <- as.character(breast_cell_types)
# Make sure annotations match the columns/cells of the scRNA-seq matrix
names(breast_cell_types) <- colnames(breast_sc_counts2)
stopifnot(ncol(breast_sc_counts2) == length(breast_cell_types))
# Convert sparse matrix to normal matrix for CIBERSORTx
breast_sc_counts2_mat <- as.matrix(breast_sc_counts2)
storage.mode(breast_sc_counts2_mat) <- "numeric"

breast_signature_cibersort <- omnideconv::build_model_cibersortx(
  single_cell_object = breast_sc_counts2_mat,
  cell_type_annotations = breast_cell_types,
  container = "apptainer",
  container_path = "/data/containers/",
  verbose = TRUE,
  input_dir = "cibersortx_inputs/breast",
  output_dir = "cibersortx_results/breast",
  display_heatmap = TRUE,
  k_max = 999,
  filter = FALSE,
  sampling = 1)
saveRDS(breast_signature_cibersort, "cibersortx_results/breast_cibersort_signature.rds")

lung_signature_cibersort <- omnideconv::build_model_cibersortx(
  single_cell_object = lung_sc_counts2,
  cell_type_annotations = lung_cell_types,
  container = c("docker", "apptainer"),
  container_path = "/data/containers/",
  verbose = TRUE,
  input_dir = "cibersortx_inputs/lung",
  output_dir = "cibersortx_results",
  display_heatmap = TRUE,
  k_max = 999,
  filter = FALSE,
  sampling = 1)
saveRDS(lung_signature_cibersort, "cibersortx_results/lung_cibersort_signature.rds")

