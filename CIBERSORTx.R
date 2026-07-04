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

### 1) Create single_cell count matrices (genes x cells)
breast_sc_counts <- breast_data@assays$RNA@counts
lung_sc_counts <- lung_data@assays$RNA@counts

### 2) Load in metadata and create CAF subtype annotations
breast_meta <- breast_data@meta.data #output: 18984 genes, 16704 cells
lung_meta <- lung_data@meta.data #output: 25765 genes, 1377 cells
# Check metadata and count matrix cell order match
stopifnot(all(rownames(breast_meta) == colnames(breast_sc_counts)))
stopifnot(all(rownames(lung_meta) == colnames(lung_sc_counts)))
# Extract CAF subtypes
breast_cell_types <- breast_meta$CAFtype
lung_cell_types <- lung_meta$cluster_ft
# Check each dataset's CAF types
table(breast_cell_types)
table(lung_cell_types)

### 3) Load simulated bulk RNA-seq count matrices
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
# Convert bulk counts to CPM/TPM-like normalised values
breast_bulk_tpm <- sweep(breast_bulk, 2, colSums(breast_bulk), "/") * 1e6
lung_bulk_tpm <- sweep(lung_bulk, 2, colSums(lung_bulk), "/") * 1e6

### 4) True CAF proportions from SimBu (for performance evaluation)
breast_truth <- read.delim("../deconv_inputs/cords-breast-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt",
header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
lung_truth <- read.delim("../deconv_inputs/cords-lung-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt",
header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

### 5) Check and Match genes between single cell and bulk data
common_genes_breast <- intersect(rownames(breast_sc_counts), rownames(breast_bulk_tpm))
breast_sc_counts2 <- breast_sc_counts[common_genes_breast, , drop = FALSE]
breast_bulk_tpm2  <- breast_bulk_tpm[common_genes_breast, , drop = FALSE]
stopifnot(all(rownames(breast_sc_counts2) == rownames(breast_bulk_tpm2)))

common_genes_lung <- intersect(rownames(lung_sc_counts), rownames(lung_bulk_tpm))
lung_sc_counts2 <- lung_sc_counts[common_genes_lung, , drop = FALSE]
lung_bulk_tpm2 <- lung_bulk_tpm[common_genes_lung, , drop = FALSE]
stopifnot(all(rownames(lung_sc_counts2) == rownames(lung_bulk_tpm2)))

### FOR CIBERSORTx SIGNATURE MATRIX CREATION
### CIBERSORTx was run manually on Lugh using the CIBERSORTx fractions container.
### The omnideconv::build_model_cibersortx() function was not used because it failed
### with an rbind.Matrix type error. Therefore, the CIBERSORTx input file was generated
### and then the signature matrix was built using a SLURM job and Singularity container.

dir.create("cibersortx/cibersortx_inputs/breast", recursive = TRUE, showWarnings = FALSE)
dir.create("cibersortx/cibersortx_results/breast", recursive = TRUE, showWarnings = FALSE)
dir.create("cibersortx/cibersortx_inputs/lung", recursive = TRUE, showWarnings = FALSE)
dir.create("cibersortx/cibersortx_results/lung", recursive = TRUE, showWarnings = FALSE)

# Function to create CIBERSORTx single-cell reference input file
write_cibersortx_input <- function(sc_counts, cell_types, output_file) {
  # CAF annotations should be plain character labels
  cell_types <- as.character(cell_types)
  cell_types <- gsub(" ", "_", cell_types)
  # Annotations should match the columns/cells of the scRNA-seq matrix
  names(cell_types) <- colnames(sc_counts)
  stopifnot(ncol(sc_counts) == length(cell_types))
  stopifnot(all(names(cell_types) == colnames(sc_counts)))
  # Convert sparse matrix to normal numeric matrix for CIBERSORTx
  sc_counts_mat <- as.matrix(sc_counts)
  storage.mode(sc_counts_mat) <- "numeric"
  # First row: GeneSymbol + CAF subtype label for each cell
  first_row <- c("GeneSymbol", cell_types)
  # Expression table: first column = gene names, remaining columns = cells
  expr_table <- cbind(
    GeneSymbol = rownames(sc_counts_mat),
    sc_counts_mat)
  # Combine phenotype row + expression matrix
  cibersortx_ref <- rbind(first_row, expr_table)
  write.table(
    cibersortx_ref,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE)}

### Create CIBERSORTx input file for breast
# this creates: cibersortx/cibersortx_inputs/breast/sample_file_for_cibersort.txt
breast_cibersortx_input_file <- "cibersortx/cibersortx_inputs/breast/sample_file_for_cibersort.txt"
if (!file.exists(breast_cibersortx_input_file)) {
  write_cibersortx_input(
    sc_counts = breast_sc_counts2,
    cell_types = breast_cell_types,
    output_file = breast_cibersortx_input_file)}

### Create CIBERSORTx input file for lung
# this creates: cibersortx/cibersortx_inputs/lung/sample_file_for_cibersort.txt
lung_cibersortx_input_file <- "cibersortx/cibersortx_inputs/lung/sample_file_for_cibersort.txt"
if (!file.exists(lung_cibersortx_input_file)) {
  write_cibersortx_input(
    sc_counts = lung_sc_counts2,
    cell_types = lung_cell_types,
    output_file = lung_cibersortx_input_file)}

### The CIBERSORTx signature matrices are generated outside R using:
### build_signature_matrix_breast.job / build_signature_matrix_lung.job
### on Lugh with the CIBERSORTx Fractions Singularity container.

### Breast CIBERSORTx signature matrix output
cibersortx_breast_signature_file <- file.path("cibersortx","cibersortx_results","breast",
"CIBERSORTx_sample_file_for_cibersort_inferred_phenoclasses.CIBERSORTx_sample_file_for_cibersort_inferred_refsample.bm.K999.txt")

if (file.exists(cibersortx_breast_signature_file)) {
  breast_signature <- read.delim(
    cibersortx_breast_signature_file,
    header = TRUE,sep = "\t",check.names = FALSE)
  
  rownames(breast_signature) <- breast_signature$NAME
  breast_signature$NAME <- NULL
  dim(breast_signature)
  head(breast_signature[, 1:5])} else {
  message("Breast CIBERSORTx signature matrix not found yet. Run the breast SLURM job first.")}

### Lung CIBERSORTx signature matrix output
cibersortx_lung_signature_file <- file.path("cibersortx","cibersortx_results","lung",
"CIBERSORTx_sample_file_for_cibersort_inferred_phenoclasses.CIBERSORTx_sample_file_for_cibersort_inferred_refsample.bm.K999.txt")

if (file.exists(cibersortx_lung_signature_file)) {
  lung_signature <- read.delim(
    cibersortx_lung_signature_file,
    header = TRUE,sep = "\t",check.names = FALSE)
  
  rownames(lung_signature) <- lung_signature$NAME
  lung_signature$NAME <- NULL
  dim(lung_signature)
  head(lung_signature[, 1:5])} else {
  message("Lung CIBERSORTx signature matrix not found yet. Run the lung SLURM job first.")}

