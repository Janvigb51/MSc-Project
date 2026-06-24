# Load the required packages
library(omnideconv)
library(SimBu)
library(Matrix)
library(Seurat)

# Load the datasets to be analysed
breast_data <- readRDS("../scRNA-seq_dataobjects/scRNA-seq/BREAST_fibro_tumour.rds")
lung_data <- readRDS("../scRNA-seq_dataobjects/scRNA-seq/LUNG_fibro_tumour.rds")
# load in metadata
breast_meta <- breast_data@meta.data #output: 18984 genes, 16704 cells
lung_meta <- lung_data@meta.data #output: 25765 genes, 1377 cells

### FOR OMNIDECONV
### 1) Create 1st object raw sc_counts (genes x cells)
breast_sc_counts <- breast_data@assays$RNA@counts
lung_sc_counts <- lung_data@assays$RNA@counts

### 2) Create 2nd object cell_types (CAF subtypes)
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
