# If deconvolution has already been run
# jump to QUICK START (line 201)

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

# SAVE RESULTS
results_dir <- "C:/Users/janvi/Desktop/MSc Project/MSc Project R/cibersortx/cibersortx_results/"
saveRDS(breast_truth, file.path(results_dir, "cibersort_breast_truth.rds"))
saveRDS(lung_truth, file.path(results_dir, "cibersort_lung_truth.rds"))

### 5) Check and Match genes between single cell and bulk data
common_genes_breast <- intersect(rownames(breast_sc_counts), rownames(breast_bulk_tpm))
breast_sc_counts2 <- breast_sc_counts[common_genes_breast, , drop = FALSE]
breast_bulk_tpm2  <- breast_bulk_tpm[common_genes_breast, , drop = FALSE]
stopifnot(all(rownames(breast_sc_counts2) == rownames(breast_bulk_tpm2)))

common_genes_lung <- intersect(rownames(lung_sc_counts), rownames(lung_bulk_tpm))
lung_sc_counts2 <- lung_sc_counts[common_genes_lung, , drop = FALSE]
lung_bulk_tpm2 <- lung_bulk_tpm[common_genes_lung, , drop = FALSE]
stopifnot(all(rownames(lung_sc_counts2) == rownames(lung_bulk_tpm2)))

### 6) FOR CIBERSORTx SIGNATURE MATRIX CREATION
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

# Create CIBERSORTx input file for breast
# this creates: cibersortx/cibersortx_inputs/breast/sample_file_for_cibersort.txt
breast_cibersortx_input_file <- "cibersortx/cibersortx_inputs/breast/sample_file_for_cibersort.txt"
if (!file.exists(breast_cibersortx_input_file)) {
  write_cibersortx_input(
    sc_counts = breast_sc_counts2,
    cell_types = breast_cell_types,
    output_file = breast_cibersortx_input_file)}

# Create CIBERSORTx input file for lung
# this creates: cibersortx/cibersortx_inputs/lung/sample_file_for_cibersort.txt
lung_cibersortx_input_file <- "cibersortx/cibersortx_inputs/lung/sample_file_for_cibersort.txt"
if (!file.exists(lung_cibersortx_input_file)) {
  write_cibersortx_input(
    sc_counts = lung_sc_counts2,
    cell_types = lung_cell_types,
    output_file = lung_cibersortx_input_file)}

### The CIBERSORTx signature matrices are generated outside R using:
### build_signature_matrix_template.job adapted to breast & lung data
### on Lugh with the CIBERSORTx Fractions Singularity container.

# Breast CIBERSORTx signature matrix output
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

# Lung CIBERSORTx signature matrix output
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

### 7) Function to create CIBERSORTx Mixture File
write_cibersortx_mixture <- function(bulk_matrix, signature_matrix, output_file) {
  common_genes <- intersect(rownames(signature_matrix), rownames(bulk_matrix))
  bulk_matrix2 <- bulk_matrix[common_genes, , drop = FALSE]
  bulk_matrix2 <- as.matrix(bulk_matrix2)
  storage.mode(bulk_matrix2) <- "numeric"
  mixture_file <- cbind(
    GeneSymbol = rownames(bulk_matrix2),
    bulk_matrix2)
  write.table(
    mixture_file,
    file = output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE)}

# Create CIBERSORTx mixture file for breast
write_cibersortx_mixture(
  bulk_matrix = breast_bulk_tpm2,
  signature_matrix = breast_signature,
  output_file = "cibersortx/cibersortx_inputs/breast/mixture_file_for_cibersort.txt")

# Create CIBERSORTx mixture file for lung
write_cibersortx_mixture(
  bulk_matrix = lung_bulk_tpm2,
  signature_matrix = lung_signature,
  output_file = "cibersortx/cibersortx_inputs/lung/mixture_file_for_cibersort.txt")

### CIBERSORTx deconvolution is performed outside R (on Lugh) using:
### run_cibersortx_deconvolution_template.job adapted to breast & lung data

############################################################
### QUICK START: Resume from saved CIBERSORTx outputs
### Use this section after CIBERSORTx has already been run on Lugh.
### This reloads the saved truth files and reads the completed
### CIBERSORTx_Adjusted.txt deconvolution outputs for evaluation/plotting.
############################################################
results_dir <- "C:/Users/janvi/Desktop/MSc Project/MSc Project R/cibersortx/cibersortx_results"
breast_truth <- readRDS(file.path(results_dir, "cibersort_breast_truth.rds"))
lung_truth <- readRDS(file.path(results_dir, "cibersort_lung_truth.rds"))

### 8) READ IN CIBERSORTx DECONVOLUTION RESULTS
breast_cibersortx_result_file <- file.path(results_dir, "breast_deconvolution", "CIBERSORTx_Adjusted.txt")
lung_cibersortx_result_file <- file.path (results_dir, "lung_deconvolution", "CIBERSORTx_Adjusted.txt")

breast_cibersortx_raw <- read.delim(
  breast_cibersortx_result_file,
  header = TRUE,sep = "\t",row.names = 1,check.names = FALSE)

lung_cibersortx_raw <- read.delim(
  lung_cibersortx_result_file,
  header = TRUE,sep = "\t",row.names = 1,check.names = FALSE)

head(breast_cibersortx_raw)
head(lung_cibersortx_raw)
colnames(breast_cibersortx_raw)
colnames(lung_cibersortx_raw)

# Clean CIBERSORTx estimated matrices (matching CAF sybtypes in truth & estimated)
breast_cibersortx_est <- breast_cibersortx_raw[
  , intersect(colnames(breast_truth), colnames(breast_cibersortx_raw)),
  drop = FALSE]
lung_cibersortx_est <- lung_cibersortx_raw[
  , intersect(colnames(lung_truth), colnames(lung_cibersortx_raw)),
  drop = FALSE]
# Convert to numeric matrices
breast_cibersortx_est <- as.matrix(breast_cibersortx_est)
storage.mode(breast_cibersortx_est) <- "numeric"
lung_cibersortx_est <- as.matrix(lung_cibersortx_est)
storage.mode(lung_cibersortx_est) <- "numeric"

### 9) Performance Evaluation Function
evaluate_deconv <- function(truth, estimated) {
  # Find shared samples and CAF types
  common_samples <- intersect(rownames(truth), rownames(estimated))
  common_celltypes <- intersect(colnames(truth), colnames(estimated))
  # Match truth and estimated matrices
  truth <- truth[common_samples, common_celltypes, drop = FALSE]
  estimated <- estimated[common_samples, common_celltypes, drop = FALSE]
  # Convert to numeric matrices
  truth <- as.matrix(truth)
  estimated <- as.matrix(estimated)
  storage.mode(truth) <- "numeric"
  storage.mode(estimated) <- "numeric"
  # Safety checks
  stopifnot(all(rownames(truth) == rownames(estimated)))
  stopifnot(all(colnames(truth) == colnames(estimated)))
  # Global Pearson correlation
  global_correlation <- cor(
    as.vector(truth),
    as.vector(estimated),
    use = "complete.obs")
  # Cell-type-specific Pearson correlation
  correlation_by_celltype <- sapply(colnames(truth), function(ct) {
    cor(truth[, ct], estimated[, ct], use = "complete.obs")})
  # Global RMSE
  global_rmse <- sqrt(mean((truth - estimated)^2, na.rm = TRUE))
  # Cell-type-specific RMSE
  rmse_by_celltype <- sapply(colnames(truth), function(ct) {
    sqrt(mean((truth[, ct] - estimated[, ct])^2, na.rm = TRUE))})
  # Return everything together
  return(list(
    truth_matched = truth,
    estimated_matched = estimated,
    global_correlation = global_correlation,
    correlation_by_celltype = correlation_by_celltype,
    global_rmse = global_rmse,
    rmse_by_celltype = rmse_by_celltype))}

# Evaluate Breast CIBERSORTx
breast_cibs_eval <- evaluate_deconv(
  truth = breast_truth,
  estimated = breast_cibersortx_est)
breast_cibs_eval$global_correlation
breast_cibs_eval$correlation_by_celltype
breast_cibs_eval$global_rmse
breast_cibs_eval$rmse_by_celltype

# Evaluate Lung CIBERSORTx
lung_cibs_eval <- evaluate_deconv(
  truth = lung_truth,
  estimated = lung_cibersortx_est)
lung_cibs_eval$global_correlation
lung_cibs_eval$correlation_by_celltype
lung_cibs_eval$global_rmse
lung_cibs_eval$rmse_by_celltype

# Save global performance results
cibersortx_global_metrics <- data.frame(
  Dataset = c("Breast", "Lung"),
  Method = "CIBERSORTx",
  Global_Correlation = c(
    breast_cibs_eval$global_correlation,
    lung_cibs_eval$global_correlation),
  Global_RMSE = c(
    breast_cibs_eval$global_rmse,
    lung_cibs_eval$global_rmse))

cibersortx_global_metrics
write.csv(cibersortx_global_metrics,
file.path(results_dir, "cibersortx_global_metrics.csv"),
row.names = FALSE)

# Save per-cell-type results
cibersortx_celltype_metrics <- rbind(
  data.frame(
    Dataset = "Breast",
    Method = "CIBERSORTx",
    CAFtype = names(breast_cibs_eval$correlation_by_celltype),
    Correlation = as.numeric(breast_cibs_eval$correlation_by_celltype),
    RMSE = as.numeric(breast_cibs_eval$rmse_by_celltype)),
  data.frame(
    Dataset = "Lung",
    Method = "CIBERSORTx",
    CAFtype = names(lung_cibs_eval$correlation_by_celltype),
    Correlation = as.numeric(lung_cibs_eval$correlation_by_celltype),
    RMSE = as.numeric(lung_cibs_eval$rmse_by_celltype)))

cibersortx_celltype_metrics
write.csv(cibersortx_celltype_metrics,
file.path(results_dir, "cibersortx_celltype_metrics.csv"),
row.names = FALSE)

### 10) Plot Results colored by CAF type
### BREAST
# Take the matching matrices
breast_truth_plot <- breast_cibs_eval$truth_matched
breast_est_plot <- breast_cibs_eval$estimated_matched

plot_breast_cibs <- data.frame(
  sample = rep(rownames(breast_truth_plot), times = ncol(breast_truth_plot)), 
  CAFtype = rep(colnames(breast_truth_plot), each = nrow(breast_truth_plot)), 
  truth = as.vector(as.matrix(breast_truth_plot)), 
  estimated = as.vector(as.matrix(breast_est_plot))) 
caf_colors_breast <- c(
  "apCAF" = "black",
  "dCAF" = "palevioletred2",
  "hsp_tpCAF" = "limegreen",
  "iCAF" = "dodgerblue3",
  "IDO_CAF" = "darkorange",
  "mCAF" = "mediumorchid3",
  "Pericyte" = "goldenrod2",
  "rCAF" = "grey60",
  "tpCAF" = "turquoise3",
  "vCAF" = "hotpink3")
plot(plot_breast_cibs$truth, 
     plot_breast_cibs$estimated, 
     xlab = "True cell-type proportions", 
     ylab = "Estimated cell-type proportions", 
     main = "CIBERSORTx deconvolution performance (breast data)", 
     pch = 16, 
     col = caf_colors_breast[plot_breast_cibs$CAFtype]) 

abline(0, 1, col = "red", lty = 2, lwd = 1.8)

### Add global correlation and RMSE to plot
metrics_text_breast <- paste0(
  "Global Pearson Correlation (r) = ", round(breast_cibs_eval$global_correlation, 3),
  "\nGlobal RMSE = ", round(breast_cibs_eval$global_rmse, 3))
text(
  x = 0.04,
  y = 0.45,
  labels = metrics_text_breast,
  adj = c(0, 1),
  cex = 0.95,
  font = 2)

legend(
  x = - 0.005,
  y = 0.47,
  legend = levels(as.factor(plot_breast_cibs$CAFtype)), 
  col = caf_colors_breast[levels(as.factor(plot_breast_cibs$CAFtype))], 
  pch = 15, 
  cex = 0.75,
  y.intersp = 0.99,
  x.intersp = 0.5,
  bty = "n")

### Add CAF labels near each cluster
label_pos_breast <- aggregate(
  cbind(truth, estimated) ~ CAFtype,
  data = plot_breast_cibs,
  FUN = median)
# Default: place labels slightly above clusters
label_pos_breast$xpos <- label_pos_breast$truth
label_pos_breast$ypos <- label_pos_breast$estimated + 0.07
# Move only certain labels
label_pos_breast$ypos[label_pos_breast$CAFtype == "hsp_tpCAF"] <- 
  label_pos_breast$estimated[label_pos_breast$CAFtype == "hsp_tpCAF"] + 0.08
label_pos_breast$xpos[label_pos_breast$CAFtype == "IDO_CAF"] <- 
  label_pos_breast$truth[label_pos_breast$CAFtype == "IDO_CAF"] + 0.035
label_pos_breast$ypos[label_pos_breast$CAFtype == "IDO_CAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "IDO_CAF"] + 0.07
label_pos_breast$xpos[label_pos_breast$CAFtype == "tpCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "tpCAF"] + 0.035
label_pos_breast$ypos[label_pos_breast$CAFtype == "tpCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "tpCAF"] + 0.035
label_pos_breast$xpos[label_pos_breast$CAFtype == "apCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "apCAF"] + 0.037
label_pos_breast$ypos[label_pos_breast$CAFtype == "apCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "apCAF"] - 0.01
label_pos_breast$xpos[label_pos_breast$CAFtype == "rCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "rCAF"] - 0.025
label_pos_breast$ypos[label_pos_breast$CAFtype == "rCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "rCAF"] - 0.01
label_pos_breast$ypos[label_pos_breast$CAFtype == "Pericyte"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "Pericyte"] + 0.09
label_pos_breast$xpos[label_pos_breast$CAFtype == "vCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "vCAF"] + 0.02
label_pos_breast$ypos[label_pos_breast$CAFtype == "vCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "vCAF"] + 0.04
label_pos_breast$xpos[label_pos_breast$CAFtype == "iCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "iCAF"] - 0.008
label_pos_breast$ypos[label_pos_breast$CAFtype == "iCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "iCAF"] + 0.085
label_pos_breast$xpos[label_pos_breast$CAFtype == "mCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "mCAF"] - 0.015
text(
  x = label_pos_breast$xpos,
  y = label_pos_breast$ypos,
  labels = label_pos_breast$CAFtype,
  cex = 0.75,
  font = 1)

### LUNG
# Take the matching matrices
lung_truth_plot <- lung_cibs_eval$truth_matched
lung_est_plot <- lung_cibs_eval$estimated_matched

plot_lung_cibs <- data.frame(
  sample = rep(rownames(lung_truth_plot), times = ncol(lung_truth_plot)), 
  CAFtype = rep(colnames(lung_truth_plot), each = nrow(lung_truth_plot)), 
  truth = as.vector(as.matrix(lung_truth_plot)), 
  estimated = as.vector(as.matrix(lung_est_plot))) 
caf_colors_lung <- c(
  "apCAF" = "black",
  "iCAF" = "dodgerblue",
  "mCAF" = "mediumorchid3",
  "Pericyte" = "goldenrod2",
  "rCAF" = "grey60",
  "tpCAF" = "turquoise3",
  "vCAF" = "hotpink3")
plot(plot_lung_cibs$truth, 
     plot_lung_cibs$estimated, 
     xlab = "True cell-type proportions", 
     ylab = "Estimated cell-type proportions", 
     main = "CIBERSORTx deconvolution performance (lung data)", 
     pch = 16, 
     col = caf_colors_lung[plot_lung_cibs$CAFtype]) 

abline(0, 1, col = "red", lty = 2, lwd = 1.8)

### Add global correlation and RMSE to plot
metrics_text_lung <- paste0(
  "Global Pearson Correlation (r) = ", round(lung_cibs_eval$global_correlation, 3),
  "\nGlobal RMSE = ", round(lung_cibs_eval$global_rmse, 3))
text(
  x = 0.05,
  y = 0.46,
  labels = metrics_text_lung,
  adj = c(0, 1),
  cex = 0.95,
  font = 2)

legend(
  x = - 0.0055,
  y = 0.475,
  legend = levels(as.factor(plot_lung_cibs$CAFtype)), 
  col = caf_colors_lung[levels(as.factor(plot_lung_cibs$CAFtype))], 
  pch = 15, 
  cex = 0.75,
  y.intersp = 0.99,
  x.intersp = 0.5,
  bty = "n")

### Add CAF labels near each cluster
label_pos_lung <- aggregate(
  cbind(truth, estimated) ~ CAFtype,
  data = plot_lung_cibs,
  FUN = median)
# Default: place labels slightly above clusters
label_pos_lung$xpos <- label_pos_lung$truth
label_pos_lung$ypos <- label_pos_lung$estimated + 0.07
# Move only certain labels
label_pos_lung$xpos[label_pos_lung$CAFtype == "rCAF"] <-
  label_pos_lung$truth[label_pos_lung$CAFtype == "rCAF"] - 0.005
label_pos_lung$ypos[label_pos_lung$CAFtype == "rCAF"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "rCAF"] + 0.038
label_pos_lung$xpos[label_pos_lung$CAFtype == "vCAF"] <-
  label_pos_lung$truth[label_pos_lung$CAFtype == "vCAF"] + 0.02
label_pos_lung$ypos[label_pos_lung$CAFtype == "vCAF"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "vCAF"] + 0.1
label_pos_lung$xpos[label_pos_lung$CAFtype == "Pericyte"] <-
  label_pos_lung$truth[label_pos_lung$CAFtype == "Pericyte"] - 0.01
label_pos_lung$ypos[label_pos_lung$CAFtype == "Pericyte"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "Pericyte"] + 0.13
label_pos_lung$xpos[label_pos_lung$CAFtype == "tpCAF"] <-
  label_pos_lung$truth[label_pos_lung$CAFtype == "tpCAF"] - 0.012
label_pos_lung$xpos[label_pos_lung$CAFtype == "iCAF"] <-
  label_pos_lung$truth[label_pos_lung$CAFtype == "iCAF"] + 0.015
label_pos_lung$xpos[label_pos_lung$CAFtype == "mCAF"] <-
  label_pos_lung$truth[label_pos_lung$CAFtype == "mCAF"] - 0.015
label_pos_lung$ypos[label_pos_lung$CAFtype == "mCAF"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "mCAF"] + 0.08
text(
  x = label_pos_lung$xpos,
  y = label_pos_lung$ypos,
  labels = label_pos_lung$CAFtype,
  cex = 0.75,
  font = 1)

# Safety checks
all.equal(
  cor(plot_breast_cibs$truth, plot_breast_cibs$estimated, use = "complete.obs"),
  breast_cibs_eval$global_correlation)
all.equal(
  cor(plot_lung_cibs$truth, plot_lung_cibs$estimated, use = "complete.obs"),
  lung_cibs_eval$global_correlation)

## 9) Faceted plot per CAF type
library(ggplot2)
## BREAST
ggplot(plot_breast_cibs, aes(x = truth, y = estimated)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, colour = "red",
              linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ CAFtype, ncol = 5) +
  coord_equal(xlim = c(0, 0.42), ylim = c(0, 0.42)) +
  theme_bw() +
  labs(
    title = "CIBERSORTx performance by CAF type (breast data)",
    x = "True cell-type proportion",
    y = "Estimated cell-type proportion")
## LUNG
ggplot(plot_lung_cibs, aes(x = truth, y = estimated)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, colour = "red",
              linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ CAFtype, ncol = 4) +
  coord_equal(xlim = c(0, 0.42), ylim = c(0, 0.42)) +
  theme_bw() +
  labs(
    title = "CIBERSORTx performance by CAF type (lung data)",
    x = "True cell-type proportion",
    y = "Estimated cell-type proportion")
