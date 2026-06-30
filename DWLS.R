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
### 6) Build Signature Matrix for DWLS
breast_signature_dwls <- omnideconv::build_model(
  single_cell_object = breast_sc_counts2,
  cell_type_annotations = breast_cell_types,
  method = "dwls",
  dwls_method = "mast_optimized",
  pval_cutoff = 0.05,
  diff_cutoff = 0.5,
  ncores = 3,
  verbose = TRUE)
saveRDS(breast_signature_dwls, "dwls_results/breast_dwls_signature_FINAL.rds")

lung_signature_dwls <- omnideconv::build_model(
  single_cell_object = lung_sc_counts2,
  cell_type_annotations = lung_cell_types,
  method = "dwls",
  dwls_method = "mast_optimized",
  pval_cutoff = 0.05,
  diff_cutoff = 0.5,
  ncores = 3,
  verbose = TRUE)
saveRDS(lung_signature_dwls, "dwls_results/lung_dwls_signature_FINAL.rds")

### 7) Run DWLS Deconvolution
common_dwls_genes_breast <- intersect(rownames(breast_signature_dwls), rownames(breast_bulk_tpm2))
breast_dwls_est <- omnideconv::deconvolute(
  bulk_gene_expression = breast_bulk_tpm2[common_dwls_genes_breast, , drop = FALSE],
  model = breast_signature_dwls[common_dwls_genes_breast, , drop = FALSE],
  method = "dwls",
  dwls_submethod = "DampenedWLS",
  normalize_results = TRUE,
  verbose = TRUE)
write.csv(breast_dwls_est, "dwls_results/breast_dwls_estimates_FINAL.csv")

common_dwls_genes_lung <- intersect(rownames(lung_signature_dwls), rownames(lung_bulk_tpm2))
lung_dwls_est <- omnideconv::deconvolute(
  bulk_gene_expression = lung_bulk_tpm2[common_dwls_genes_lung, , drop = FALSE],
  model = lung_signature_dwls[common_dwls_genes_lung, , drop = FALSE],
  method = "dwls",
  dwls_submethod = "DampenedWLS",
  normalize_results = TRUE,
  verbose = TRUE)
write.csv(lung_dwls_est, "dwls_results/lung_dwls_estimates_FINAL.csv")

### 8) Performance Evaluation Function
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

### 9) Evaluate Breast Data
breast_dwls_eval <- evaluate_deconv(
  truth = breast_truth,
  estimated = breast_dwls_est)
breast_dwls_eval$global_correlation
breast_dwls_eval$correlation_by_celltype
breast_dwls_eval$global_rmse
breast_dwls_eval$rmse_by_celltype

### 10) Evaluate Lung Data
lung_dwls_eval <- evaluate_deconv(
  truth = lung_truth,
  estimated = lung_dwls_est)
lung_dwls_eval$global_correlation
lung_dwls_eval$correlation_by_celltype
lung_dwls_eval$global_rmse
lung_dwls_eval$rmse_by_celltype

### 11) Plot Results colored by CAF type
### BREAST
plot_breast <- data.frame(
  sample = rep(rownames(breast_truth), times = ncol(breast_truth)), 
  CAFtype = rep(colnames(breast_truth), each = nrow(breast_truth)), 
  truth = as.vector(as.matrix(breast_truth)), 
  estimated = as.vector(as.matrix(breast_dwls_est))) 
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
plot(plot_breast$truth, 
     plot_breast$estimated, 
     xlab = "True cell-type proportions", 
     ylab = "Estimated cell-type proportions", 
     main = "DWLS deconvolution performance (breast data)", 
     pch = 16, 
     col = caf_colors_breast[plot_breast$CAFtype]) 

abline(0, 1, col = "red", lty = 2, lwd = 1.8)

### Add global correlation and RMSE to plot
metrics_text_breast <- paste0(
  "Global Pearson Correlation (r) = ", round(breast_dwls_eval$global_correlation, 3),
  "\nGlobal RMSE = ", round(breast_dwls_eval$global_rmse, 3))
text(
  x = 0.09,
  y = 0.45,
  labels = metrics_text_breast,
  adj = c(0, 1),
  cex = 0.95,
  font = 2)

legend(
  x = 0.265,
  y = 0.465,
  legend = levels(as.factor(plot_breast$CAFtype)), 
  col = caf_colors_breast[levels(as.factor(plot_breast$CAFtype))], 
  pch = 15, 
  cex = 0.75,
  y.intersp = 0.99,
  x.intersp = 0.5,
  bty = "n")

### Add CAF labels near each cluster
label_pos_breast <- aggregate(
  cbind(truth, estimated) ~ CAFtype,
  data = plot_breast,
  FUN = median)
# Default: place labels slightly above clusters
label_pos_breast$xpos <- label_pos_breast$truth
label_pos_breast$ypos <- label_pos_breast$estimated + 0.07
# Move only certain labels
label_pos_breast$xpos[label_pos_breast$CAFtype == "apCAF"] <- 
  label_pos_breast$truth[label_pos_breast$CAFtype == "apCAF"] + 0.035
label_pos_breast$ypos[label_pos_breast$CAFtype == "apCAF"] <- 
  label_pos_breast$estimated[label_pos_breast$CAFtype == "apCAF"] + 0.005
label_pos_breast$xpos[label_pos_breast$CAFtype == "tpCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "tpCAF"] + 0.037
label_pos_breast$ypos[label_pos_breast$CAFtype == "tpCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "tpCAF"] + 0.01
label_pos_breast$xpos[label_pos_breast$CAFtype == "rCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "rCAF"] - 0.02
label_pos_breast$ypos[label_pos_breast$CAFtype == "rCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "rCAF"] + 0.03
label_pos_breast$xpos[label_pos_breast$CAFtype == "hsp_tpCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "hsp_tpCAF"] - 0.03
label_pos_breast$xpos[label_pos_breast$CAFtype == "Pericyte"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "Pericyte"] - 0.01
label_pos_breast$ypos[label_pos_breast$CAFtype == "vCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "vCAF"] + 0.11
label_pos_breast$xpos[label_pos_breast$CAFtype == "mCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "mCAF"] - 0.025
label_pos_breast$ypos[label_pos_breast$CAFtype == "mCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "mCAF"] + 0.09
text(
  x = label_pos_breast$xpos,
  y = label_pos_breast$ypos,
  labels = label_pos_breast$CAFtype,
  cex = 0.75,
  font = 1)

### LUNG
plot_lung <- data.frame(
  sample = rep(rownames(lung_truth), times = ncol(lung_truth)), 
  CAFtype = rep(colnames(lung_truth), each = nrow(lung_truth)), 
  truth = as.vector(as.matrix(lung_truth)), 
  estimated = as.vector(as.matrix(lung_dwls_est))) 
caf_colors_lung <- c(
  "apCAF" = "black",
  "iCAF" = "dodgerblue",
  "mCAF" = "mediumorchid3",
  "Pericyte" = "goldenrod2",
  "rCAF" = "grey60",
  "tpCAF" = "turquoise3",
  "vCAF" = "hotpink3")
plot(plot_lung$truth, 
     plot_lung$estimated, 
     xlab = "True cell-type proportions", 
     ylab = "Estimated cell-type proportions", 
     main = "DWLS deconvolution performance (lung data)", 
     pch = 16, 
     col = caf_colors_lung[plot_lung$CAFtype]) 

abline(0, 1, col = "red", lty = 2, lwd = 1.8)

### Add global correlation and RMSE to plot
metrics_text_lung <- paste0(
  "Global Pearson Correlation (r) = ", round(lung_dwls_eval$global_correlation, 3),
  "\nGlobal RMSE = ", round(lung_dwls_eval$global_rmse, 3))
text(
  x = 0.03,
  y = 0.44,
  labels = metrics_text_lung,
  adj = c(0, 1),
  cex = 0.95,
  font = 2)

legend(
  x = 0.31,
  y = 0.44,
  legend = levels(as.factor(plot_lung$CAFtype)), 
  col = caf_colors_lung[levels(as.factor(plot_lung$CAFtype))], 
  pch = 15, 
  cex = 0.75,
  y.intersp = 0.99,
  x.intersp = 0.5,
  bty = "n")

### Add CAF labels near each cluster
label_pos_lung <- aggregate(
  cbind(truth, estimated) ~ CAFtype,
  data = plot_lung,
  FUN = median)
# Default: place labels slightly above clusters
label_pos_lung$xpos <- label_pos_lung$truth
label_pos_lung$ypos <- label_pos_lung$estimated + 0.07
# Move only certain labels
label_pos_lung$ypos[label_pos_lung$CAFtype == "apCAF"] <- 
  label_pos_lung$estimated[label_pos_lung$CAFtype == "apCAF"] + 0.05
label_pos_lung$xpos[label_pos_lung$CAFtype == "vCAF"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "vCAF"] + 0.055
label_pos_lung$ypos[label_pos_lung$CAFtype == "vCAF"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "vCAF"] + 0.005
label_pos_lung$xpos[label_pos_lung$CAFtype == "rCAF"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "rCAF"] - 0.305
label_pos_lung$ypos[label_pos_lung$CAFtype == "rCAF"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "rCAF"] - 0.01
label_pos_lung$ypos[label_pos_lung$CAFtype == "Pericyte"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "Pericyte"] + 0.08
label_pos_lung$xpos[label_pos_lung$CAFtype == "tpCAF"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "tpCAF"] - 0.1
text(
  x = label_pos_lung$xpos,
  y = label_pos_lung$ypos,
  labels = label_pos_lung$CAFtype,
  cex = 0.75,
  font = 1)                         
