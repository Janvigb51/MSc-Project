# If deconvolution has already been run
# jump to QUICK START (line 105)

install.packages("gitcreds")
library(gitcreds)

gitcreds::gitcreds_delete()
Sys.unsetenv("GITHUB_PAT")
Sys.unsetenv("GITHUB_TOKEN")

# install InstaPrism once
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("Biobase", quietly = TRUE))
  BiocManager::install("Biobase")

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github(
  "humengying0907/InstaPrism",
  auth_token = NULL
)
library(InstaPrism)

### 1) Load in the Data
### We need the pseudobulk counts (genes across samples) into matrix with genes as row names
breast_bulk_expr <- read.table("C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-breast-100samples-2000cells-20240618/bulk_counts_2000cells_100samps.txt",
header = TRUE, sep = "\t", check.names = FALSE)
rownames(breast_bulk_expr) <- breast_bulk_expr$Gene
breast_bulk_expr <- as.matrix(breast_bulk_expr[, -1])
storage.mode(breast_bulk_expr) <- "numeric"

lung_bulk_expr <- read.table("C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-lung-100samples-2000cells-20240618/bulk_counts_2000cells_100samps.txt",
header = TRUE, sep = "\t", check.names = FALSE)
rownames(lung_bulk_expr) <- lung_bulk_expr$Gene
lung_bulk_expr <- as.matrix(lung_bulk_expr[, -1])
storage.mode(lung_bulk_expr) <- "numeric"

# check
dim(breast_bulk_expr) # 18093 genes in breast
breast_bulk_expr[1:5, 1:5]
dim(lung_bulk_expr) # 19142 genes in lung
lung_bulk_expr[1:5, 1:5]

### 2) Prepare the reference
set.seed(123)
breast_refPhi_obj <- InstaPrism::refPrepare(sc_Expr = sparse_matrix_breast, 
cell.type.labels = as.character(breast_annotation$cell_type),
cell.state.labels = as.character(breast_annotation$cell_type)) 
breast_refPhi_obj
set.seed(123)
lung_refPhi_obj <- InstaPrism::refPrepare(sc_Expr = sparse_matrix_lung, 
cell.type.labels = as.character(lung_annotation$cell_type),
cell.state.labels = as.character(lung_annotation$cell_type)) 
lung_refPhi_obj

### 3) Deconvolution with InstaPrism
set.seed(123)
breast_deconv = InstaPrism(bulk_Expr = breast_bulk_expr, refPhi_cs = breast_refPhi_obj)
set.seed(123)
lung_deconv = InstaPrism(bulk_Expr = lung_bulk_expr, refPhi_cs = lung_refPhi_obj)
#check output
class(breast_deconv)
slotNames(breast_deconv)

# Predicted cell type proportions
breast_est = t(breast_deconv@Post.ini.ct@theta)
head(breast_est)
lung_est = t(lung_deconv@Post.ini.ct@theta)
head(lung_est)

# SAVE RESULTS
results_dir <- "C:/Users/janvi/Desktop/MSc Project/instaprism_results"

saveRDS(breast_refPhi_obj, file.path(results_dir, "instaprism_breast_refphi_obj.rds"))
saveRDS(lung_refPhi_obj, file.path(results_dir, "instaprism_lung_refphi_obj.rds"))

saveRDS(breast_deconv, file.path(results_dir, "instaprism_breast_deconv.rds"))
saveRDS(lung_deconv, file.path(results_dir, "instaprism_lung_deconv.rds"))

saveRDS(breast_est, file.path(results_dir, "instaprism_breast_est.rds"))
saveRDS(lung_est, file.path(results_dir, "instaprism_lung_est.rds"))

# Ground truth cell type proportions
breast_truth <- read.table("C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-breast-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt",
header = TRUE, row.names = 1, check.names = FALSE)
breast_truth <- as.matrix(breast_truth)
storage.mode(breast_truth) <- "numeric"
dim(breast_truth)
head(breast_truth)

lung_truth <- read.table("C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-lung-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt",
header = TRUE, row.names = 1, check.names = FALSE)
lung_truth <- as.matrix(lung_truth)
storage.mode(lung_truth) <- "numeric"
dim(lung_truth)
head(lung_truth)

# SAVE RESULTS
saveRDS(breast_truth, file.path(results_dir, "instaprism_breast_truth.rds"))
saveRDS(lung_truth, file.path(results_dir, "instaprism_lung_truth.rds"))

############################################################
### QUICK START: Resume from saved InstaPrism results
### Use this section if deconvolution has already been run
### and the .rds result files already exist.
### This avoids rerunning InstaPrism deconvolution every time
### for performance evaluation plotting & benchmarking.
### (Comment out the below commands)
############################################################
# results_dir <- "C:/Users/janvi/Desktop/MSc Project/MSc Project R/instaprism_results"
# breast_est <- readRDS(file.path(results_dir, "instaprism_breast_est.rds"))
# lung_est <- readRDS(file.path(results_dir, "instaprism_lung_est.rds"))
# breast_truth <- readRDS(file.path(results_dir, "instaprism_breast_truth.rds"))
# lung_truth <- readRDS(file.path(results_dir, "instaprism_lung_truth.rds"))
############################################################

## 4) Performance Evaluation Function
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

## Evaluate Breast Data
breast_instaprism_eval <- evaluate_deconv(
  truth = breast_truth,
  estimated = breast_est)
breast_instaprism_eval$global_correlation
breast_instaprism_eval$correlation_by_celltype
breast_instaprism_eval$global_rmse
breast_instaprism_eval$rmse_by_celltype

## Evaluate Lung Data
lung_instaprism_eval <- evaluate_deconv(
  truth = lung_truth,
  estimated = lung_est)
lung_instaprism_eval$global_correlation
lung_instaprism_eval$correlation_by_celltype
lung_instaprism_eval$global_rmse
lung_instaprism_eval$rmse_by_celltype

### 11) Plot Results colored by CAF type
### BREAST
### Match truth and estimated matrices first
common_samples_breast <- intersect(rownames(breast_truth), rownames(breast_est))
common_celltypes_breast <- intersect(colnames(breast_truth), colnames(breast_est))
breast_truth_plot <- breast_truth[common_samples_breast, common_celltypes_breast, drop = FALSE]
breast_est_plot <- breast_est[common_samples_breast, common_celltypes_breast, drop = FALSE]
stopifnot(all(rownames(breast_truth_plot) == rownames(breast_est_plot)))
stopifnot(all(colnames(breast_truth_plot) == colnames(breast_est_plot)))

plot_breast <- data.frame(
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
plot(plot_breast$truth, 
     plot_breast$estimated, 
     xlab = "True cell-type proportions", 
     ylab = "Estimated cell-type proportions", 
     main = "InstaPrism deconvolution performance (breast data)", 
     pch = 16, 
     col = caf_colors_breast[plot_breast$CAFtype]) 

abline(0, 1, col = "red", lty = 2, lwd = 1.8)

### Add global correlation and RMSE to plot
metrics_text_breast <- paste0(
  "Global Pearson Correlation (r) = ", round(breast_instaprism_eval$global_correlation, 3),
  "\nGlobal RMSE = ", round(breast_instaprism_eval$global_rmse, 3))
text(
  x = 0.08,
  y = 0.4,
  labels = metrics_text_breast,
  adj = c(0, 1),
  cex = 0.95,
  font = 2)

legend(
  x = - 0.01,
  y = 0.42,
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
label_pos_breast$xpos[label_pos_breast$CAFtype == "mCAF"] <- 
  label_pos_breast$truth[label_pos_breast$CAFtype == "mCAF"] - 0.02
label_pos_breast$ypos[label_pos_breast$CAFtype == "mCAF"] <- 
  label_pos_breast$estimated[label_pos_breast$CAFtype == "mCAF"] + 0.05
label_pos_breast$xpos[label_pos_breast$CAFtype == "vCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "vCAF"] + 0.01
label_pos_breast$ypos[label_pos_breast$CAFtype == "vCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "vCAF"] + 0.04
label_pos_breast$xpos[label_pos_breast$CAFtype == "Pericyte"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "Pericyte"] - 0.01
label_pos_breast$ypos[label_pos_breast$CAFtype == "Pericyte"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "Pericyte"] + 0.04
label_pos_breast$xpos[label_pos_breast$CAFtype == "apCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "apCAF"] + 0.03
label_pos_breast$ypos[label_pos_breast$CAFtype == "apCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "apCAF"] + 0.035
label_pos_breast$xpos[label_pos_breast$CAFtype == "hsp_tpCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "hsp_tpCAF"] + 0.04
label_pos_breast$ypos[label_pos_breast$CAFtype == "hsp_tpCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "hsp_tpCAF"] + 0.002
label_pos_breast$xpos[label_pos_breast$CAFtype == "tpCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "tpCAF"] + 0.03
label_pos_breast$ypos[label_pos_breast$CAFtype == "tpCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "tpCAF"] + 0.03
label_pos_breast$ypos[label_pos_breast$CAFtype == "dCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "dCAF"] + 0.08
label_pos_breast$xpos[label_pos_breast$CAFtype == "rCAF"] <-
  label_pos_breast$truth[label_pos_breast$CAFtype == "rCAF"] - 0.01
label_pos_breast$ypos[label_pos_breast$CAFtype == "rCAF"] <-
  label_pos_breast$estimated[label_pos_breast$CAFtype == "rCAF"] + 0.01
text(
  x = label_pos_breast$xpos,
  y = label_pos_breast$ypos,
  labels = label_pos_breast$CAFtype,
  cex = 0.75,
  font = 1)

### LUNG
### Match truth and estimated matrices first
common_samples_lung <- intersect(rownames(lung_truth), rownames(lung_est))
common_celltypes_lung <- intersect(colnames(lung_truth), colnames(lung_est))
lung_truth_plot <- lung_truth[common_samples_lung, common_celltypes_lung, drop = FALSE]
lung_est_plot <- lung_est[common_samples_lung, common_celltypes_lung, drop = FALSE]
stopifnot(all(rownames(lung_truth_plot) == rownames(lung_est_plot)))
stopifnot(all(colnames(lung_truth_plot) == colnames(lung_est_plot)))

plot_lung <- data.frame(
  sample = rep(rownames(lung_truth_plot), times = ncol(lung_truth_plot)),
  CAFtype = rep(colnames(lung_truth_plot), each = nrow(lung_truth_plot)),
  truth = as.vector(as.matrix(lung_truth_plot)),
  estimated = as.vector(as.matrix(lung_est_plot)))
caf_colors_lung <- c(
  "apCAF" = "black",
  "iCAF" = "dodgerblue3",
  "mCAF" = "mediumorchid3",
  "Pericyte" = "goldenrod2",
  "rCAF" = "grey60",
  "tpCAF" = "turquoise3",
  "vCAF" = "hotpink3")
plot(plot_lung$truth,
     plot_lung$estimated,
     xlab = "True cell-type proportions",
     ylab = "Estimated cell-type proportions",
     main = "InstaPrism deconvolution performance (lung data)",
     pch = 16,
     col = caf_colors_lung[plot_lung$CAFtype])

abline(0, 1, col = "red", lty = 2, lwd = 1.8)

### Add global correlation and RMSE to plot
metrics_text_lung <- paste0(
  "Global Pearson Correlation (r) = ", round(lung_instaprism_eval$global_correlation, 3),
  "\nGlobal RMSE = ", round(lung_instaprism_eval$global_rmse, 3))
text(
  x = 0.06,
  y = 0.52,
  labels = metrics_text_lung,
  adj = c(0, 1),
  cex = 0.95,
  font = 2)

legend(
  x = - 0.01,
  y = 0.54,
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
label_pos_lung$ypos[label_pos_lung$CAFtype == "vCAF"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "vCAF"] + 0.05
label_pos_lung$xpos[label_pos_lung$CAFtype == "Pericyte"] <-
  label_pos_lung$truth[label_pos_lung$CAFtype == "Pericyte"] - 0.01
label_pos_lung$ypos[label_pos_lung$CAFtype == "Pericyte"] <-
  label_pos_lung$estimated[label_pos_lung$CAFtype == "Pericyte"] + 0.045
label_pos_lung$xpos[label_pos_lung$CAFtype == "iCAF"] <-
  label_pos_lung$truth[label_pos_lung$CAFtype == "iCAF"] + 0.027
text(
  x = label_pos_lung$xpos,
  y = label_pos_lung$ypos,
  labels = label_pos_lung$CAFtype,
  cex = 0.75,
  font = 1)

## 9) Faceted plot per CAF type
library(ggplot2)
## BREAST
ggplot(plot_breast, aes(x = truth, y = estimated)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, colour = "red",
              linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ CAFtype) +
  coord_equal(xlim = c(0, 0.42), ylim = c(0, 0.42)) +
  theme_bw() +
  labs(
    title = "InstaPrism performance by CAF type (breast data)",
    x = "True cell-type proportion",
    y = "Estimated cell-type proportion")
## LUNG
ggplot(plot_lung, aes(x = truth, y = estimated)) +
  geom_point(size = 1.8, alpha = 0.7) +
  geom_abline(intercept = 0, slope = 1, colour = "red",
              linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ CAFtype) +
  coord_equal(xlim = c(0, 0.42), ylim = c(0, 0.42)) +
  theme_bw() +
  labs(
    title = "InstaPrism performance by CAF type (lung data)",
    x = "True cell-type proportion",
    y = "Estimated cell-type proportion")
