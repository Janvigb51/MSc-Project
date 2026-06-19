# InstaPrism Deconvolution with changed proportions
# Change scenarios in SimBu::simulate_bulk

# 1) Collect all CAF types
celltypes <- sort(unique(annot_table$cell_type))
celltypes

# 2) Collect original CAF proportions
base_props <- prop.table(table(annot_table$cell_type))
base_props <- base_props[celltypes]
base_props

# 3) Varying proportions with CUSTOM
target_cell <- "vCAF"
target_props <- seq(0.02, 0.50, by = 0.04)
replicates <- 10
custom_rows <- list()
for (p in target_props) {
  for (rep in seq_len(replicates)) {
    other_cells <- setdiff(celltypes, target_cell)
    
    # distribute the remaining proportion across other CAF types
    other_base <- base_props[other_cells]
    other_base <- other_base / sum(other_base)
    prop_vec <- setNames(rep(0, length(celltypes)), celltypes)
    prop_vec[target_cell] <- p
    prop_vec[other_cells] <- (1 - p) * other_base
    sample_name <- paste0(
      target_cell,
      "_",
      sprintf("%.2f", p),
      "_rep",
      rep
    )
    custom_rows[[sample_name]] <- prop_vec
  }
}
custom_scenario_data <- as.data.frame(do.call(rbind, custom_rows))

# 4) Simulation
simulation_vCAF_custom <- SimBu::simulate_bulk(
  data = cords_lung,
  scenario = "custom",
  custom_scenario_data = custom_scenario_data,
  ncells = 2000,
  nsamples = nrow(custom_scenario_data),
  scaling_factor = "read_number",
  seed = 20240618
)

# 5) Check Results
bulk_expr_vCAF <- SummarizedExperiment::assays(simulation_vCAF_custom$bulk)[["bulk_counts"]]
dim(bulk_expr_vCAF)
dim(simulation_vCAF_custom$cell_fractions)
head(simulation_vCAF_custom$cell_fractions)

# 6) Deconvolution
bulk_expr_vCAF <- as.matrix(bulk_expr_vCAF)
deconv_vCAF_custom <- InstaPrism(bulk_Expr = bulk_expr_vCAF, refPhi_cs = refPhi_obj)

# Predicted cell type proportions
estimated_vCAF_custom <- t(deconv_vCAF_custom@Post.ini.ct@theta)
head(estimated_vCAF_custom)

# Ground truth cell type proportions
truth_vCAF_custom <- simulation_vCAF_custom$cell_fractions
truth_vCAF_custom <- truth_vCAF_custom[
  rownames(estimated_vCAF_custom),
  colnames(estimated_vCAF_custom)]

# 7) Plot results
plot(
  truth_vCAF_custom[, "vCAF"],
  estimated_vCAF_custom[, "vCAF"],
  xlab = "True vCAF proportion",
  ylab = "Estimated vCAF proportion",
  main = "Deconvolution performance with varying vCAF",
  pch = 16)

abline(0, 1, col = "red", lty = 2, lwd = 2)

# 8) Performance Evaluation
# Pearson Correlation
vCAF_cor <- cor(truth_vCAF_custom[, "vCAF"],estimated_vCAF_custom[, "vCAF"])
vCAF_cor
# RMSE
vCAF_rmse <- sqrt(mean((truth_vCAF_custom[, "vCAF"] - estimated_vCAF_custom[, "vCAF"])^2))
vCAF_rmse
# Add metrics to plot
text(
  x = 0.045,
  y = 0.52,
  labels = paste0(
    "Pearson Correlation = ", round(vCAF_cor, 4),
    "\nRMSE = ", round(vCAF_rmse, 4)
  ),
  adj = 0
)

# 9) Faceted plots

## Option 1: Estimated proportion of other CAFs ################################
plot_vCAF_df <- data.frame(
  sample = rep(rownames(truth_vCAF_custom), times = ncol(truth_vCAF_custom)),
  CAFtype = rep(colnames(truth_vCAF_custom), each = nrow(truth_vCAF_custom)),
  
  # x-axis: the manipulated vCAF truth, repeated for every CAF panel
  true_vCAF = rep(truth_vCAF_custom[, "vCAF"], times = ncol(truth_vCAF_custom)),
  
  # y-axis: estimated proportion for each CAF type
  estimated = as.vector(as.matrix(estimated_vCAF_custom)),
  
  # optional: the actual true proportion of each CAF type in that sample
  true_CAF = as.vector(as.matrix(truth_vCAF_custom)))

library(ggplot2)
ggplot(plot_vCAF_df, aes(x = true_vCAF, y = estimated)) +
  geom_point(size = 1.8, alpha = 0.7) +
  facet_wrap(~ CAFtype, ncol = 5) +
  theme_bw() +
  labs(
    title = "Effect of changing vCAF on other CAFs",
    x = "True vCAF proportion",
    y = "Estimated CAF proportion")

## Option 2: True + Estimated proportion of other CAFs #########################

plot_vCAF_compare_df <- rbind(
  data.frame(
    sample = rep(rownames(truth_vCAF_custom), times = ncol(truth_vCAF_custom)),
    CAFtype = rep(colnames(truth_vCAF_custom), each = nrow(truth_vCAF_custom)),
    true_vCAF = rep(truth_vCAF_custom[, "vCAF"], times = ncol(truth_vCAF_custom)),
    value = as.vector(as.matrix(truth_vCAF_custom)),
    Source = "Truth"
  ),
  data.frame(
    sample = rep(rownames(estimated_vCAF_custom), times = ncol(estimated_vCAF_custom)),
    CAFtype = rep(colnames(estimated_vCAF_custom), each = nrow(estimated_vCAF_custom)),
    true_vCAF = rep(truth_vCAF_custom[, "vCAF"], times = ncol(estimated_vCAF_custom)),
    value = as.vector(as.matrix(estimated_vCAF_custom)),
    Source = "InstaPrism"
  )
)

ggplot(plot_vCAF_compare_df, aes(x = true_vCAF, y = value, colour = Source)) +
  geom_point(size = 1.6, alpha = 0.7) +
  stat_summary(aes(group = Source), fun = mean, geom = "line", linewidth = 0.9) +
  facet_wrap(~ CAFtype, ncol = 4) +
  theme_bw() +
  theme(
    legend.position = c(0.82, 0.32),
    legend.justification = c(0, 1),
    legend.background = element_blank()) +
  labs(
    title = "True vs estimated CAFs when vCAF is varied",
    x = "True vCAF proportion",
    y = "CAF proportion"
  )

############ Full plot #########################################################

# highest_vCAF <- max(truth_vCAF_custom[, "vCAF"])
# keep_samples <- rownames(truth_vCAF_custom)[
#   truth_vCAF_custom[, "vCAF"] == highest_vCAF]
# 
# # Full plot dataframe (all CAF types together)
# truth_high <- truth_vCAF_custom[keep_samples, , drop = FALSE]
# estimated_high <- estimated_vCAF_custom[keep_samples, , drop = FALSE]
# 
# plot_df_high <- data.frame(
#   truth = as.vector(as.matrix(truth_high)),
#   estimated = as.vector(as.matrix(estimated_high)),
#   CAFtype = rep(colnames(truth_high), each = nrow(truth_high))
# )
# 
# global_cor <- cor(plot_df_high$truth, plot_df_high$estimated)
# global_rmse <- sqrt(mean((plot_df_high$truth - plot_df_high$estimated)^2))
# 
# # Plot coloured by CAF type
# plot(
#   plot_df_high$truth,
#   plot_df_high$estimated,
#   xlab = "True cell-type proportions",
#   ylab = "Estimated cell-type proportions",
#   main = paste0(
#     "InstaPrism deconvolution performance\n(highest vCAF setting = ",
#     highest_vCAF, ")"
#   ),
#   pch = 16,
#   col = as.factor(plot_df_high$CAFtype)
# )
# 
# abline(0, 1, col = "red", lty = 2, lwd = 1.8)
# 
# legend(
#   "topleft",
#   legend = levels(as.factor(plot_df_high$CAFtype)),
#   col = seq_along(levels(as.factor(plot_df_high$CAFtype))),
#   pch = 15,
#   cex = 0.7
# )
# 
# text(
#   x = max(plot_df_high$truth) * 0.18,
#   y = max(plot_df_high$estimated) * 0.92,
#   labels = paste0(
#     "Global Pearson Correlation = ", round(global_cor, 4),
#     "\nRMSE = ", round(global_rmse, 4)
#   ),
#   adj = 0
# )
