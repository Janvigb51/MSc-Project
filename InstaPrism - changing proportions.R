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
target_cell <- "IDO_CAF"
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
simulation_IDO_CAF_custom <- SimBu::simulate_bulk(
  data = cords_breast,
  scenario = "custom",
  custom_scenario_data = custom_scenario_data,
  ncells = 2000,
  nsamples = nrow(custom_scenario_data),
  scaling_factor = "read_number",
  seed = 20240618
)

# 5) Check Results
bulk_expr_IDO_CAF <- SummarizedExperiment::assays(simulation_IDO_CAF_custom$bulk)[["bulk_counts"]]
dim(bulk_expr_IDO_CAF)
dim(simulation_IDO_CAF_custom$cell_fractions)
head(simulation_IDO_CAF_custom$cell_fractions)

# 6) Deconvolution
bulk_expr_IDO_CAF <- as.matrix(bulk_expr_IDO_CAF)
deconv_IDO_CAF_custom <- InstaPrism(bulk_Expr = bulk_expr_IDO_CAF, refPhi_cs = refPhi_obj)

# Predicted cell type proportions
estimated_IDO_CAF_custom <- t(deconv_IDO_CAF_custom@Post.ini.ct@theta)
head(estimated_IDO_CAF_custom)

# Ground truth cell type proportions
truth_IDO_CAF_custom <- simulation_IDO_CAF_custom$cell_fractions
truth_IDO_CAF_custom <- truth_IDO_CAF_custom[
  rownames(estimated_IDO_CAF_custom),
  colnames(estimated_IDO_CAF_custom)]

# 7) Plot results
plot(
  truth_IDO_CAF_custom[, "IDO_CAF"],
  estimated_IDO_CAF_custom[, "IDO_CAF"],
  xlab = "True IDO_CAF proportion",
  ylab = "Estimated IDO_CAF proportion",
  main = "Deconvolution performance with varying IDO_CAF",
  pch = 16)

abline(0, 1, col = "red", lty = 2, lwd = 2)

# 8) Performance Evaluation
# Pearson Correlation
IDO_CAF_cor <- cor(truth_IDO_CAF_custom[, "IDO_CAF"],estimated_IDO_CAF_custom[, "IDO_CAF"])
IDO_CAF_cor
# RMSE
IDO_CAF_rmse <- sqrt(mean((truth_IDO_CAF_custom[, "IDO_CAF"] - estimated_IDO_CAF_custom[, "IDO_CAF"])^2))
IDO_CAF_rmse
# Add metrics to plot
text(
  x = 0.045,
  y = 0.52,
  labels = paste0(
    "Pearson Correlation = ", round(IDO_CAF_cor, 4),
    "\nRMSE = ", round(IDO_CAF_rmse, 4)
  ),
  adj = 0
)

# 9) Faceted plots

## Option 1: Estimated proportion of other CAFs ################################
plot_IDO_CAF_df <- data.frame(
  sample = rep(rownames(truth_IDO_CAF_custom), times = ncol(truth_IDO_CAF_custom)),
  CAFtype = rep(colnames(truth_IDO_CAF_custom), each = nrow(truth_IDO_CAF_custom)),
  
  # x-axis: the manipulated IDO_CAF truth, repeated for every CAF panel
  true_IDO_CAF = rep(truth_IDO_CAF_custom[, "IDO_CAF"], times = ncol(truth_IDO_CAF_custom)),
  
  # y-axis: estimated proportion for each CAF type
  estimated = as.vector(as.matrix(estimated_IDO_CAF_custom)),
  
  # optional: the actual true proportion of each CAF type in that sample
  true_CAF = as.vector(as.matrix(truth_IDO_CAF_custom)))

library(ggplot2)
ggplot(plot_IDO_CAF_df, aes(x = true_IDO_CAF, y = estimated)) +
  geom_point(size = 1.8, alpha = 0.7) +
  facet_wrap(~ CAFtype, ncol = 5) +
  theme_bw() +
  labs(
    title = "Effect of changing IDO_CAF on other CAFs",
    x = "True IDO_CAF proportion",
    y = "Estimated CAF proportion")

## Option 2: True + Estimated proportion of other CAFs #########################

plot_IDO_CAF_compare_df <- rbind(
  data.frame(
    sample = rep(rownames(truth_IDO_CAF_custom), times = ncol(truth_IDO_CAF_custom)),
    CAFtype = rep(colnames(truth_IDO_CAF_custom), each = nrow(truth_IDO_CAF_custom)),
    true_IDO_CAF = rep(truth_IDO_CAF_custom[, "IDO_CAF"], times = ncol(truth_IDO_CAF_custom)),
    value = as.vector(as.matrix(truth_IDO_CAF_custom)),
    Source = "Truth"
  ),
  data.frame(
    sample = rep(rownames(estimated_IDO_CAF_custom), times = ncol(estimated_IDO_CAF_custom)),
    CAFtype = rep(colnames(estimated_IDO_CAF_custom), each = nrow(estimated_IDO_CAF_custom)),
    true_IDO_CAF = rep(truth_IDO_CAF_custom[, "IDO_CAF"], times = ncol(estimated_IDO_CAF_custom)),
    value = as.vector(as.matrix(estimated_IDO_CAF_custom)),
    Source = "InstaPrism"
  )
)

ggplot(plot_IDO_CAF_compare_df, aes(x = true_IDO_CAF, y = value, colour = Source)) +
  geom_point(size = 1.6, alpha = 0.7) +
  stat_summary(aes(group = Source), fun = mean, geom = "line", linewidth = 0.9) +
  facet_wrap(~ CAFtype, ncol = 5) +
  theme_bw() +
  labs(
    title = "True vs estimated CAFs when IDO_CAF is varied",
    x = "True IDO_CAF proportion",
    y = "CAF proportion"
  )

############ Full plot #########################################################

# highest_IDO_CAF <- max(truth_IDO_CAF_custom[, "IDO_CAF"])
# keep_samples <- rownames(truth_IDO_CAF_custom)[
#   truth_IDO_CAF_custom[, "IDO_CAF"] == highest_IDO_CAF]
# 
# # Full plot dataframe (all CAF types together)
# truth_high <- truth_IDO_CAF_custom[keep_samples, , drop = FALSE]
# estimated_high <- estimated_IDO_CAF_custom[keep_samples, , drop = FALSE]
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
#     "InstaPrism deconvolution performance\n(highest IDO_CAF setting = ",
#     highest_IDO_CAF, ")"
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
