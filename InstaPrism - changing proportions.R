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
target_cell <- "mCAF"
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
simulation_mCAF_custom <- SimBu::simulate_bulk(
  data = cords_breast,
  scenario = "custom",
  custom_scenario_data = custom_scenario_data,
  ncells = 2000,
  nsamples = nrow(custom_scenario_data),
  scaling_factor = "read_number",
  seed = 20240618
)

# 5) Check Results
bulk_expr_mCAF <- SummarizedExperiment::assays(simulation_mCAF_custom$bulk)[["bulk_counts"]]
dim(bulk_expr_mCAF)
dim(simulation_mCAF_custom$cell_fractions)
head(simulation_mCAF_custom$cell_fractions)

# 6) Deconvolution
bulk_expr_mCAF <- as.matrix(bulk_expr_mCAF)
deconv_mCAF_custom <- InstaPrism(bulk_Expr = bulk_expr_mCAF, refPhi_cs = refPhi_obj)

# Predicted cell type proportions
estimated_mCAF_custom <- t(deconv_mCAF_custom@Post.ini.ct@theta)
head(estimated_mCAF_custom)

# Ground truth cell type proportions
truth_mCAF_custom <- simulation_mCAF_custom$cell_fractions
truth_mCAF_custom <- truth_mCAF_custom[
  rownames(estimated_mCAF_custom),
  colnames(estimated_mCAF_custom)]

# 7) Plot results
plot(
  truth_mCAF_custom[, "mCAF"],
  estimated_mCAF_custom[, "mCAF"],
  xlab = "True mCAF proportion",
  ylab = "Estimated mCAF proportion",
  main = "Deconvolution performance with varying mCAF",
  pch = 16)

abline(0, 1, col = "red", lty = 2, lwd = 2)

# 8) Performance Evaluation
# Pearson Correlation
mCAF_cor <- cor(truth_mCAF_custom[, "mCAF"],estimated_mCAF_custom[, "mCAF"])
mCAF_cor
# RMSE
mCAF_rmse <- sqrt(mean((truth_mCAF_custom[, "mCAF"] - estimated_mCAF_custom[, "mCAF"])^2))
mCAF_rmse
