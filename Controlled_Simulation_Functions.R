## SHARED FUNCTIONS SCRIPT ##
# This script contains reusable functions for the controlled sensitivity simulations.

# make_scenario()
# Creates one custom SimBu proportion table where one CAF type is varied.

# make_all_scenarios()
# Repeats make_scenario() for every CAF type in a dataset.

# simulate_spikein_bulk()
# Uses SimBu to generate pseudobulk samples from one custom scenario.

# simulate_all_scenarios()
# Runs SimBu simulations for all CAF-type scenarios in a dataset.

# run_instaprism_deconv()
# Runs InstaPrism deconvolution on one simulated pseudobulk dataset.

# run_dwls_deconv()
# Runs DWLS deconvolution on one simulated pseudobulk dataset.

# format_estimates()
# Checks the deconvolution output orientation and transposes it if needed.

# align_truth_estimated()
# Aligns true and estimated CAF proportion matrices by sample and CAF type.

# evaluate_controlled_result()
# Calculates correlation, RMSE, and error for the controlled experiment.

# run_instaprism_scenario()
# Runs the full controlled workflow for one CAF type using InstaPrism.

# run_dwls_scenario()
# Runs the full controlled workflow for one CAF type using DWLS.

# plot_target_spikein()
# Plots true vs estimated proportions for the CAF type being varied.

# plot_all_caf_behaviour()
# Plots how all CAF types behave as the target CAF type is varied.

# save_controlled_result()
# Saves result objects, metric tables, and plots for one controlled experiment.

########################################################################################################

### MAKE CONTROLLED PSEUDOBULKS FUNCTION (where each CAF is varied)
make_scenario <- function(annot_table, target_cell, dataset_name = "dataset", 
                          target_props = seq(0.02, 0.50, by = 0.04), replicates = 10,
                          cell_type_col = "cell_type") {
  # Check that the annotation table has the correct CAF subtype column
  if (!cell_type_col %in% colnames(annot_table)) {
    stop(cell_type_col, " column was not found in annot_table.")}
  # Collect CAF types from whichever dataset you provide
  celltypes <- sort(unique(annot_table[[cell_type_col]]))
  # Check that the target CAF exists in this dataset
  if (!target_cell %in% celltypes) {
    stop(target_cell, " is not present in this dataset.")}
  # Check target proportions are valid
  if (any(target_props < 0 | target_props > 1)) {
    stop("All target proportions must be between 0 and 1.")}
  # Original CAF proportions in this dataset
  base_props <- prop.table(table(annot_table[[cell_type_col]]))
  base_props <- base_props[celltypes]
  custom_rows <- list()
  for (p in target_props) {
    for (rep in seq_len(replicates)) {
      other_cells <- setdiff(celltypes, target_cell)
      # Keep the relative balance of the other CAF types
      other_base <- base_props[other_cells]
      other_base <- other_base / sum(other_base)
      prop_vec <- setNames(rep(0, length(celltypes)), celltypes)
      prop_vec[target_cell] <- p
      prop_vec[other_cells] <- (1 - p) * other_base
      sample_name <- paste0(
        dataset_name, "_",
        target_cell, "_",
        sprintf("%.2f", p), "_",
        "rep", rep)
      custom_rows[[sample_name]] <- prop_vec}}
  custom_scenario_data <- as.data.frame(do.call(rbind, custom_rows))
  # Check every pseudobulk sample sums to 1
  if (any(abs(rowSums(custom_scenario_data) - 1) > 1e-8)) {
    stop("Some rows do not sum to 1.")}
  return(custom_scenario_data)}

## WRAPPER FUNCTION (controlled pseudobulks where all CAFs are varied)
make_all_scenarios <- function(annot_table, dataset_name,
target_props = seq(0.02, 0.50, by = 0.04), replicates = 10,
cell_type_col = "cell_type") {
  caf_types <- sort(unique(annot_table[[cell_type_col]]))
  scenario_list <- lapply(caf_types, function(ct) {
    make_scenario(
      annot_table = annot_table,
      target_cell = ct,
      dataset_name = dataset_name,
      target_props = target_props,
      replicates = replicates,
      cell_type_col = cell_type_col)})
  names(scenario_list) <- caf_types
  return(scenario_list)}

## SIMULATION FUNCTION
simulate_spikein_bulk <- function(sc_data, custom_scenario_data, ncells = 2000, seed = 20240618) {
  simulation_obj <- SimBu::simulate_bulk(
    data = sc_data,
    scenario = "custom",
    custom_scenario_data = custom_scenario_data,
    ncells = ncells,
    nsamples = nrow(custom_scenario_data),
    scaling_factor = "read_number",
    seed = seed)
  return(simulation_obj)}

## WRAPPER FUNCTION (for all scenarios)
simulate_all_scenarios <- function(sc_data, scenario_list, ncells = 2000, seed = 20240618) {
  simulation_list <- lapply(seq_along(scenario_list), function(i) {
    caf_name <- names(scenario_list)[i]
    message("Simulating: ", caf_name)
    simulate_spikein_bulk(
      sc_data = sc_data,
      custom_scenario_data = scenario_list[[caf_name]],
      ncells = ncells,
      seed = seed + i)})
  names(simulation_list) <- names(scenario_list)
  return(simulation_list)}

## INSTAPRISM DECONVOLUTION FUNCTION
run_instaprism_deconv <- function(simulation_obj, refPhi_obj) {
  # Extract simulated bulk counts
  bulk_expr <- SummarizedExperiment::assays(simulation_obj$bulk)[["bulk_counts"]]
  bulk_expr <- as.matrix(bulk_expr)
  # Run InstaPrism
  deconv_result <- InstaPrism(
    bulk_Expr = bulk_expr,
    refPhi_cs = refPhi_obj)
  # Extract estimated CAF proportions
  estimated <- t(deconv_result@Post.ini.ct@theta)
  return(list(
    bulk_expr = bulk_expr,
    deconv_result = deconv_result,
    estimated = estimated))}

## DWLS DECONVOLUTION FUNCTION
run_dwls_deconv <- function(simulation_obj, signature_matrix, verbose = TRUE) {
  # Extract simulated bulk counts
  bulk_counts <- SummarizedExperiment::assays(simulation_obj$bulk)[["bulk_counts"]]
  bulk_counts <- as.matrix(bulk_counts)
  mode(bulk_counts) <- "numeric"
  # Convert counts to CPM/TPM-like scale, same idea as original DWLS script
  bulk_tpm <- sweep(bulk_counts, 2, colSums(bulk_counts), "/") * 1e6
  bulk_tpm[is.na(bulk_tpm)] <- 0
  # Match genes between DWLS signature and simulated bulk
  common_genes <- intersect(rownames(signature_matrix), rownames(bulk_tpm))
  bulk_use <- bulk_tpm[common_genes, , drop = FALSE]
  model_use <- signature_matrix[common_genes, , drop = FALSE]
  # Run DWLS deconvolution
  deconv_result <- omnideconv::deconvolute(
    bulk_gene_expression = bulk_use,
    model = model_use,
    method = "dwls",
    dwls_submethod = "DampenedWLS",
    normalize_results = TRUE,
    verbose = verbose)
  estimated <- as.matrix(deconv_result)
  mode(estimated) <- "numeric"
  return(list(
    bulk_expr = bulk_use,
    deconv_result = deconv_result,
    estimated = estimated))}

## FORMAT DECONVOLUTION RESULTS FUNCTION
format_estimates <- function(estimated_raw, truth_matrix) {
  estimated <- as.matrix(estimated_raw)
  mode(estimated) <- "numeric"
  if (is.null(rownames(estimated)) || is.null(colnames(estimated))) {
    stop("Estimated matrix has missing rownames or colnames. Check DWLS output format.")}
  # Case 1: already samples x CAF types
  if (all(rownames(truth_matrix) %in% rownames(estimated)) &&
      all(colnames(truth_matrix) %in% colnames(estimated))) {
    return(estimated)}
  # Case 2: flipped, CAF types x samples
  if (all(rownames(truth_matrix) %in% colnames(estimated)) &&
      all(colnames(truth_matrix) %in% rownames(estimated))) {
    return(t(estimated))}
  stop("Estimated matrix does not match truth matrix. Check rownames/colnames.")}

## ALIGN TRUTH & ESTIMATED MATRICES FUNCTION
align_truth_estimated <- function(simulation_obj, estimated) {
  truth <- simulation_obj$cell_fractions
  truth <- as.matrix(truth)
  estimated <- as.matrix(estimated)
  common_samples <- intersect(rownames(truth), rownames(estimated))
  common_celltypes <- intersect(colnames(truth), colnames(estimated))
  if (length(common_samples) == 0) {
    stop("No matching sample names between truth and estimated matrices.")}
  if (length(common_celltypes) == 0) {
    stop("No matching CAF type names between truth and estimated matrices.")}
  truth_aligned <- truth[common_samples, common_celltypes, drop = FALSE]
  estimated_aligned <- estimated[common_samples, common_celltypes, drop = FALSE]
  estimated_aligned <- estimated_aligned[
    rownames(truth_aligned),
    colnames(truth_aligned),
    drop = FALSE]
  return(list(
    truth = truth_aligned,
    estimated = estimated_aligned))}

## EVALUATE CONTROLLED EXPERIMENT
safe_cor <- function(x, y) {
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    return(NA_real_)}
  cor(x, y, use = "complete.obs")}

evaluate_controlled_result <- function(truth, estimated, target_cell, dataset_name, method_name) {
  if (!target_cell %in% colnames(truth)) {
    stop(target_cell, " is not present in truth matrix.")}
  # Target CAF only
  target_df <- data.frame(
    Method = method_name,
    Dataset = dataset_name,
    Sample = rownames(truth),
    Target_CAF = target_cell,
    True = as.numeric(truth[, target_cell]),
    Estimated = as.numeric(estimated[, target_cell]))
  target_df$Error <- target_df$Estimated - target_df$True
  target_metrics <- data.frame(
    Method = method_name,
    Dataset = dataset_name,
    Target_CAF = target_cell,
    Correlation = safe_cor(target_df$True, target_df$Estimated),
    RMSE = sqrt(mean((target_df$Estimated - target_df$True)^2)),
    Mean_Error = mean(target_df$Error))
  # All CAF type metrics
  celltype_metrics <- bind_rows(lapply(colnames(truth), function(ct) {
    true_values <- as.numeric(truth[, ct])
    estimated_values <- as.numeric(estimated[, ct])
    data.frame(
      Method = method_name,
      Dataset = dataset_name,
      Target_CAF = target_cell,
      Measured_CAF = ct,
      Correlation = safe_cor(true_values, estimated_values),
      RMSE = sqrt(mean((estimated_values - true_values)^2)),
      Mean_Error = mean(estimated_values - true_values))}))
  # Long dataframe for plotting behaviour of all CAFs
  truth_long <- as.data.frame(truth) %>%
    rownames_to_column("Sample") %>%
    pivot_longer(
      cols = -Sample,
      names_to = "CAFtype",
      values_to = "Proportion") %>%
    mutate(Source = "Truth")
  estimated_long <- as.data.frame(estimated) %>%
    rownames_to_column("Sample") %>%
    pivot_longer(
      cols = -Sample,
      names_to = "CAFtype",
      values_to = "Proportion") %>%
    mutate(Source = method_name)
  target_map <- data.frame(
    Sample = rownames(truth),
    Target_True = as.numeric(truth[, target_cell]))
  eval_long <- bind_rows(truth_long, estimated_long) %>%
    left_join(target_map, by = "Sample") %>%
    mutate(
      Method = method_name,
      Dataset = dataset_name,
      Target_CAF = target_cell)
  return(list(
    target_df = target_df,
    target_metrics = target_metrics,
    celltype_metrics = celltype_metrics,
    eval_long = eval_long))}

## WRAPPER FUNCTION FOR 1 INSTAPRISM EXPERIMENT
run_instaprism_scenario <- function(custom_scenario_data, sc_data, refPhi_obj, target_cell,
                                    dataset_name, ncells = 2000, seed = 20240618) {
  # 1. Simulate bulk
  simulation_obj <- simulate_spikein_bulk(
    sc_data = sc_data,
    custom_scenario_data = custom_scenario_data,
    ncells = ncells,
    seed = seed)
  # 2. Run InstaPrism
  deconv <- run_instaprism_deconv(
    simulation_obj = simulation_obj,
    refPhi_obj = refPhi_obj)
  # 3. Align truth and estimates
  aligned <- align_truth_estimated(
    simulation_obj = simulation_obj,
    estimated = deconv$estimated)
  # 4. Evaluate
  evaluation <- evaluate_controlled_result(
    truth = aligned$truth,
    estimated = aligned$estimated,
    target_cell = target_cell,
    dataset_name = dataset_name,
    method_name = "InstaPrism")
  return(list(
    method = "InstaPrism",
    dataset = dataset_name,
    target_cell = target_cell,
    simulation = simulation_obj,
    bulk_expr = deconv$bulk_expr,
    deconv_result = deconv$deconv_result,
    truth = aligned$truth,
    estimated = aligned$estimated,
    target_df = evaluation$target_df,
    target_metrics = evaluation$target_metrics,
    celltype_metrics = evaluation$celltype_metrics,
    eval_long = evaluation$eval_long))}

## WRAPPER FUNCTION FOR 1 DWLS EXPERIMENT
run_dwls_scenario <- function(custom_scenario_data, sc_data, signature_matrix, target_cell,
                              dataset_name, ncells = 2000, seed = 20240618) {
  # 1. Simulate bulk
  simulation_obj <- simulate_spikein_bulk(
    sc_data = sc_data,
    custom_scenario_data = custom_scenario_data,
    ncells = ncells,
    seed = seed)
  # 2. Run DWLS
  deconv <- run_dwls_deconv(
    simulation_obj = simulation_obj,
    signature_matrix = signature_matrix,
    verbose = TRUE)
  # 3. Format DWLS estimates
  truth_raw <- as.matrix(simulation_obj$cell_fractions)
  estimated <- format_estimates(
    estimated_raw = deconv$estimated,
    truth_matrix = truth_raw)
  # 4. Align truth and estimates
  aligned <- align_truth_estimated(
    simulation_obj = simulation_obj,
    estimated = estimated)
  # 5. Evaluate
  evaluation <- evaluate_controlled_result(
    truth = aligned$truth,
    estimated = aligned$estimated,
    target_cell = target_cell,
    dataset_name = dataset_name,
    method_name = "DWLS")
  return(list(
    method = "DWLS",
    dataset = dataset_name,
    target_cell = target_cell,
    simulation = simulation_obj,
    bulk_expr = deconv$bulk_expr,
    deconv_result = deconv$deconv_result,
    truth = aligned$truth,
    estimated = aligned$estimated,
    target_df = evaluation$target_df,
    target_metrics = evaluation$target_metrics,
    celltype_metrics = evaluation$celltype_metrics,
    eval_long = evaluation$eval_long))}

## MAIN TARGET CAF SCATTER PLOT
plot_target_spikein <- function(scenario_result) {
  cor_value <- round(scenario_result$target_metrics$Correlation, 3)
  rmse_value <- round(scenario_result$target_metrics$RMSE, 3)
  ggplot(scenario_result$target_df, aes(x = True, y = Estimated)) +
    geom_point(alpha = 0.7, colour = "black") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
      label = paste0("Pearson Correlation = ", cor_value, "\nRMSE = ", rmse_value)) +
    labs(
      title = paste0("Deconvolution performance with varying ", scenario_result$target_cell),
      x = paste0("True ", scenario_result$target_cell, " proportion"),
      y = paste0("Estimated ", scenario_result$target_cell, " proportion")) +
    theme_bw()}

## BEHAVIOR OF ALL CAF TYPES
plot_all_caf_behaviour <- function(scenario_result) {
  ggplot(scenario_result$eval_long,
         aes(x = Target_True, y = Proportion, colour = Source)) +
    geom_point(alpha = 0.4) +
    stat_summary(aes(group = Source), fun = mean, geom = "line", linewidth = 0.8) +
    facet_wrap(~ CAFtype, scales = "free_y", ncol = 4) +
    scale_colour_manual(
      values = c(
        "InstaPrism" = "#F8766D",
        "DWLS" = "#F8766D",
        "Truth" = "#00BFC4")) +
    labs(
      title = paste0("True vs estimated CAFs when ", scenario_result$target_cell, " is varied"),
      x = paste0("True ", scenario_result$target_cell, " proportion"),
      y = "CAF proportion") +
    theme_bw()}

## SAVE CONTROLLED RESULT
save_controlled_result <- function(scenario_result, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dataset <- tolower(scenario_result$dataset)
  method <- tolower(scenario_result$method)
  target <- scenario_result$target_cell
  prefix <- paste(dataset, method, target, sep = "_")
  saveRDS(
    scenario_result,
    file.path(output_dir, paste0(prefix, "_full_result.rds")))
  write_csv(
    scenario_result$target_df,
    file.path(output_dir, paste0(prefix, "_target_values.csv")))
  write_csv(
    scenario_result$target_metrics,
    file.path(output_dir, paste0(prefix, "_target_metrics.csv")))
  write_csv(
    scenario_result$celltype_metrics,
    file.path(output_dir, paste0(prefix, "_all_celltype_metrics.csv")))
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_target_plot.png")),
    plot = plot_target_spikein(scenario_result),
    width = 6,
    height = 5,
    dpi = 300)
  ggsave(
    filename = file.path(output_dir, paste0(prefix, "_all_caf_behaviour_plot.png")),
    plot = plot_all_caf_behaviour(scenario_result),
    width = 10,
    height = 8,
    dpi = 300)}
