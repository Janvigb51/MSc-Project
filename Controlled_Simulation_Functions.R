## SHARED FUNCTIONS SCRIPT ##

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
run_dwls_deconv <- function(simulation_obj, signature_matrix, container = "apptainer",
container_path = "/data/containers/", verbose = TRUE) {
  # Extract simulated bulk counts
  bulk_expr <- SummarizedExperiment::assays(simulation_obj$bulk)[["bulk_counts"]]
  bulk_expr <- as.matrix(bulk_expr)
  # Run DWLS
  deconv_result <- omnideconv::deconvolute(
    bulk_gene_expression = bulk_expr,
    signature = signature_matrix,
    method = "dwls",
    container = container,
    container_path = container_path,
    verbose = verbose)
  return(list(
    bulk_expr = bulk_expr,
    deconv_result = deconv_result))}

## FORMAT DECONVOLUTION RESULTS FUNCTION
format_estimates <- function(estimated_raw, truth_matrix) {
  estimated <- as.matrix(estimated_raw)
  mode(estimated) <- "numeric"
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
