## CONTROLLED SENSITIVITY SIMULATION ##
# DWLS deconvolution with changed proportions

suppressPackageStartupMessages({
  library(SimBu)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(omnideconv)})

# Load shared functions
source("Controlled_Simulation_Functions.R")

# Load SimBu-ready data
cords_breast <- readRDS("simbu_inputs/cords_breast_simbu.rds")
cords_lung <- readRDS("simbu_inputs/cords_lung_simbu.rds")

# Load controlled scenarios
breast_scenarios <- readRDS("controlled_scenarios/breast_all_controlled_scenarios.rds")
lung_scenarios <- readRDS("controlled_scenarios/lung_all_controlled_scenarios.rds")

# Load DWLS signatures from original DWLS benchmark
breast_signature_dwls <- readRDS("dwls_results/breast_dwls_signature_FINAL.rds")
lung_signature_dwls <- readRDS("dwls_results/lung_dwls_signature_FINAL.rds")

# Output folder
dwls_controlled_dir <- "dwls_results/controlled_spikein"
dir.create(dwls_controlled_dir, recursive = TRUE, showWarnings = FALSE)

## Optional: Change to TRUE and test 1 controlled scenario first
## If it runs successfully, change back to FALSE.
if (FALSE) {
  lung_mCAF_dwls <- run_dwls_scenario(
    custom_scenario_data = lung_scenarios[["mCAF"]],
    sc_data = cords_lung,
    signature_matrix = lung_signature_dwls,
    target_cell = "mCAF",
    dataset_name = "Lung",
    ncells = 2000,
    seed = 20240618)
  # Check result
  lung_mCAF_dwls$target_metrics
  head(lung_mCAF_dwls$target_df)
  dim(lung_mCAF_dwls$truth)
  dim(lung_mCAF_dwls$estimated)
  # Plot result
  plot_target_spikein(lung_mCAF_dwls)
  plot_all_caf_behaviour(lung_mCAF_dwls)
  # Save result
  save_controlled_result(
    scenario_result = lung_mCAF_dwls,
    output_dir = dwls_controlled_dir)}

## Run all DWLS controlled scenarios
run_all_dwls_scenarios <- function(scenario_list, sc_data, signature_matrix, dataset_name, 
                                   output_dir, ncells = 2000, seed = 20240618) {
  results <- list()
  for (target_cell in names(scenario_list)) {
    message("Running DWLS ", dataset_name, " target CAF: ", target_cell)
    scenario_result <- run_dwls_scenario(
      custom_scenario_data = scenario_list[[target_cell]],
      sc_data = sc_data,
      signature_matrix = signature_matrix,
      target_cell = target_cell,
      dataset_name = dataset_name,
      ncells = ncells,
      seed = seed)
    save_controlled_result(
      scenario_result = scenario_result,
      output_dir = output_dir)
    results[[target_cell]] <- scenario_result}
  return(results)}

# FOR BREAST
breast_dwls_controlled <- run_all_dwls_scenarios(
  scenario_list = breast_scenarios,
  sc_data = cords_breast,
  signature_matrix = breast_signature_dwls,
  dataset_name = "Breast",
  output_dir = dwls_controlled_dir,
  ncells = 2000,
  seed = 20240618)

# FOR LUNG
lung_dwls_controlled <- run_all_dwls_scenarios(
  scenario_list = lung_scenarios,
  sc_data = cords_lung,
  signature_matrix = lung_signature_dwls,
  dataset_name = "Lung",
  output_dir = dwls_controlled_dir,
  ncells = 2000,
  seed = 20240618)

# Save full lists
saveRDS(
  breast_dwls_controlled,
  file.path(dwls_controlled_dir, "breast_dwls_all_controlled_results.rds"))
saveRDS(
  lung_dwls_controlled,
  file.path(dwls_controlled_dir, "lung_dwls_all_controlled_results.rds"))

# Combine DWLS controlled metrics
dwls_target_metrics <- do.call(
  rbind, c(
    lapply(breast_dwls_controlled, function(x) x$target_metrics),
    lapply(lung_dwls_controlled, function(x) x$target_metrics)))
dwls_target_metrics <- as.data.frame(dwls_target_metrics)

dwls_celltype_metrics <- do.call(
  rbind, c(
    lapply(breast_dwls_controlled, function(x) x$celltype_metrics),
    lapply(lung_dwls_controlled, function(x) x$celltype_metrics)))
dwls_celltype_metrics <- as.data.frame(dwls_celltype_metrics)

write.csv(
  dwls_target_metrics,
  file.path(dwls_controlled_dir, "dwls_controlled_target_metrics.csv"),
  row.names = FALSE)

write.csv(
  dwls_celltype_metrics,
  file.path(dwls_controlled_dir, "dwls_controlled_all_celltype_metrics.csv"),
  row.names = FALSE)
