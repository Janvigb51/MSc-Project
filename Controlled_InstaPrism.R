## CONTROLLED SENSITIVITY SIMULATION ##
# InstaPrism deconvolution with changed proportions

suppressPackageStartupMessages({
  library(SimBu)
  library(SummarizedExperiment)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(InstaPrism)})

# Load shared functions
source("Controlled_Simulation_Functions.R")

# Load SimBu-ready data
cords_breast <- readRDS("simbu_inputs/cords_breast_simbu.rds")
cords_lung <- readRDS("simbu_inputs/cords_lung_simbu.rds")

# Load controlled scenarios
breast_scenarios <- readRDS("controlled_scenarios/breast_all_controlled_scenarios.rds")
lung_scenarios <- readRDS("controlled_scenarios/lung_all_controlled_scenarios.rds")

# Load InstaPrism references
breast_refPhi_obj <- readRDS("instaprism_results/instaprism_breast_refphi_obj.rds")
lung_refPhi_obj <- readRDS("instaprism_results/instaprism_lung_refphi_obj.rds")

# Output folder
instaprism_controlled_dir <- "instaprism_results/controlled_spikein"
dir.create(instaprism_controlled_dir, recursive = TRUE, showWarnings = FALSE)

## Optional: Change to TRUE and test 1 controlled scenario first (lung mCAF)
## If it runs sucessfully change to FALSE to ignore later.
if (FALSE) {
  lung_mCAF_instaprism <- run_instaprism_scenario(
    custom_scenario_data = lung_scenarios[["mCAF"]],
    sc_data = cords_lung,
    refPhi_obj = lung_refPhi_obj,
    target_cell = "mCAF",
    dataset_name = "Lung",
    ncells = 2000,
    seed = 20240618)
# Check result
  lung_mCAF_instaprism$target_metrics
  head(lung_mCAF_instaprism$target_df)
  dim(lung_mCAF_instaprism$truth)
  dim(lung_mCAF_instaprism$estimated)
# Plot result
  plot_target_spikein(lung_mCAF_instaprism)
  plot_all_caf_behaviour(lung_mCAF_instaprism)
# Save Result
  save_controlled_result(
    scenario_result = lung_mCAF_instaprism,
    output_dir = instaprism_controlled_dir)}

## Run all InstaPrism controlled scenarios
run_all_instaprism_scenarios <- function(scenario_list, sc_data, refPhi_obj, dataset_name,
                                         output_dir, ncells = 2000, seed = 20240618) {
  results <- list()
  for (target_cell in names(scenario_list)) {
    message("Running InstaPrism ", dataset_name, " target CAF: ", target_cell)
    scenario_result <- run_instaprism_scenario(
      custom_scenario_data = scenario_list[[target_cell]],
      sc_data = sc_data,
      refPhi_obj = refPhi_obj,
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
breast_instaprism_controlled <- run_all_instaprism_scenarios(
  scenario_list = breast_scenarios,
  sc_data = cords_breast,
  refPhi_obj = breast_refPhi_obj,
  dataset_name = "Breast",
  output_dir = instaprism_controlled_dir,
  ncells = 2000,
  seed = 20240618)

# FOR LUNG
lung_instaprism_controlled <- run_all_instaprism_scenarios(
  scenario_list = lung_scenarios,
  sc_data = cords_lung,
  refPhi_obj = lung_refPhi_obj,
  dataset_name = "Lung",
  output_dir = instaprism_controlled_dir,
  ncells = 2000,
  seed = 20240618)

# Save full lists
saveRDS(
  breast_instaprism_controlled,
  file.path(instaprism_controlled_dir, "breast_instaprism_all_controlled_results.rds"))
saveRDS(
  lung_instaprism_controlled,
  file.path(instaprism_controlled_dir, "lung_instaprism_all_controlled_results.rds"))

# Combine InstaPrism controlled metrics
instaprism_target_metrics <- do.call(
  rbind, c(
    lapply(breast_instaprism_controlled, function(x) x$target_metrics),
    lapply(lung_instaprism_controlled, function(x) x$target_metrics)))
instaprism_target_metrics <- as.data.frame(instaprism_target_metrics)

instaprism_celltype_metrics <- do.call(
  rbind, c(
    lapply(breast_instaprism_controlled, function(x) x$celltype_metrics),
    lapply(lung_instaprism_controlled, function(x) x$celltype_metrics)))
instaprism_celltype_metrics <- as.data.frame(instaprism_celltype_metrics)

write_csv(
  instaprism_target_metrics,
  file.path(instaprism_controlled_dir, "instaprism_controlled_target_metrics.csv"))

write_csv(
  instaprism_celltype_metrics,
  file.path(instaprism_controlled_dir, "instaprism_controlled_all_celltype_metrics.csv"))
