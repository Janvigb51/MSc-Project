## CONTROLLED SENSITIVITY SIMULATION ##
# This script should be run once"
# Change scenarios in SimBu::simulate_bulk
# Use functions from Controlled_Simulation_Functions.R"
source("Controlled_Simulation_Functions.R")

### 1) Load the data/annotations
breast_annot <- readRDS("simbu_inputs/breast_annotation.rds")
lung_annot <- readRDS("simbu_inputs/lung_annotation.rds")

### 2) Make controlled scenarios for every CAF type in each dataset
breast_scenarios <- make_all_scenarios(
  annot_table = breast_annot,
  dataset_name = "Breast",
  target_props = seq(0.02, 0.50, by = 0.04),
  replicates = 10)

lung_scenarios <- make_all_scenarios(
  annot_table = lung_annot,
  dataset_name = "Lung",
  target_props = seq(0.02, 0.50, by = 0.04),
  replicates = 10)

## SAVE CONTROLLED SCENARIOS
controlled_dir <- "controlled_scenarios"
dir.create(controlled_dir, showWarnings = FALSE)
saveRDS(breast_scenarios, file.path(controlled_dir, "breast_all_controlled_scenarios.rds"))
saveRDS(lung_scenarios, file.path(controlled_dir, "lung_all_controlled_scenarios.rds"))

### 3) Simulate bulk samples with SimBu
## Load controlled scenarios
breast_scenarios <- readRDS("controlled_scenarios/breast_all_controlled_scenarios.rds")
lung_scenarios <- readRDS("controlled_scenarios/lung_all_controlled_scenarios.rds")
# Check available CAF types
names(breast_scenarios)
names(lung_scenarios)
# For example (get the controlled proportion tables where mCAF is gradually changed)
breast_scenarios[["mCAF"]]
lung_scenarios[["mCAF"]]
## Load simbu-ready datasets
cords_breast <- readRDS("simbu_inputs/cords_breast_simbu.rds")
cords_lung <- readRDS("simbu_inputs/cords_lung_simbu.rds")

### Run Simulation for all CAF scenarios automatically
breast_simulations <- simulate_all_scenarios(
  sc_data = cords_breast,
  scenario_list = breast_scenarios,
  ncells = 2000,
  seed = 20240618)

lung_simulations <- simulate_all_scenarios(
  sc_data = cords_lung,
  scenario_list = lung_scenarios,
  ncells = 2000,
  seed = 20240618)

## SAVE SIMULATED BULK OBJECTS
sim_dir <- "controlled_simulations"
dir.create(sim_dir, showWarnings = FALSE)
saveRDS(breast_simulations,
        file.path(sim_dir, "breast_controlled_simulations.rds"))
saveRDS(lung_simulations,
        file.path(sim_dir, "lung_controlled_simulations.rds"))

## Minor checks
# For each CAF we have:
# Simulated bulk expression to feed into deconvolution tools ($bulk)
# And true known CAF proportions to compare against deconvolution estimates ($cell_fractions)
breast_simulations[["mCAF"]]
lung_simulations[["mCAF"]]
