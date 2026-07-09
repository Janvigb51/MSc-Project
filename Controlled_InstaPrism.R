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
})

# Load shared functions
source("Controlled_Simulation_Functions.R")

# Load SimBu-ready data
cords_breast <- readRDS("simbu_inputs/cords_breast_simbu.rds")
cords_lung <- readRDS("simbu_inputs/cords_lung_simbu.rds")

# Load controlled scenarios
breast_scenarios <- readRDS("controlled_scenarios/breast_all_controlled_scenarios.rds")
lung_scenarios <- readRDS("controlled_scenarios/lung_all_controlled_scenarios.rds")

# Load InstaPrism references
breast_refPhi_obj <- readRDS("instaprism_results/breast_refPhi_obj.rds")
lung_refPhi_obj <- readRDS("instaprism_results/lung_refPhi_obj.rds")

# Output folder
instaprism_controlled_dir <- "instaprism_results/controlled_spikein"
dir.create(instaprism_controlled_dir, recursive = TRUE, showWarnings = FALSE)