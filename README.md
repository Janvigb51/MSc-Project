# Benchmarking cell-type deconvolution methods for cancer-associated fibroblast (CAF) subtypes

MSc research project, MSc Genomic Data Science, University of Galway (2025–2026).

This repository benchmarks three cell-type deconvolution methods, **DWLS**, **InstaPrism** and **CIBERSORTx**, for estimating the abundance of closely related cancer-associated fibroblast (CAF) subtypes from bulk RNA-seq data in **breast** and **lung** cancer. The methods are first evaluated on simulated pseudo-bulk samples with known cell-type proportions, then stress-tested in a controlled sensitivity analysis in which the abundance of each CAF subtype is increased step by step from 2% to 50%.

## Key findings

- **DWLS performed best** in both cancer types, followed by InstaPrism, with CIBERSORTx the least accurate.
- The controlled sensitivity analysis supported this ranking: DWLS tracked gradual changes in CAF abundance most faithfully, particularly in lung cancer.
- Accuracy depended strongly on CAF subtype and dataset. Breast iCAFs and lung mCAFs were estimated well, whereas the closely related tpCAF / hsp tpCAF pair and vCAFs / pericytes were consistently underestimated and hard to separate.
- A high global correlation does not guarantee accurate estimates for every subtype, so global and per-subtype metrics should be read together.

### Primary benchmark (100 pseudo-bulk samples per cancer type)

| Method | Dataset | Global Pearson r | Global RMSE |
|---|---|---|---|
| DWLS | Breast | 0.923 | 0.049 |
| DWLS | Lung | 0.956 | 0.038 |
| InstaPrism | Breast | 0.793 | 0.060 |
| InstaPrism | Lung | 0.848 | 0.082 |
| CIBERSORTx | Breast | 0.606 | 0.083 |
| CIBERSORTx | Lung | 0.540 | 0.109 |

CIBERSORTx was excluded from the controlled sensitivity analysis because of its weakest performance in the primary benchmark.

## Study design

**Reference data.** Single-cell RNA-seq datasets of tumour fibroblasts from [Cords et al., *Nature Communications* 14, 4294 (2023)](https://doi.org/10.1038/s41467-023-39762-1), stored as Seurat objects:

- **Breast:** 16,704 cells from 14 patients (14,315 CAFs in 9 subtypes plus 2,389 pericytes); labels in the `CAFtype` metadata column; 18,984 genes.
- **Lung:** 1,377 cells (6 CAF subtypes plus pericytes); labels in `cluster_ft`; 25,765 genes.

The data files are not included in this repository.

**Primary benchmark.** Pseudo-bulk samples were generated with [SimBu](https://github.com/omnideconv/SimBu): 100 samples per cancer type, 2,000 cells per sample, `mirror_db` scenario (proportions mirror those of the source dataset), random seed 20240618. The known proportions serve as ground truth.

**Controlled sensitivity analysis.** Each CAF subtype is targeted in turn and raised from 2% to 50% in steps of 4 percentage points (13 abundance levels) with 10 replicates per level, i.e. 130 pseudo-bulk samples per target population (2,000 cells each, SimBu `custom` scenario, seed 20240618). The remaining cell types are rescaled so that every composition sums to one. This yields 1,300 pseudo-bulk profiles for breast cancer (10 populations) and 910 for lung cancer (7 populations), analysed with InstaPrism and DWLS.

**Evaluation metrics.** Pearson correlation and RMSE between estimated and true proportions, reported globally and per cell type. The controlled analysis additionally reports mean error (positive = overestimation, negative = underestimation) for each target population as its abundance increases.

**Methods.** DWLS (DampenedWLS, run through [omnideconv](https://github.com/omnideconv/omnideconv)), InstaPrism (fast, derandomised implementation of BayesPrism) and CIBERSORTx.

## Repository structure

| File / folder | Purpose |
|---|---|
| `SimBu.R` | Prepares the single-cell datasets and generates the primary breast and lung pseudo-bulk samples |
| `InstaPrism.R` | Builds InstaPrism references, runs deconvolution and evaluates the estimates |
| `DWLS.R` | Builds DWLS signature matrices, runs deconvolution and evaluates the estimates |
| `CIBERSORTx.R` | Prepares CIBERSORTx inputs, processes its outputs and evaluates the estimates |
| `Compare_all_methods.R` | Combines the performance metrics of the three methods into one comparison |
| `Prepare_Controlled_Simulations.R` | Generates the controlled CAF-abundance scenarios and their pseudo-bulk simulations |
| `Controlled_Simulation_Functions.R` | Reusable functions for scenario generation, simulation, deconvolution, evaluation and plotting |
| `Controlled_InstaPrism.R` | Runs the controlled sensitivity workflow with InstaPrism |
| `Controlled_DWLS.R` | Runs the controlled sensitivity workflow with DWLS |
| `cibersortx/` | CIBERSORTx inputs, container-related files and outputs |
| `instaprism_results/` | InstaPrism reference objects, estimates, metrics and figures |
| `dwls_results/` | DWLS signature matrices, estimates, metrics and figures |
| `combined_results/` | Combined benchmarking metrics (`combined_global_metrics.csv`, `combined_celltype_metrics.csv`) |
| `figures/controlled_analysis/` | Selected figures from the controlled sensitivity experiments |

## How to run

Run the scripts in this order:

1. `SimBu.R`: generates the primary pseudo-bulk samples and the ground-truth proportions.
2. `DWLS.R` and `InstaPrism.R`: run each method on the primary benchmark (any order).
3. `CIBERSORTx.R`: prepares the CIBERSORTx inputs and processes its results (see below).
4. `Compare_all_methods.R`: combines the metrics from steps 2 and 3 into `combined_results/`.
5. `Prepare_Controlled_Simulations.R`, then `Controlled_InstaPrism.R` and `Controlled_DWLS.R`: the controlled sensitivity analysis.

Update the file paths at the top of each script to match your own directory layout before running.

**CIBERSORTx.** CIBERSORTx was run outside R, on the Lugh HPC cluster (University of Galway), using the CIBERSORTx Fractions Singularity container submitted as SLURM jobs. Reproducing this step requires a CIBERSORTx registration and access token from the [CIBERSORTx website](https://cibersortx.stanford.edu/).

## Requirements

- R and the RStudio project file `MSc Project.Rproj`
- CRAN / Bioconductor: `Seurat`, `SeuratObject`, `SimBu`, `SummarizedExperiment`, `BiocParallel`, `Biobase`, `Matrix`, `dplyr`, `tidyr`, `tibble`, `readr`, `ggplot2`
- `omnideconv` (used for DWLS and for data handling)
- `InstaPrism`, installed from GitHub: `devtools::install_github("humengying0907/InstaPrism")`

## Limitations

- The same single-cell datasets were used to build both the references and the pseudo-bulk mixtures, so performance on independent bulk RNA-seq samples (other patients, platforms or studies) may be lower.
- The breast and lung datasets differ in size, number of cell types and source, so differences between them cannot be attributed to cancer type alone.
- The pseudo-bulk samples contain only fibroblast-lineage populations, without malignant, immune or endothelial cells.
- CIBERSORTx was not included in the controlled sensitivity analysis, so its ability to track gradual changes was not tested.
- Because estimated fractions sum to one, underestimation of one population is often a side effect of overestimation of another.

## Reference

Cords L, Tietscher S, Anzeneder T, et al. Cancer-associated fibroblast classification in single-cell and spatial proteomics data. *Nature Communications* 14, 4294 (2023). https://doi.org/10.1038/s41467-023-39762-1

## Author

Janvi Gianchandani Baxani, MSc Genomic Data Science, University of Galway.
