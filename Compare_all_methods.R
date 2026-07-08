# After all deconvolution tools have been run
# Save and compare evaluation performance metrics/results
library(dplyr)
library(readr)
library(tidyr)

instaprism_dir <- "instaprism_results"
dwls_dir <- "dwls_results"
cibersortx_dir <- "cibersortx/cibersortx_results"
combined_dir <- "combined_results"
dir.create(combined_dir, showWarnings = FALSE)

# Save Global Metrics
global_metrics <- bind_rows(
  read_csv(file.path(dwls_dir, "dwls_global_metrics.csv"), show_col_types = FALSE),
  read_csv(file.path(instaprism_dir, "instaprism_global_metrics.csv"), show_col_types = FALSE),
  read_csv(file.path(cibersortx_dir, "cibersortx_global_metrics.csv"), show_col_types = FALSE))
global_metrics <- global_metrics %>%
  mutate(
    Method = factor(Method, levels = c("DWLS", "InstaPrism", "CIBERSORTx")),
    Dataset = factor(Dataset, levels = c("Breast", "Lung"))) %>%
  arrange(Method, Dataset) %>%
  select(Method, Dataset, Global_Correlation, Global_RMSE)
global_metrics
write_csv(global_metrics,
          file.path(combined_dir, "combined_global_metrics.csv"))

# Save per CAF type Metrics
celltype_metrics <- bind_rows(
  read_csv(file.path(dwls_dir, "dwls_celltype_metrics.csv"), show_col_types = FALSE),
  read_csv(file.path(instaprism_dir, "instaprism_celltype_metrics.csv"), show_col_types = FALSE),
  read_csv(file.path(cibersortx_dir, "cibersortx_celltype_metrics.csv"), show_col_types = FALSE))
celltype_metrics <- celltype_metrics %>%
  mutate(
    Method = factor(Method, levels = c("DWLS", "InstaPrism", "CIBERSORTx")),
    Dataset = factor(Dataset, levels = c("Breast", "Lung"))) %>%
  arrange(Method, Dataset, CAFtype) %>%
  select(Method, Dataset, CAFtype, Correlation, RMSE)
print(celltype_metrics, n = Inf)
write_csv(celltype_metrics,
          file.path(combined_dir, "combined_celltype_metrics.csv"))
