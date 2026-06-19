# Load Packages & Libraries
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("SimBu")
BiocManager::install("SummarizedExperiment")
BiocManager::install("BiocParallel")

install.packages("Seurat")
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(Matrix)
  library(harmony)
  library(SimBu)
  library(tibble)
  library(tidyr)
  library(SummarizedExperiment)
})

## Introduction
### The purpose of this Rmd is to generate 100 pseudobulk datasets of CAF subpopulations using the `SimBu` R pacakge.
### Thanks to John O'Grady for suggesting this R package.
lung_data <- readRDS("C:/Users/janvi/Desktop/MSc Project/scRNA-seq_dataobjects/scRNA-seq/LUNG_fibro_tumour.rds")
# don't read in with header as R does not like duplicate colnames
# load in metadata
meta <- lung_data@meta.data

### 1) Create the Sparse Count Matrix (genes x cells)
sparse_matrix <- lung_data@assays$RNA@counts
# Output: 18,984 genes and 16,704 cells

# Make sure metadata order matches count matrix columns
meta <- meta[colnames(sparse_matrix), ]
# check match
stopifnot(all(rownames(meta) == colnames(sparse_matrix)))
# Create new simple cell IDs
new_cell_ids <- paste0("cell", seq_len(ncol(sparse_matrix)))

# Save old-to-new mapping in case we need it later
cell_id_map <- data.frame(
  old_id = colnames(sparse_matrix),
  new_id = new_cell_ids
)

# Rename columns of count matrix
colnames(sparse_matrix) <- new_cell_ids
# Rename rows of metadata
rownames(meta) <- new_cell_ids

### 2) Create the Annotation Table (cells x info)
annot_table <- data.frame(ID = rownames(meta), cell_type = meta$cluster_ft)
view(annot_table)
# check match
stopifnot(all(annot_table$ID == colnames(sparse_matrix)))

### 3) Simulation
cords_lung <- SimBu::dataset(
  annotation = annot_table, #cell_id belonging to a cell type
  count_matrix = sparse_matrix, #genes per unique cell_ids
  #tpm_matrix = NULL,
  name = "CAFtype"
)

print("simulating data")
start.time <- Sys.time()
simulation_mirror_db <- SimBu::simulate_bulk(
  data = cords_lung,
  scenario = "mirror_db",
  ncells = 2000,
  nsamples = 100,
  scaling_factor = "read_number", # use number of reads as scaling factor
  BPPARAM = BiocParallel::MulticoreParam(workers = 3), # this will use 4 threads to run the simulation
  run_parallel = TRUE,
  balance_even_mirror_scenario = 0.01,
  seed = 20240618
)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

bulk_counts_df <- as.data.frame(SummarizedExperiment::assays(simulation_mirror_db$bulk)[["bulk_counts"]]) %>% rownames_to_column(var = "Gene")
print("writing to file")
write.table(bulk_counts_df, file = "C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-lung-100samples-2000cells-20240618/bulk_counts_2000cells_100samps.txt", quote = F, sep = "\t", row.names = F)
write.table(simulation_mirror_db$cell_fractions, file = "C:/Users/janvi/Desktop/MSc Project/deconv_inputs/cords-lung-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt", quote = F, sep = "\t", row.names = T)
# simulation_mirror_db$bulk is our pseudo bulk samples that we feed in Insta Prism
# simulation_mirror_db$cell_fractions is our ground truth data that we use to compare deconv results (the estimated)
