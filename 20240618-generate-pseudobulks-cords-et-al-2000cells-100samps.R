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
})

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("SummarizedExperiment")
BiocManager::install("BiocParallel")

## Introduction
### The purpose of this Rmd is to generate 100 pseudobulk datasets of CAF subpopulations using the `SimBu` R pacakge. Thanks to John O'Grady for suggesting this R package.

breast_data <- readRDS("C:/Users/janvi/Desktop/MSc Project/scRNA-seq_dataobjects/scRNA-seq/BREAST_fibro_tumour.rds")
# don't read in with header as R does not like duplicate colnames
# Sparse matrix
sparse_matrix <- breast_data@assays$RNA@counts
df <- data.frame(ID = rownames(meta), cell_type = meta$CAFtype)
df

cords_breast_CAFtype <- SimBu::dataset(
  annotation = df, #cell_id belonging to a cell type
  count_matrix = sparse_matrix, #genes per unique cell_ids
  #tpm_matrix = NULL,
  name = "cords_breast_CAFtype"
)

print("simulating data")
start.time <- Sys.time()
simulation_mirror_db <- SimBu::simulate_bulk(
  data = cords_breast_CAFtype,
  scenario = "mirror_db",
  ncells = 2000,
  nsamples = 100,
  scaling_factor = "read_number", # use number of reads as scaling factor
  BPPARAM = BiocParallel::MulticoreParam(workers = 3), # this will use 4 threads to run the simulation
  run_parallel = TRUE,
  balance_even_mirror_scenario = 0.01,
  seed = 20240618
)

# View/Study the Seurat Object + View the Expression Matrix
class(breast_data) # our file contains a seurat obejct
Assays(breast_data) # the object contains 2 versions of the expression data (data layers)
DefaultAssay(breast_data) # seurat will use SCT by default
GetAssayData(breast_data, assay = "SCT", layer = "counts")[1:5, 1:5] # tiny preview of the expression matrix, 5 first genes and 5 first cells

# View the Metadata from Seurat
view(breast_data@meta.data)
# Make a smaller annotation matrix
meta <- breast_data@meta.data
annot_small <- meta[, c("phenotype_final", "CAFtype"), drop = FALSE]
# check how many cells belong to each label
table(annot_small$phenotype_final)
table(annot_small$CAFtype)

# Extract the count matrix from Seurat
counts <- LayerData(breast_data, assay = "RNA", layer = "counts")
# make sure annotation order matches count matrix columns
annot_small <- annot_small[colnames(counts), ]
# check that they match
all(rownames(annot_small) == colnames(counts))

# Make genes x cells matrix
avg_by_phenotype <- sapply(
  split(seq_len(ncol(counts)), annot_small$phenotype_final),
  function(cell_ids) {
    Matrix::rowMeans(counts[, cell_ids, drop = FALSE])
  }
)
# Check the final dimensions
dim(avg_by_phenotype)
# Output: 18,984 genes and 16,704 cells

# Simulation SimBu
#1 Choose the phenotype_final annotation
celltype_col <- "phenotype_final"
#2 Remove cells with missing labels
keep <- !is.na(annot_small[[celltype_col]]) & annot_small[[celltype_col]] != ""
#3 Apply to both counts and annotation
counts_sparse <- counts[, keep]
annot_use <- annot_small[keep, , drop = FALSE]
#4 Create SimBU annotation table
annotation <- data.frame(
  ID = rownames(annot_use),
  cell_type = annot_use[[celltype_col]],
  stringsAsFactors = FALSE
)
# Final safety checks
stopifnot(ncol(counts_sparse) == nrow(annotation))
stopifnot(all(colnames(counts_sparse) == annotation$ID))
# Check labels
table(annotation$cell_type)

######################################################################################################################
# filter out some genes during simulation instead with SimBu
breast_data_counts <- read.table("~/Documents/PhD/Projects/caf-bc/data/deconvolution/cibersort_inputs/cords_ref_20240220/deconv/cords_breast_caf_data_labelled_filter_bayesprism.txt")
counts_sparse <- as.sparse(breast_data_counts[2:nrow(breast_data_counts),2:ncol(breast_data_counts)])
rownames(counts_sparse) <- breast_data_counts[2:nrow(breast_data_counts),1]
colnames(counts_sparse) <- paste0("cell_", rep(1:ncol(counts_sparse)))
cell_types <- unlist(breast_data_counts[1,2:ncol(breast_data_counts)])
names(cell_types) <- NULL
IDs <- colnames(counts_sparse)
annotation <- data.frame(
  ID = IDs,
  cell_type = cell_types,
  row.names = NULL
)
rm(breast_data_counts)
gc()
#lung_data <- readRDS("~/Documents/PhD/Projects/caf-bc/data/scRNA-seq-cords/LUNG_fibro_tumour.rds")


# simulate dataset
stopifnot(ncol(counts_sparse) == nrow(annotation))
ds <- SimBu::dataset(
  annotation = annotation, #cell_id belonging to a cell type
  count_matrix = counts_sparse, #genes per unique cell_ids
  #tpm_matrix = NULL,
  name = "cords_dataset"
)

print("simulating data")
start.time <- Sys.time()
simulation_mirror_db <- SimBu::simulate_bulk(
  data = ds,
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
write.table(bulk_counts_df, file = "~/Documents/PhD/Projects/caf-bc/data/deconvolution/cibersort_inputs/benchmark/deconv/cords-breast-100samples-2000cells-20240618/bulk_counts_2000cells_100samps.txt", quote = F, sep = "\t", row.names = F)
write.table(simulation_mirror_db$cell_fractions, file = "~/Documents/PhD/Projects/caf-bc/data/deconvolution/cibersort_inputs/benchmark/deconv/cords-breast-100samples-2000cells-20240618/bulk_props_2000cells_100samps.txt", quote = F, sep = "\t", row.names = T)
