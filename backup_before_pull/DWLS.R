# 1) Install/load omnideconv
install.packages("Bioconductor")
install.packages("pak")
pak::pkg_install("omnideconv/omnideconv", dependencies = TRUE)

library(omnideconv)

### 1) Create 1st object raw sc_counts (genes x cells)
sc_counts <- lung_data@assays$RNA@counts
# Output: 18,984 genes and 16,704 cells

### 2) Create 2nd object cell_types (CAF subtypes)
cell_types <- as.character(meta$CAFtype)
# Clear any spaces if needed
cell_types <- gsub(" ", "_", caf_labels)

### 3) Create 3rd object bulk_tpm (genes x bulk samples, TPM-normalised)
bulk_tpm <- SummarizedExperiment::assays(simulation_mirror_db$bulk)[["bulk_counts"]]