library(Seurat)

# Loading our data
breast_data <- readRDS("C:/Users/janvi/Desktop/MSc Project/scRNA-seq_dataobjects/scRNA-seq/BREAST_fibro_tumour.rds")

# Create sparse count matrix: genes x cells
sparse_matrix <- breast_data@assays$RNA@counts

# Load metadata
meta <- breast_data@meta.data

# Make sure metadata order matches count matrix columns
meta <- meta[colnames(sparse_matrix), ]
# check match
stopifnot(all(rownames(meta) == colnames(sparse_matrix)))

# Get CAF subtype labels
caf_labels <- meta$CAFtype   # change this if your column has a different name

# Rename cell columns so CIBERSORTx knows the cell type
colnames(sparse_matrix) <- paste0(caf_labels, "_", seq_along(caf_labels))

# Optional but recommended: downsample cells per CAF subtype first
# set.seed(20240618)
# 
# cells_to_keep <- unlist(
#   tapply(seq_along(caf_labels), caf_labels, function(x) {
#     sample(x, min(length(x), 200))
#   })
# )
# 
# sparse_matrix_small <- sparse_matrix[, cells_to_keep]

# Convert to a data frame for writing
cibersortx_ref <- data.frame(
  GeneSymbol = rownames(sparse_matrix),
  as.matrix(sparse_matrix),
  check.names = FALSE
)

# Write tab-delimited file
write.table(
  cibersortx_ref,
  file = "CIBERSORTx_CAF_single_cell_reference.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)