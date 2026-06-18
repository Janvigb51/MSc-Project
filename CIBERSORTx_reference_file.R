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
# Cell names (columns) belonging to each CAF type
caf_labels <- as.character(meta$CAFtype)

# Clear any spaces if needed
caf_labels <- gsub(" ", "_", caf_labels)

# Build the first row:
# first column = GeneSymbol
# remaining columns = CAF subtype labels
first_row <- c("GeneSymbol", caf_labels)

# Build the expression table:
# first column = gene names
# remaining columns = expression values
expr_table <- cbind(
  GeneSymbol = rownames(sparse_matrix),
  as.matrix(sparse_matrix))

# Combine first row + expression table
cibersortx_ref <- rbind(first_row, expr_table)

# Write tab-delimited file without column names
write.table(
  cibersortx_ref,
  file = "CIBERSORTx_scref.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE)
