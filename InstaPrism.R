install.packages("gitcreds")
library(gitcreds)

gitcreds::gitcreds_delete()
Sys.unsetenv("GITHUB_PAT")
Sys.unsetenv("GITHUB_TOKEN")

# install InstaPrism once
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("Biobase", quietly = TRUE))
  BiocManager::install("Biobase")

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

devtools::install_github(
  "humengying0907/InstaPrism",
  auth_token = NULL
)
library(InstaPrism)

# example deconvolution

## 1) We need the pseudobulk counts (genes across samples) into matrix with genes as row names
bulk_expr <- SummarizedExperiment::assays(simulation_mirror_db$bulk)[["bulk_counts"]]
#check
dim(bulk_expr)
bulk_expr[1:5, 1:5]
#convert to R matrix
bulk_expr <- as.matrix(bulk_expr)

## 2) Prepare the reference
refPhi_obj <- InstaPrism::refPrepare(sc_Expr = sparse_matrix, 
                       cell.type.labels = as.character(annot_table$cell_type),
                       cell.state.labels = as.character(annot_table$cell_type)) 
refPhi_obj

## 3) Deconvolution with InstaPrism
deconv_res = InstaPrism(bulk_Expr = bulk_expr, refPhi_cs = refPhi_obj)
#check output
class(deconv_res)
slotNames(deconv_res)

# Predicted cell type proportions
estimated_frac = t(deconv_res@Post.ini.ct@theta)
head(estimated_frac)

# Ground truth cell type proportions
truth_frac <- simulation_mirror_db$cell_fractions
dim(truth_frac)
head(truth_frac)

## 4) Align rows/columns 
truth_frac <- truth_frac[rownames(estimated_frac), colnames(estimated_frac)]

## 5) Performance Evaluation
## 5.1) Global Pearson Correlation
correlation <- cor(as.vector(as.matrix(truth_frac)),
                  as.vector(as.matrix(estimated_frac)))
correlation
# Cell-type-specific correlation
correlation_by_celltype <- sapply(colnames(truth_frac),
function(ct) {cor(truth_frac[, ct],
estimated_frac[, ct])})
correlation_by_celltype

## 5.2) Root Mean Squared Error
rmse <- sqrt(mean((as.matrix(truth_frac) - as.matrix(estimated_frac))^2))
rmse
# Cell-type-specific RMSE
rmse_by_celltype <- sapply(colnames(truth_frac),
function(ct) {sqrt(mean((truth_frac[, ct] - estimated_frac[, ct])^2))})
rmse_by_celltype

# deconvolved gene expression for 1 cell type
Z = get_Z_array(deconv_res) # a sample by gene by cell-type array
head(Z[,1:10,'malignant'])
