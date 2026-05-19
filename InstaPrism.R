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

# example deconvolution

## 1) We need the pseudobulk counts (genes across samples)
view(bulk_counts_df) #or
bulk_counts_df[1:5,1:5]

## 2) Prepare the reference
OV_ref = InstaPrism_reference('OV') 

# deconvolution with InstaPrism
deconv_res = InstaPrism(bulk_Expr = bulk_counts_df,refPhi_cs = OV_ref)

# Predicted cell type proportions
estimated_frac = t(deconv_res@Post.ini.ct@theta)
head(estimated_frac)

# deconvolved gene expression for 1 cell type
Z = get_Z_array(deconv_res) # a sample by gene by cell-type array
head(Z[,1:10,'malignant'])
