library(dplyr)
library(stringr)
library(glmmTMB)

sessionInfo()
cat("Running on node:", Sys.info()[["nodename"]], "\n")
for (i in c(10, 11, 20, 21)) {
  source(sprintf("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/LASI-DAD/4.%d-EWAS_model%d.R", i, i))
}

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/")
pca_result = readRDS(snakemake@input[["pca_result"]])
lasi = readRDS(snakemake@input[["data"]])
# lasi = readRDS(paste0("data4/meth_pheno_data.",chunk,".RDS"))
manifest = readRDS(snakemake@input[["manifest"]])
lasi$abeta_ratio_res = lasi$w1abeta42_final_res/lasi$w1abeta40_final_res

model = paste0("model",snakemake@params[["model"]])
biomarker = snakemake@params[["biomarker"]]

output_name = snakemake@output[["results"]]
# output_name = paste0("results4_", model, "/results.", model, ".",biomarker,".", chunk, ".RDS")

model_fun <- get(model)
model_fun(model, lasi, manifest, biomarker, pca_result, output_name)
