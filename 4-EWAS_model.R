library(dplyr)
library(stringr)
library(glmmTMB)

for (i in 1:6) {
  source(sprintf("4.%d-EWAS_model%d.R", i, i))
}

pca_result = readRDS(snakemake@output[["pca_result"]])
lasi = readRDS(snakemake@output[["data"]])
manifest = readRDS(snakemake@output[["manifest"]])

model = paste0("model", snakemake@params[["model"]])
biomarker = snakemake@params[["biomarker"]]

output_name = snakemake@output[["results"]]

model_fun <- get(paste0("model", model))
model_fun(model, lasi, manifest, biomarker, pca_result, output_name)