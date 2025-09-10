library(dplyr)
library(stringr)
library(glmmTMB)

for (i in 1:6) {
  source(sprintf("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/LASI-DAD/4.%d-EWAS_model%d.R", i, i))
}

pca_result = readRDS(snakemake@input[["pca_result"]])
lasi = readRDS(snakemake@input[["data"]])
manifest = readRDS(snakemake@input[["manifest"]])

model = paste0("model",snakemake@params[["model"]])
biomarker = snakemake@params[["biomarker"]]

output_name = snakemake@output[["results"]]

model_fun <- get(model)
model_fun(model, lasi, manifest, biomarker, pca_result, output_name)
