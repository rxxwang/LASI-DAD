library(dplyr)
library(readxl)

# Read all residuals file
nchunks <- snakemake@params[["nchunks"]]
results_list = list()
for(i in 1:nchunks){
  results_list[[i]] <- readRDS(snakemake@input[["res_chunk"]][i])
}
combined_result <- do.call(rbind, results_list)

manifest_excludeX <- readRDS(snakemake@input[["manifest"]])

# get each residual data frame based on chr and save
for(chr in 1:22){
  manifest_chr = manifest_excludeX %>% filter(Chr == paste0('chr', chr))
  data_filtered_chr = combined_result[rownames(combined_result) %in% manifest_chr$CpG,]
  print(paste0("chr", chr, ": ",dim(data_filtered_chr)))
  saveRDS(data_filtered_chr, file = snakemake@output[["residuals"]][chr])
}

