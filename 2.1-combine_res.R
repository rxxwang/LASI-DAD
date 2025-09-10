library(dplyr)
library(readxl)

# Read all residuals file
nchunks <- snakemake@params[["nchunks"]]
results_list = list()
for(i in 1:nchunks){
  # results_list[[i]] <- readRDS(paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/null_model/null_residuals.",i,".RDS"))
  results_list[[i]] <- readRDS(snakemake@input[["res_chunk"]][i])
  print(i)
}
combined_result <- do.call(rbind, results_list)
print("finished")

# manifest_excludeX <- readRDS("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/manifest.RDS")
manifest_excludeX <- readRDS(snakemake@input[["manifest"]])

# get each residual data frame based on chr and save
for(chr in 1:22){
  manifest_chr = manifest_excludeX %>% filter(Chr == paste0('chr', chr))
  data_filtered_chr = combined_result[combined_result$cpg %in% manifest_chr$CpG,]
  print(paste0("chr", chr, ": ",dim(data_filtered_chr)[1]))
  saveRDS(data_filtered_chr, file = snakemake@output[["residuals"]][chr])
  # saveRDS(data_filtered_chr, file = paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/null_model/null_residuals.chr",i,".RDS"))
}

