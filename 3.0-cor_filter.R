library(dplyr)
library(psych)
library(irlba)

data_filtered_chr <- readRDS(snakemake@input[["residuals"]])
data_filtered_chr = t(as.matrix(data_filtered_chr))
cpg = colnames(data_filtered_chr)
res_sd <- apply(data_filtered_chr, 2, sd, na.rm = TRUE)
cpg_list <- cpg[order(res_sd, decreasing = TRUE)]
keep_list <- c()

n <- dim(data_filtered_chr)[1]
cor_thresh = snakemake@params[["cor_threshold"]]
while(length(cpg_list) > 0){
  next_cpg <- cpg_list[1]
  keep_list <- c(keep_list, next_cpg)
  cpg_list <- cpg_list[-1]
  cors <- cor(as.matrix(as.data.frame(data_filtered_chr)[[next_cpg]]), 
              as.matrix(as.data.frame(data_filtered_chr)[,cpg_list]))
  cpg_list <- cpg_list[abs(cors) < cor_thresh]
  print(length(cpg_list))
}
data_filtered_chr = as.matrix(as.data.frame(data_filtered_chr[,keep_list]))
saveRDS(data_filtered_chr, file = snakemake@output[["residuals_filter"]])

