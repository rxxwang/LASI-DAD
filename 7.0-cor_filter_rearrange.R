library(dplyr)
library(data.table)
library(psych)
library(irlba)

sessionInfo()

#setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/")

adbio = c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")
# data = list()
# for(i in 1:5){
#   #data[[i]] = readRDS(paste0("results3_fillin_nan/model2.fillin.nan.",adbio[i],".RDS")) %>% dplyr::select(CpG, p)
#   data[[i]] = readRDS(snakemake@input[["data"]][i]) %>% dplyr::select(CpG, p)
# }
# dt_list <- lapply(data, as.data.table)
# for (i in seq_along(dt_list)) {
#   setnames(dt_list[[i]], old = "p", new = adbio[i])
# }
# merged <- Reduce(function(x, y) merge(x, y, by = "CpG", all = FALSE), dt_list)
# cpg = merged$CpG
# merged = as.matrix(merged[,-1])
# rownames(merged) = cpg
# saveRDS(merged, file = "pval.model2.RDS")
merged = readRDS(snakemake@input[["pval"]]) 
mins <- apply(merged, 1, min)
print("aaa")
mins_sorted <- sort(mins)
print("aaaa")
cpg_sorted = names(mins_sorted)

print(1)

manifest_excludeX <- readRDS(snakemake@input[["manifest"]])
chr = snakemake@params[["chr"]]
manifest_chr = manifest_excludeX %>% filter(Chr == paste0('chr', chr))
cpg_sorted_chr = cpg_sorted[cpg_sorted %in% manifest_chr$CpG]

print(2)

data_filtered_chr <- readRDS(snakemake@input[["residuals"]])
data_filtered_chr$ord <- match(data_filtered_chr$cpg, cpg_sorted_chr)
data_filtered_chr <- data_filtered_chr[order(is.na(data_filtered_chr$ord), data_filtered_chr$ord), ]
data_filtered_chr$ord <- NULL
cpg = data_filtered_chr$cpg
data_filtered_chr = t(as.matrix(data_filtered_chr %>% dplyr::select(-cpg)))
colnames(data_filtered_chr) = cpg
cpg_list <- cpg
keep_list <- c()

print(3)

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

