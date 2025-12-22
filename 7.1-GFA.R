# library(GFA)
library(ggplot2)
library(dplyr)
library(data.table)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/")
adbio = c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")

data = list()
for(i in 1:5){
  # data[[i]] = readRDS(paste0("results4_fillin_nan/model2.fillin.nan.",adbio[i],".RDS")) %>%
  #   mutate(z = beta/SE) %>% dplyr::select(CpG, z)
  data[[i]] = readRDS(snakemake@input[["data"]][i]) %>%
    mutate(z = beta/SE) %>% dplyr::select(CpG, z)
}
dt_list <- lapply(data, as.data.table)
for (i in seq_along(dt_list)) {
  setnames(dt_list[[i]], old = "z", new = adbio[i])
}
merged <- Reduce(function(x, y) merge(x, y, by = "CpG", all = FALSE), dt_list)
cpg = merged$CpG
merged = as.matrix(merged[,-1])
rownames(merged) = cpg

results_list = list()
for(i in 1:22){
  # results_list[[i]] <- readRDS(paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/residuals_chr4_rearrange/keep_cpg.", i,".RDS"))
  results_list[[i]] <- readRDS(snakemake@input[["residuals_filter"]][i])
}
combined_result <- do.call(cbind, results_list)
cpg_prune = colnames(combined_result)
merged = merged[rownames(merged) %in% cpg_prune,]
cpg = rownames(merged)
colnames(merged) = adbio

print(1)
pheno = readRDS(snakemake@input[["pheno"]])
# pheno = readRDS("data4/meth_pheno_data.1.RDS")
pheno$abeta_ratio = pheno$w1abeta42_final/pheno$w1abeta40_final
biomarker = pheno[,adbio]
N = unname(colSums(!is.na(biomarker)))
cor_matrix = cor(as.matrix(biomarker), use = "pairwise.complete.obs")

# save(merged, N, cor_matrix, file = "gfadata.Rda")
save(merged, N, cor_matrix, file = snakemake@output[["gfadata"]])

# print(2)
# gfa_model = gfa_fit(Z_hat = merged, 
#                     N = N, 
#                     R = cor_matrix)
# save(gfa_model, file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250924_Diagnostics/gfa_model.Rda")


# load("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250924_Diagnostics/gfa_model.Rda")
# data = list()
# for(i in 1:5){
#   data[[i]] = readRDS(paste0("results3_fillin_nan/model2.fillin.nan.",adbio[i],".RDS")) %>% dplyr::select(CpG, beta)
# }
# dt_list <- lapply(data, as.data.table)
# for (i in seq_along(dt_list)) {
#   setnames(dt_list[[i]], old = "beta", new = adbio[i])
# }
# merged_beta <- Reduce(function(x, y) merge(x, y, by = "CpG", all = FALSE), dt_list)
# dim(merged_beta)
# cpg = merged_beta$CpG
# merged_beta = as.matrix(merged_beta[,-1])
# rownames(merged_beta) = cpg
# colnames(merged_beta) = adbio

# data = list()
# for(i in 1:5){
#   data[[i]] = readRDS(paste0("results3_fillin_nan/model2.fillin.nan.",adbio[i],".RDS")) %>% dplyr::select(CpG, SE)
# }
# dt_list <- lapply(data, as.data.table)
# for (i in seq_along(dt_list)) {
#   setnames(dt_list[[i]], old = "SE", new = adbio[i])
# }
# merged_se <- Reduce(function(x, y) merge(x, y, by = "CpG", all = FALSE), dt_list)
# dim(merged_se)
# cpg = merged_se$CpG
# merged_se = as.matrix(merged_se[,-1])
# rownames(merged_se) = cpg
# colnames(merged_se) = adbio

# gfa_loading = gfa_loadings_gls(merged_beta, merged_se, gfa_model)
# saveRDS(gfa_loading, file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250924_Diagnostics/gfa_loading.RDS")

# factor_matrix = gfa_model$F_hat
# rownames(factor_matrix) = c("GfAP", "NfL", "p-Tau 181", "Total Tau", "Abeta Ratio")
# df <- as.data.frame(as.table(factor_matrix))
# colnames(df) <- c("Biomarker", "Factor", "value")
# ggplot(df, aes(x = Factor, y = Biomarker, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient2(
#     low = "red",
#     mid = "white",
#     high = "midnightblue",
#     midpoint = 0,
#     limits = c(-1, 1)   # optional but recommended
#   )+
#   theme_minimal(base_size = 14)
# ggsave("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250924_Diagnostics/factor_plot.png")


