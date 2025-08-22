library(dplyr)
library(psych)
library(irlba)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3")

files <- list.files(path = "residuals_chr", pattern = paste0("^keep_cpg_.*\\.Rda$"))
results_list <- lapply(files, function(f) {
  load(paste0("residuals_chr/",f))               # loads the object into the environment
  residuals
  get(ls()[ls() != "f"])  # retrieve the object just loaded
})
combined_result <- do.call(cbind, results_list)

cpg_prune = colnames(combined_result)
pca_result <- prcomp_irlba(combined_result, n = 20, center = TRUE)

sdev <- pca_result$sdev
var_explained <- sdev^2 / sum(sdev^2)
pdf("PCA_ld_result.pdf")
plot(var_explained[1:20], type = "b", pch = 19,
     xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "Scree Plot")
dev.off()
save(cpg_prune, pca_result, file = "pca_result_ld.Rda")

load("pca_result_ld.Rda")
load("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/all_data/all_data_1.Rda")

lasi$abeta_ratio = lasi$w1abeta42_final/lasi$w1abeta40_final
cor.result = corr.test(pca_result$x[,1:10], lasi[,c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")])
sig.cor = cor.result$ci %>% filter(lower * upper > 0)

r_mat <- cor.result$r
p_mat <- cor.result$p
formatted <- matrix(nrow = nrow(r_mat), ncol = ncol(r_mat))
for (i in 1:nrow(r_mat)) {
  for (j in 1:ncol(r_mat)) {
    r_val <- round(r_mat[i, j], 2)
    if (!is.na(p_mat[i, j]) && p_mat[i, j] < 0.05) {
      formatted[i, j] <- paste0("\\textbf{", r_val, "}")
    } else {
      formatted[i, j] <- as.character(r_val)
    }
  }
}
rownames(formatted) <- rownames(r_mat)
colnames(formatted) <- colnames(r_mat)
formatted_df <- as.data.frame(formatted)

library(xtable)
print(xtable(formatted_df), sanitize.text.function = identity, include.rownames = TRUE)
