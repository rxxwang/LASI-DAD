library(dplyr)
library(psych)
library(irlba)
library(xtable)

# Read all filtered residuals file
results_list = list()
for(i in 1:22){
  results_list[[i]] <- readRDS(snakemake@input[["residuals_filter"]][i])
}
combined_result <- do.call(rbind, results_list)

# PCA
cpg_prune = colnames(combined_result)
pca_result <- prcomp_irlba(combined_result, n = 20, center = TRUE)

#Scree plot
sdev <- pca_result$sdev
var_explained <- sdev^2 / sum(sdev^2)
png(snakemake@output[["scree_plot"]])
plot(var_explained[1:20], type = "b", pch = 19,
     xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "Scree Plot")
dev.off()
saveRDS(pca_result, file = snakemake@output[["pca_result"]])

# Correlation of PCs with biomarkers
lasi = readRDS(snakemake@input[["phenos"]])
lasi$abeta_ratio = lasi$w1abeta42_final/lasi$w1abeta40_final
cor.result = corr.test(pca_result$x[,1:10], lasi[,c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")])
sig.cor = cor.result$ci %>% filter(lower * upper > 0)

# Make correlation table
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


print(xtable(formatted_df), sanitize.text.function = identity, include.rownames = TRUE,
      file = snakemake@output[["cor_table"]])

# PCA loading plot
loading = pca_result$rotation
manifest <- readRDS(snakemake@input[["manifest"]])
PCA_loading = data.frame(CpG = cpg_prune) %>% left_join(manifest) %>%
  mutate(Chr = as.numeric(gsub("chr", "", Chr)))

pdf(snakemake@output[["pca_loading"]], width = 12, height = 8)
par(mfrow = c(2, 3))
for(i in 1:6){
  PCA_loading$loading = loading[,i]
  qqman::manhattan(
    x = PCA_loading,               
    chr = "Chr",                    
    bp = "Position",                
    p = "loading",                   
    snp = "CpG",                     
    logp = FALSE,                  
    annotateTop = FALSE,           
    ylim = c(0,ceiling(max(loading[,1:6])*50)/50),
    ylab = paste0("PC", i, " loading"),
    annotatePval = 0.1,
    main = paste0("Manhattan Plot of PC", i, " loading")
  )
}
dev.off()
