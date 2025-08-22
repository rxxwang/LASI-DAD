library(dplyr)
library(psych)
library(irlba)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3")

load("pca_result_ld.Rda")

loading = pca_result$rotation

load('/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024.rda')
manifest = Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024[,c(9,10,12)]
manifest$CpG = rownames(manifest)
rm(Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024)
colnames(manifest) = c('Chr', 'Position', 'Gene', 'CpG')


PCA_loading = data.frame(CpG = cpg_prune) %>% left_join(manifest) %>%
  mutate(Chr = as.numeric(gsub("chr", "", Chr)))

pdf("PCA_loading/PCA_loading.pdf", width = 12, height = 8)
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
    ylim = c(0,0.14),
    ylab = paste0("PC", i, " loading"),
    annotatePval = 0.1,
    main = paste0("Manhattan Plot of PC", i, " loading")
  )
}
dev.off()