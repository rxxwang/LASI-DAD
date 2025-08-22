library(dplyr)
library(psych)
library(irlba)
library(readxl)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/residuals")

args <- commandArgs(trailingOnly = TRUE)
chr <- as.numeric(args[1])

print(0)

load("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/residuals/Combined_result.Rda")

print(1)

load('/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024.rda')
HRS_Gap_Probes_v1 <- read_excel("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250528_EWAS1/HRS Gap Probes v1.xlsx")
manifest = Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024[,c(9,10,12)]
manifest$CpG = rownames(manifest)
rm(Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024)
colnames(manifest) = c('Chr', 'Position', 'Gene', 'CpG')
manifest_excludeX = manifest %>% filter(Chr != 'chrX' & Chr != 'chrY' & Chr != 'chrM') %>%
  mutate(CPG = sub("_.*", "", CpG)) %>%
  filter(!(CPG %in% HRS_Gap_Probes_v1$CPG))

manifest_chr = manifest_excludeX %>% filter(Chr == paste0('chr', chr))
data_filtered_chr = combined_result[rownames(combined_result) %in% manifest_chr$CpG,]
print(dim(data_filtered_chr))

data_filtered_chr = t(as.matrix(data_filtered_chr))
cpg = colnames(data_filtered_chr)
res_sd <- apply(data_filtered_chr, 2, sd, na.rm = TRUE)
cpg_list <- cpg[order(res_sd, decreasing = TRUE)]
keep_list <- c()

n <- dim(data_filtered_chr)[1]
cor_thresh = 0.1
while(length(cpg_list) > 0){
  next_cpg <- cpg_list[1]
  keep_list <- c(keep_list, next_cpg)
  cpg_list <- cpg_list[-1]
  cors <- cor(as.matrix(as.data.frame(data_filtered_chr)[[next_cpg]]), 
              as.matrix(as.data.frame(data_filtered_chr)[,cpg_list]))
  cpg_list <- cpg_list[abs(cors) < cor_thresh]
  print(length(cpg_list))
}
data_filtered_chr2 = as.matrix(as.data.frame(data_filtered_chr[,keep_list]))
save(data_filtered_chr2, file = paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/residuals_chr/keep_cpg_", chr, ".Rda"))

