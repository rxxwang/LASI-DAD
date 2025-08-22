library(glmmTMB)
library(dplyr)
library(readxl)
library(haven)
print(1)
setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250528_EWAS1")
args <- commandArgs(trailingOnly = TRUE)
array_id <- as.numeric(args[1])
print(2)

adbio_col = c("w1abeta40_final", "w1abeta42_final", "w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final")
adbio_gen <- read_dta("/net/orion/skardia_lab/clubhouse/research/projects/LASI/Alzheimers_Disease/Updated LASI-DAD AD biomarker data/lasidad_w12adbio_final_genomicID.dta")
adbio_gen <- adbio_gen %>% filter(inw1_adbio == 1) %>% 
  dplyr::select(newid, w1abeta40_final, w1abeta42_final, w1gfap_final, w1nfl_final, w1ptau_final, w1totaltau_final) %>%
  filter(!is.na(w1abeta40_final) & newid != "")
adbio_meth <- read_dta("/net/orion/skardia_lab/clubhouse/research/projects/LASI/Alzheimers_Disease/Updated LASI-DAD AD biomarker data/lasidad_w12adbio_final_methylID.dta")
adbio_meth <- adbio_meth %>% filter(inw1_adbio == 1 & barcode_w1 != "") %>% 
  dplyr::select(barcode_w1, w1abeta40_final, w1abeta42_final, w1gfap_final, w1nfl_final, w1ptau_final, w1totaltau_final) %>%
  filter(!if_all(all_of(adbio_col), is.na)) %>%
  mutate(Sample_name = barcode_w1)

load('/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/LASI_betas_filtered_928074x2340_beadlt4NA_best_detP_Reps.RData')
load('/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024.rda')
HRS_Gap_Probes_v1 <- read_excel("HRS Gap Probes v1.xlsx")

# Exclude CpGs on chr X, chr Y and M
manifest = Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024[,c(9,10,12)]
manifest$CpG = rownames(manifest)
rm(Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024)
colnames(manifest) = c('Chr', 'Position', 'Gene', 'CpG')
manifest_excludeX = manifest %>% filter(Chr != 'chrX' & Chr != 'chrY' & Chr != 'chrM') %>%
  mutate(CPG = sub("_.*", "", CpG)) %>%
  filter(!(CPG %in% HRS_Gap_Probes_v1$CPG))

methx = beta2[rownames(beta2) %in% manifest_excludeX$CpG,]
rm(beta2)
methy = data.frame(sample=colnames(methx), t(methx), row.names = NULL)
rm(methx)
cpg = colnames(methy)[-1]
 
load("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250528_EWAS1/Methylation.Rda")
Meth = t(Meth[rownames(Meth) %in% cpg,])
Meth = data.frame(sample=rownames(Meth), Meth, row.names = NULL)
Unmeth = t(Unmeth[rownames(Unmeth) %in% cpg,])
Unmeth = data.frame(sample=rownames(Unmeth), Unmeth, row.names = NULL)

load(file='/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/LASI_phenos_2340.RData')
phenos = phenos[which(phenos$wave=='wave1'),]
phenos$MedGenome_Sample_ID = as.numeric(phenos$MedGenome_Sample_ID)
colnames(phenos)[3] = "sample"
phenos$r1smoken = as.factor(ifelse(phenos$r1smoken=='.m:Missing', NA, phenos$r1smoken))

phenos$smoke = as.factor(
  ifelse(phenos$r1smokev=='0.No', 0,
         ifelse(phenos$r1smokev=='1.Yes' & phenos$r1smoken=='0.No', 1,
                ifelse(phenos$r1smoken=='1.Yes', 2, NA))))

phenos = phenos %>% mutate(ragender = case_when(ragender == "1.Man" ~ 0, ragender == "2.Woman" ~ 1)) %>%
  inner_join(adbio_meth, by = c("Sample_name"))

wbc = read.csv('/net/orion/skardia_lab/clubhouse/research/projects/LASI/Methylation_Nov14_2024/Scott_QC_test/LASI_WBC.csv')
phenos = merge(phenos, wbc, by='sample')
phenos = phenos[-which(duplicated(phenos$MedGenome_Sample_ID)),]
lasi = merge(methy, phenos, by='sample')
lasi_meth = merge(Meth, phenos, by = 'sample')
lasi_unmeth = merge(Unmeth, phenos, by = 'sample')

pheno_data = lasi[,(ncol(lasi)-100):ncol(lasi)]
save(cpg, lasi, lasi_meth, lasi_unmeth, manifest_excludeX, file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/all_data.Rda")
save(pheno_data,  file = "/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/pheno_data.Rda")
