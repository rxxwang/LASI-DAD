library(dplyr)
library(readxl)
library(tidyr)
library(haven)
library(lme4)

adbio_col = c("w1abeta40_final", "w1abeta42_final", "w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final")

adbio_meth <- read_dta(snakemake@input[["biomarker"]])
adbio_meth <- adbio_meth %>% filter(inw1_adbio == 1 & barcode_w1 != "") %>% 
  dplyr::select(barcode_w1, w1abeta40_final, w1abeta42_final, w1gfap_final, w1nfl_final, w1ptau_final, w1totaltau_final,
                w1abeta40_plate, w1abeta42_plate, w1gfap_plate, w1nfl_plate, w1ptau_plate, w1totaltau_plate) %>%
  filter(!if_all(all_of(adbio_col), is.na)) %>%
  mutate(Sample_name = barcode_w1)

# loads a matrix called beta2
# load(snakemake@input[["beta"]])

#loads Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024 which is a dataframe
load(snakemake@input[["manifest"]])
HRS_Gap_Probes_v1 <- read_excel(snakemake@input[["gap"]])

# Exclude CpGs on chr X, chr Y and M
manifest = Modified_EPICv2_manifest_in_TJ.Peters_BMCgenomics_2024[,c(9,10,12)] 
manifest$CpG = rownames(manifest)
colnames(manifest) = c('Chr', 'Position', 'Gene', 'CpG')
manifest_excludeX = manifest %>% filter(Chr != 'chrX' & Chr != 'chrY' & Chr != 'chrM') %>%
  mutate(CPG = sub("_.*", "", CpG)) %>%
  filter(!(CPG %in% HRS_Gap_Probes_v1$CPG))

print(-1)

# loads matrix Meth and Unmeth
load(snakemake@input[["meth"]])
load(snakemake@input[["unmeth"]])
cpgs <- manifest_excludeX$CpG
keep_meth_ix <- which(rownames(Methx) %in% cpgs)
keep_cpgs <- rownames(Methx)[keep_meth_ix]
samples = colnames(Methx)

stopifnot(all(rownames(Methx) == rownames(Unmethx)))
stopifnot(all(colnames(Methx) == colnames(Unmethx)))

print(0)

# TODO: Check that betas derived from Meth and Unmeth match beta2 from Scott
Meth = t(Methx[keep_meth_ix,])
Meth = data.frame(sample=samples, Meth, row.names = NULL)
colnames(Meth) <- c("sample", paste0(keep_cpgs, "_M"))

Unmeth = t(Unmethx[keep_meth_ix,])
Unmeth = data.frame(sample=samples, Unmeth, row.names = NULL)
colnames(Unmeth) <- c("sample", paste0(keep_cpgs, "_U"))

#loads data frame phenos
problem_id = 330565
load(file=snakemake@input[["pheno"]])
phenos = phenos[which(phenos$wave=='wave1'),]
phenos$MedGenome_Sample_ID = as.numeric(phenos$MedGenome_Sample_ID)
phenos = phenos %>% filter(MedGenome_Sample_ID != problem_id)
colnames(phenos)[3] = "sample"
phenos$r1smoken = as.factor(ifelse(phenos$r1smoken=='.m:Missing', NA, phenos$r1smoken))

phenos$smoke = as.factor(
  ifelse(phenos$r1smokev=='0.No', 0,
         ifelse(phenos$r1smokev=='1.Yes' & phenos$r1smoken=='0.No', 1,
                ifelse(phenos$r1smoken=='1.Yes', 2, NA))))

phenos = phenos %>% mutate(ragender = case_when(ragender == "1.Man" ~ 0, ragender == "2.Woman" ~ 1)) %>%
  inner_join(adbio_meth, by = c("Sample_name"))



wbc = read.csv(snakemake@input[["wbc"]])

phenos = inner_join(phenos, wbc, by='sample')

phenos <- phenos %>%
  inner_join(Meth, by = "sample") %>%
  inner_join(Unmeth, by = "sample")

# TODO: Find out what is going on with duplicates? Are they the same data twice or different data?
# phenos = phenos[-which(duplicated(phenos$MedGenome_Sample_ID)),]

# calculate the residuals of lmm for biomarkers and plate, add them to the model
col_biomarker <-data.frame(plate = c("w1abeta40_plate", "w1abeta42_plate", "w1gfap_plate", "w1nfl_plate", "w1ptau_plate", "w1totaltau_plate"),
                           biomarker = c("w1abeta40_final", "w1abeta42_final", "w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final"))
for(i in 1:6){
  data = data.frame(biomarker = phenos[,col_biomarker[i,2]],
                    plate = phenos[,col_biomarker[i,1]],
                    rownum = 1:nrow(phenos))
  data = data[complete.cases(data),]
  colnames(data) = c("biomarker", "plate", "rownum")
  data$plate = ifelse(data$plate < 2, 1, data$plate)
  data$plate = factor(data$plate)
  model <- lmer(log(biomarker) ~ (1 | plate), data)
  res <- data.frame(res = scale(residuals(model)), rownum = data$rownum) %>%
    complete(rownum = 1:max(rownum))
  phenos[[paste0(col_biomarker[i,2], "_res")]] = res$res
}

nchunks <- snakemake@params[["nchunks"]]
cpgs_per_chunk <- ceiling(length(keep_cpgs)/nchunks)
chunk <- rep(1:nchunks, each = cpgs_per_chunk)[1:length(keep_cpgs)]

pheno_data = phenos %>%
  select( -matches("(_M|_U)$") )
print(3)

for(i in 1:nchunks){
  cpgs_chunk <- keep_cpgs[chunk == i]
  data_chunk <- select(phenos, all_of(c("sample", paste0(cpgs_chunk, "_M"), paste0(cpgs_chunk, "_U")))) %>%
    inner_join(pheno_data, by = "sample")
  print(dim(data_chunk))
  saveRDS(data_chunk, file = snakemake@output[["merged_data"]][i])
}

saveRDS(keep_cpgs, file = snakemake@output[["cpg_list"]])
saveRDS(manifest_excludeX, file = snakemake@output[["manifest"]])
