library(dplyr)
library(psych)
library(ggplot2)
library(irlba)
library(glmmTMB)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3")

args <- commandArgs(trailingOnly = TRUE)
array_id <- as.numeric(args[1])
adbio_n <- as.numeric(args[2])


model = "model4"

load(paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/all_data/all_data_",array_id,".Rda"))
load("pca_result_ld.Rda")
print(3)
adbio <- c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")
biomarker = adbio[adbio_n]

n = length(cpg)
results = matrix(0, nrow = n, ncol = 21)
lasi_meth$abeta_ratio = lasi_meth$w1abeta42_final/lasi_meth$w1abeta40_final
print(Sys.time())



lasi_meth[, biomarker] <- log(lasi_meth[, biomarker])  # use log1p(df) if df has no zero values
lasi_meth[, biomarker] <- scale(lasi_meth[, biomarker])
lasi_meth[, biomarker] = ifelse(lasi_meth[, biomarker] > 3, 3, ifelse(lasi_meth[, biomarker] < -3, -3, lasi_meth[, biomarker]))

for(i in 1:n) {
  beta = cpg[i]
  temp_meth = na.omit(data.frame(lasi_meth[, c('r1hagey','ragender','batch','Plate','r1hmbmi','MedGenome_Sample_ID',
                                               'Sample_Section','smoke','CD8T','CD4T','NK','Bcell','Mono')],
                                 meth = lasi_meth[, beta], unmeth = lasi_unmeth[, beta],
                                 biomarker = lasi_meth[,biomarker],
                                 PC1 = scale(pca_result$x[,1]), PC2 = scale(pca_result$x[,2]), PC3 = scale(pca_result$x[,3]), PC4 = scale(pca_result$x[,4]))) %>%
    mutate(
      row = factor(as.numeric(substr(Sample_Section, 3, 3))),
      m_u = meth + unmeth + 1,
      #biomarker_log = log(biomarker+1),
      meth = meth + 1
      #age_scale = scale(r1hagey),
      
    )
  # temp = na.omit(data.frame(lasi[, c(beta,'r1hagey','ragender','batch','Plate','Sample_Section','r1smokev')],
  #                           beta = lasi[,beta],
  #                           biomarker = lasi[,biomarker])) %>% filter(r1smokev != ".m:Missing") %>%
  #   mutate(
  #     row = factor(as.numeric(substr(Sample_Section, 3, 3))),
  #     biomarker_log = log(biomarker+1)
  #   )
  # MODEL_biomarker = glmmTMB(
  #   log(beta) ~ biomarker_log + log(r1hagey) + ragender + r1smokev + batch + (1 | Plate) + (1 | row),
  #   data = temp,
  #   family = gaussian()
  # )
  # MODEL = glmmTMB(
  #   meth ~ log(r1hagey) + ragender + smoke + r1hmbmi + CD8T + CD4T + NK + Bcell + Mono + batch + (1 | Plate) + (1 | row) + offset(log(m_u)),
  #   data = temp_meth,
  #   family = Gamma(link = "log")
  # )
  MODEL_biomarker = glmmTMB(
    meth ~ biomarker + log(r1hagey) + ragender + smoke + r1hmbmi + CD8T + CD4T + NK + Bcell + Mono + PC1 + PC2 + PC3 + PC4 + batch + (1 | Plate) + (1 | row) + offset(log(m_u)),
    data = temp_meth,
    family = Gamma(link = "log")
  )
  # lmer_result = summary(MODEL_biomarker)
  # gamma_result = summary(MODEL)
  gamma_result = summary(MODEL_biomarker)
  if(is.nan(gamma_result$coefficients$cond[2,4])) print(paste(i, beta, "NaN"))
  if(i %% 100 == 0){
    print(paste(Sys.time(), i, beta))
  }
  # residual_gamma = data.frame(res = residuals(MODEL), sample = temp_meth$MedGenome_Sample_ID)
  # patients_all = patients %>% left_join(residual_gamma)
  # patients_all$res = ifelse(is.na(patients_all$res), 0, patients_all$res)
  # residuals[ifelse(i%%n == 0, n, i%%n),] = as.vector(patients_all$res)
  results[i,] =
    unname(c(gamma_result$coefficients$cond[2,], gamma_result$coefficients$cond[3,],
             gamma_result$coefficients$cond[4,], gamma_result$coefficients$cond[6,],
             gamma_result$AICtab[1], MODEL_biomarker$fit$convergence == 0, beta,
             manifest_excludeX[manifest_excludeX$CpG == beta,1],
             manifest_excludeX[manifest_excludeX$CpG == beta,2]))
}
save(results, file = paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/results_",model,"/results_", model, "_", biomarker, "_",array_id,".Rda"))
