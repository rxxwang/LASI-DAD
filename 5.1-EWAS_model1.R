library(glmmTMB)
library(dplyr)
library(readxl)
print(1)

args <- commandArgs(trailingOnly = TRUE)
array_id <- as.numeric(args[1])
adbio_n <- as.numeric(args[2])
print(2)

model = "model1"

load(paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/all_data/all_data_",array_id,".Rda"))
print(3)
adbio <- c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")
biomarker = adbio[adbio_n]

n = length(cpg)
results = matrix(0, nrow = n, ncol = 21)
lasi_meth$abeta_ratio = lasi_meth$w1abeta42_final/lasi_meth$w1abeta40_final
print(Sys.time())

for(i in 1:n) {
  beta = cpg[i]
  temp_meth = na.omit(data.frame(lasi_meth[, c('r1hagey','ragender','batch','Plate','r1hmbmi','MedGenome_Sample_ID',
                                               'Sample_Section','smoke','CD8T','CD4T','NK','Bcell','Mono')],
                                 meth = lasi_meth[, beta], unmeth = lasi_unmeth[, beta],
                                 biomarker = lasi_meth[,biomarker])) %>%
    mutate(
      row = factor(as.numeric(substr(Sample_Section, 3, 3))),
      m_u = meth + unmeth + 1,
      meth = meth + 1
    )
  MODEL = glmmTMB(
    meth ~ log(biomarker) + log(r1hagey) + ragender + smoke + r1hmbmi + CD8T + CD4T + NK + Bcell + Mono + batch + (1 | Plate) + (1 | row) + offset(log(m_u)),
    data = temp_meth,
    family = Gamma(link = "log")
  )
  gamma_result = summary(MODEL)
  # gamma_result = summary(MODEL_biomarker)
  if(i %% 1000 == 0){
    print(paste(Sys.time(), i, beta))
  }
  results[ifelse(i%%n == 0, n, i%%n),] =
    unname(c(gamma_result$coefficients$cond[2,], gamma_result$coefficients$cond[3,],
             gamma_result$coefficients$cond[4,], gamma_result$coefficients$cond[6,],
             gamma_result$AICtab[1], MODEL$fit$convergence == 0, beta,
             manifest_excludeX[manifest_excludeX$CpG == beta,1],
             manifest_excludeX[manifest_excludeX$CpG == beta,2]))
  
}
save(results, file = paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/results_",model,"/results_", model, "_", biomarker, "_",array_id,".Rda"))
