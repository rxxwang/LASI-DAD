library(glmmTMB)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
array_id <- as.numeric(args[1])

load(paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/all_data/all_data_",array_id,".Rda"))
n = length(cpg)
residuals = matrix(0, nrow = n, ncol = nrow(lasi))
patients = data.frame(sample = lasi$MedGenome_Sample_ID, n = rep(0, nrow(lasi)))
print(Sys.time())

for(i in 1:n) {
  beta = cpg[i]
  temp_meth = na.omit(data.frame(lasi_meth[, c('r1hagey','ragender','batch','Plate','r1hmbmi','MedGenome_Sample_ID',
                                               'Sample_Section','smoke','CD8T','CD4T','NK','Bcell','Mono')],
                                 meth = lasi_meth[, beta], unmeth = lasi_unmeth[, beta])) %>%
    mutate(
      row = factor(as.numeric(substr(Sample_Section, 3, 3))),
      meth = meth + 1,
      m_u = meth + unmeth + 1,
      meth = meth
    )
  MODEL = glmmTMB(
    meth ~ log(r1hagey) + ragender + smoke + r1hmbmi + CD8T + CD4T + NK + Bcell + Mono + batch + (1 | Plate) + (1 | row) + offset(log(m_u)),
    data = temp_meth,
    family = Gamma(link = "log")
  )
  gamma_result = summary(MODEL)
  if(i %% 1000 == 0){
    print(paste(Sys.time(), i, beta))
  }
  residual_gamma = data.frame(res = residuals(MODEL), sample = temp_meth$MedGenome_Sample_ID)
  patients_all = patients %>% left_join(residual_gamma)
  patients_all$res = ifelse(is.na(patients_all$res), 0, patients_all$res)
  residuals[ifelse(i%%n == 0, n, i%%n),] = as.vector(patients_all$res)
}
rownames(residuals) = cpg
save(residuals, file = paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/residuals/residuals_nonbiomarkers_", array_id, ".Rda"))

