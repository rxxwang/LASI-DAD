library(glmmTMB)
library(dplyr)
library(stringr)

# lasi <- readRDS("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/data/meth_pheno_data.1.RDS")
lasi <- readRDS(snakemake@input[["data"]])
cpgs <- str_subset(colnames(lasi), "_M$") %>% str_replace("_M$", "")

n = length(cpgs)
residuals = matrix(0, nrow = n, ncol = nrow(lasi))
patients = data.frame(sample = lasi$MedGenome_Sample_ID, n = rep(0, nrow(lasi)))


for(i in 1:n) {
  
  # TODO: Understand zero methylation/unmethylation values
  
  Mname <- paste0(cpgs[i], "_M")
  Uname <- paste0(cpgs[i], "_U")
  formula <- as.formula(paste0(Mname, "~ log(r1hagey) + ragender + smoke + r1hmbmi + CD8T + CD4T + NK + Bcell + Mono + batch + (1 | Plate) + (1 | row) + offset(log(m_u))"))
  # May use UQ(Mname) instead of !!Mname to unquote
  temp_meth = filter(lasi, !is.na((!!Mname)) & !is.na((!!Uname))) %>%
    mutate(
      row = factor(as.numeric(substr(Sample_Section, 3, 3))),
      m_u = !!sym(Mname) + !!sym(Uname)
    ) %>% 
    dplyr::select(MedGenome_Sample_ID, !!Mname, !!Uname, r1hagey, ragender, smoke, r1hmbmi, CD8T, CD4T, NK, Bcell, Mono, batch, Plate, row, m_u) %>%
    filter(if_all(everything(), ~ !is.na(.)))
  
  MODEL = glmmTMB(
    formula,
    data = temp_meth,
    family = Gamma(link = "log")
  )
  gamma_result = summary(MODEL)
  
  if(MODEL$fit$convergence != 0) print(paste(i, cpgs[i], "No convergence"))
  
  residual_gamma = data.frame(res = residuals(MODEL), sample = temp_meth$MedGenome_Sample_ID)
  patients_all = patients %>% left_join(residual_gamma)
  patients_all$res = ifelse(is.na(patients_all$res), 0, patients_all$res)
  residuals[i,] = as.vector(patients_all$res)
}
residuals = as.data.frame(residuals) %>% mutate(cpg = cpgs)
saveRDS(residuals, file = snakemake@output[["out"]])

