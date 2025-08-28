library(glmmTMB)
library(dplyr)

model = paste0("model", snakemake@params[["model"]])

lasi = readRDS(snakemake@output[["data"]])
manifest = readRDS(snakemake@output[["manifest"]])
biomarker = snakemake@params[["biomarker"]]

cpgs <- str_subset(colnames(lasi, "_M$")) %>% str_replace("_M$", "")
n = length(cpgs)
results = matrix(0, nrow = n, ncol = 21)
lasi$abeta_ratio = lasi$w1abeta42_final/lasi$w1abeta40_final

invnorm = function(x){
  x_order <- order(x)
  p <- seq(0, 1, length.out = length(x) + 2)[ 2: (length(x) + 1)]
  normal_quantiles <- qnorm(p = p)
  x_inv_normalized <- normal_quantiles[x_order]
  return(x_inv_normalized)
}
lasi_meth[,biomarker] = invnorm(lasi_meth[,biomarker])

for(i in 1:n) {
  Mname <- paste0(cpgs[i], "_M")
  Uname <- paste0(cpgs[i], "_U")
  formula <- as.formula(paste0(Mname, "~ log(", biomarker, ") + log(r1hagey) + ragender + smoke + r1hmbmi + CD8T + CD4T + NK + Bcell + Mono + batch + (1 | Plate) + (1 | row) + offset(log(m_u))"))
  
  temp_meth = filter(lasi, !is.na((!!Mname)) & !is.na((!!Uname))) %>%
    mutate(
      row = factor(as.numeric(substr(Sample_Section, 3, 3))),
      m_u = !!Mname + !!Uname
    )
  
  MODEL = glmmTMB(
    formula,
    data = temp_meth,
    family = Gamma(link = "log")
  )
  gamma_result = summary(MODEL)
  if(i %% 1000 == 0){
    print(paste(Sys.time(), i, cpg[i]))
  }
  results[ifelse(i%%n == 0, n, i%%n),] =
    unname(c(gamma_result$coefficients$cond[2,], gamma_result$coefficients$cond[3,],
             gamma_result$coefficients$cond[4,], gamma_result$coefficients$cond[6,],
             gamma_result$AICtab[1], MODEL$fit$convergence == 0, cpgs[i],
             manifest[manifest$CpG == cpgs[i],1],
             manifest[manifest$CpG == cpgs[i],2]))
}
saveRDS(results, file = snakemake@output[["results"]])
