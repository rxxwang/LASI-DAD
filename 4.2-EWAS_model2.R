library(dplyr)
library(stringr)
library(glmmTMB)

model2 <- function(model, lasi, manifest, biomarker, pca_result, output_name){
  # model = paste0("model", snakemake@params[["model"]])
  # 
  # lasi = readRDS(snakemake@output[["data"]])
  # manifest = readRDS(snakemake@output[["manifest"]])
  # pca_result = readRDS(snakemake@output[["pca_result"]])
  # biomarker = snakemake@params[["biomarker"]]
  
  cpgs <- str_subset(colnames(lasi), "_M$") %>% str_replace("_M$", "")
  n = length(cpgs)
  results = matrix(0, nrow = n, ncol = 21)
  lasi$abeta_ratio = lasi$w1abeta42_final/lasi$w1abeta40_final
  
  print(Sys.time())
  
  for(i in 1:n) {
    Mname <- paste0(cpgs[i], "_M")
    Uname <- paste0(cpgs[i], "_U")
    formula <- as.formula(paste0(Mname, "~ log(", biomarker, ") + log(r1hagey) + ragender + smoke + r1hmbmi + CD8T + CD4T + NK + Bcell + Mono + PC1 + PC2 + PC3 + PC4 + batch + (1 | Plate) + (1 | row) + offset(log(m_u))"))
    temp_meth = filter(lasi, !is.na((!!Mname)) & !is.na((!!Uname))) %>%
      mutate(
        row = factor(as.numeric(substr(Sample_Section, 3, 3))),
        m_u = !!sym(Mname) + !!sym(Uname),
        PC1 = scale(pca_result$x[,1]), PC2 = scale(pca_result$x[,2]), PC3 = scale(pca_result$x[,3]), PC4 = scale(pca_result$x[,4])
      )
    MODEL = glmmTMB(
      formula,
      data = temp_meth,
      family = Gamma(link = "log")
    )
    gamma_result = summary(MODEL)
    if(is.nan(gamma_result$coefficients$cond[2,4])) print(paste(i, cpgs[i], "NaN"))
    if(i %% 1000 == 0){
      print(paste(Sys.time(), i, cpgs[i]))
    }
    results[ifelse(i%%n == 0, n, i%%n),] =
      unname(c(gamma_result$coefficients$cond[2,], gamma_result$coefficients$cond[3,],
               gamma_result$coefficients$cond[4,], gamma_result$coefficients$cond[6,],
               gamma_result$AICtab[1], MODEL$fit$convergence == 0, cpgs[i],
               manifest[manifest$CpG == cpgs[i],1],
               manifest[manifest$CpG == cpgs[i],2]))
  }
  saveRDS(results, file =  output_name)
}

