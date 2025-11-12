library(dplyr)
library(stringr)

sessionInfo()
model_sub = as.numeric(snakemake@params[["model"]]) %% 10
model = paste0("model",snakemake@params[["model"]])
biomarker = snakemake@params[["biomarker"]]

nchunks <- snakemake@params[["nchunks"]]
results_list = list()
for(i in 1:nchunks){
  results_list[[i]] <- readRDS(snakemake@input[["data"]][i])
  print(i)
}
combined_result <- do.call(rbind, results_list)
print("finished")
dim(combined_result)

data = data.frame(CpG = combined_result[,19], 
                  chr = as.numeric(str_sub(combined_result[,20], 4)), 
                  pos = as.numeric(combined_result[,21]),
                  beta = as.numeric(combined_result[,1]),
                  SE = as.numeric(combined_result[,2]),
                  p = as.numeric(combined_result[,4]),
                  p.adj = p.adjust(as.numeric(combined_result[,4]), method = "fdr"),
                  plate_varcor = as.numeric(combined_result[,22]),
                  row_varcor = ifelse(rep(model_sub, nrow(combined_result)) == rep(0,nrow(combined_result)), as.numeric(combined_result[,23]), NA),
                  Chr = combined_result[,20]) 
saveRDS(data, file = snakemake@output[["combined_result"]])
