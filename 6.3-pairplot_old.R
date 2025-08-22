library(dplyr)
library(ggplot2)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
adbio_n <- as.numeric(args[1])
model = paste0("model", as.numeric(args[2]))

adbio_old <- c("r1gfap","r1nfl","r1ptau181","r1totaltau", "abeta_ratio")
if(model == "model2"){
  biomarker = adbio_old[adbio_n]
  setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250702_EWAS2/")
  files <- list.files(pattern = paste0("^results_PCA_ld_",biomarker,".*\\.Rda$"))
  length(files)
  results_list <- lapply(files, function(f) {
    load(f)               # loads the object into the environment
    get(ls()[ls() != "f"])  # retrieve the object just loaded
  })
  combined_result1 <- do.call(rbind, results_list)
}

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3")

adbio <- c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")
biomarker = adbio[adbio_n]
files <- list.files(path = paste0("results_",model),
                    pattern = paste0("^results_", model, "_", biomarker, ".*\\.Rda$"),
                    full.names = TRUE)
length(files)
results_list <- lapply(files, function(f) {
  load(f)
  get(ls()[ls() != "f"])
})
combined_result2 <- do.call(rbind, results_list)

pvals <- data.frame(oldp = as.numeric(combined_result1[,4]), CpG = combined_result1[,19]) %>% 
  inner_join(data.frame(newp = as.numeric(combined_result2[,4]), CpG = combined_result2[,19]))
scatterplot = ggplot(pvals, aes(x = -log10(oldp), y = -log10(newp))) + geom_point(alpha = 0.5) + theme_bw()
ggsave(scatterplot, file = paste0("pairplot/scatterplot_", biomarker,"_old.png"))

