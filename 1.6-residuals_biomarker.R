library(haven)
library(ggplot2)
library(patchwork)
library(lme4)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/")

biomarker_data = read_dta("/net/orion/skardia_lab/clubhouse/research/projects/LASI/Alzheimers_Disease/Updated LASI-DAD AD biomarker data/lasidad_w12adbio_final_methylID.dta") %>% filter(inw1_adbio == 1)
meth_pheno_data <- readRDS("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/data3/meth_pheno_data.1.RDS")
biomarker_data = biomarker_data[biomarker_data$barcode_w1 %in% meth_pheno_data$Sample_name,]
col_biomarker <-data.frame(plate = c("w1abeta40_plate", "w1abeta42_plate", "w1gfap_plate", "w1nfl_plate", "w1ptau_plate", "w1totaltau_plate"),
                           biomarker = c("w1abeta40_final", "w1abeta42_final", "w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final"))
plot = list()
for(i in 1:6){
  data = data.frame(biomarker = biomarker_data[,col_biomarker[i,2]],
                    plate = biomarker_data[,col_biomarker[i,1]])
  data = data[complete.cases(data),]
  colnames(data) = c("biomarker", "plate")
  data$plate = ifelse(data$plate < 2, 1, data$plate)
  data$plate = factor(data$plate)
  plot[[i]] <- ggplot(data, aes(x = plate, y = biomarker)) +
    geom_boxplot() + ylab(col_biomarker[i,2])
  model <- lmer(biomarker ~ (1 | plate), data)
  res <- data.frame(res = residuals(model), rownum = data$rownum) %>%
    complete(rownum = 1:max(rownum))
}
combined <- plot[[1]] + plot[[2]] + plot[[3]] + 
  plot[[4]] + plot[[5]] + plot[[6]] +
  plot_layout(ncol = 3)  # arrange in 2 rows, 3 columns
ggsave("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250924_Diagnostics/boxplot_biomarker_plate.png", combined, width = 12, height = 8)

