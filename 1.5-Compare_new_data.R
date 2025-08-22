library(dplyr)
library(tidyr)
library(ggplot2)
library(gtsummary)

load("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/pheno_data.Rda")
load(paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250528_EWAS1/all_data/all_data_1.Rda"))
pheno_data_old = lasi[,9001:9064]
pheno_data_old$abeta_ratio_old = pheno_data_old$r1abeta42/pheno_data_old$r1abeta40
pheno_data_old = pheno_data_old %>% filter(!if_all(all_of(c("r1totaltau", "abeta_ratio_old", "r1nfl","r1gfap", "r1ptau181")), is.na))
pheno_data = pheno_data[,66:101]
pheno_data$abeta_ratio = pheno_data$w1abeta42_final/pheno_data$w1abeta40_final

pheno_data_old[,c("r1hagey", "ragender", "r1hmbmi", "smoke","r1totaltau", "abeta_ratio_old", "r1nfl","r1gfap", "r1ptau181")] %>% 
  mutate(smoke = factor(smoke)) %>%
  tbl_summary(
    by = ragender,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",   # for numeric variables
      all_categorical() ~ "{n} ({p}%)"      # for factor/character variables
    ),
    digits = all_continuous() ~ 2           # optional: round to 2 decimal places
  ) %>%
  as_kable(format = "latex", booktabs = TRUE, escape = FALSE) %>%
  cat()

pheno_data[,c("r1hagey", "ragender", "r1hmbmi", "smoke", "w1totaltau_final","abeta_ratio", "w1nfl_final", "w1gfap_final","w1ptau_final")] %>% 
  mutate(smoke = factor(smoke)) %>%
  tbl_summary(
    by = ragender,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",   # for numeric variables
      all_categorical() ~ "{n} ({p}%)"      # for factor/character variables
    ),
    digits = all_continuous() ~ 2           # optional: round to 2 decimal places
  ) %>%
  as_kable(format = "latex", booktabs = TRUE, escape = FALSE) %>%
  cat()


phenos = pheno_data_old[,c("MedGenome_Sample_ID","r1hagey", "ragender", "r1hmbmi", "smoke","r1totaltau", "abeta_ratio_old", "r1nfl","r1gfap", "r1ptau181")] %>% 
  full_join(pheno_data[,c("MedGenome_Sample_ID","r1hagey", "ragender", "r1hmbmi", "smoke", "w1totaltau_final","abeta_ratio", "w1nfl_final", "w1gfap_final","w1ptau_final")], by = c("MedGenome_Sample_ID"))

all(phenos$r1hagey.x[!is.na(phenos$r1hagey.x) & !is.na(phenos$r1hagey.y)] ==
      phenos$r1hagey.y[!is.na(phenos$r1hagey.x) & !is.na(phenos$r1hagey.y)])

adbio_name = data.frame(old = c("r1totaltau", "abeta_ratio_old", "r1nfl","r1gfap", "r1ptau181"),
                        new = c("w1totaltau_final","abeta_ratio", "w1nfl_final", "w1gfap_final","w1ptau_final"))
for(i in 1:5){
  stats <- phenos[,as.character(adbio_name[i,])] %>%
    summarise(
      total_non_na  = sum(!is.na(!!sym(adbio_name[i, 1])) & !is.na(!!sym(adbio_name[i, 2]))),
      num_same      = sum(!!sym(adbio_name[i, 1]) == !!sym(adbio_name[i, 2]), na.rm = TRUE),
      na1_val2      = sum(is.na(!!sym(adbio_name[i, 1])) & !is.na(!!sym(adbio_name[i, 2]))),
      val1_na2      = sum(!is.na(!!sym(adbio_name[i, 1])) & is.na(!!sym(adbio_name[i, 2])))
    )
  print(paste0(adbio_name[i,1], " Same: ", stats$num_same, "/", stats$total_non_na))
  if(stats$num_same != stats$total_non_na){
    diff = phenos[which(phenos[,as.character(adbio_name[i,1])] != phenos[,as.character(adbio_name[i,2])]),]
    print(diff[,c("MedGenome_Sample_ID", as.character(adbio_name[i,]))])
  }
  print(paste0("NA in old data but has value in new data: ", stats$na1_val2))
  print(paste0("Has value in old data but is NA in new data: ", stats$val1_na2))
}

df_long <- phenos[,c("r1totaltau", "w1totaltau_final")] %>%
  pivot_longer(cols = c(r1totaltau, w1totaltau_final), names_to = "column", values_to = "value") %>%
  mutate(Totaltau = case_when(column == "r1totaltau" ~ "old", column == "w1totaltau_final" ~ "new"))
ggplot(df_long, aes(x = value, fill = Totaltau)) +
  geom_density(alpha = 0.4) +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0(
             "Same: ", stats$num_same, "/", stats$total_non_na
           ),
           hjust = 1.1, vjust = 1.5,
           size = 5
  ) +
  theme_bw()
