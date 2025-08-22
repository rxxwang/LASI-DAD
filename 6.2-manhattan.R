library(readr)
library(ggplot2)
library(stringr)
library(patchwork)

setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3")
model = "model2"
adbio_n = 1
adbio <- c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")
biomarker = adbio[adbio_n]
load(paste0("results_", model, "_", biomarker, ".Rda"))

threshold = 5e-8
chr_lengths <- data %>%
  group_by(chr) %>%
  summarise(chr_len = max(pos)) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0))  # where each chr starts

data = data %>% left_join(chr_lengths, by = "chr") %>%
  mutate(pos_cum = pos + chr_start) %>% filter(chr > 0)
axis_ticks <- data %>%
  group_by(chr) %>%
  summarize(center = mean(pos_cum))
odd_ticks <- axis_ticks %>% filter(chr %% 2 == 1)
even_ticks <- axis_ticks %>% filter(chr %% 2 == 0)
man = data %>%
  mutate(log10pval = -log10(p),
         color = case_when(p > threshold & chr %% 2 == 1 ~ 1, 
                           p > threshold & chr %% 2 == 0 ~ 2,
                           p <= threshold ~ 3)) %>%
  ggplot(aes(
    x = pos_cum,
    y = log10pval,
    color = as.factor(color)
  )) +
  geom_point(size = 0.1) + 
  geom_hline(yintercept = -log10(threshold), color = "red", linetype = "dashed", size = 0.2) +
  scale_color_manual(values = c("lightgrey", "darkgrey", "blue")) +theme_bw()+
  scale_x_continuous(labels = axis_ticks$chr, breaks = axis_ticks$center)+
  labs(y = "-log10(p-values)", x = "", title = paste0("Manhattan Plot of Model 2 ", biomarker)) +
  scale_y_continuous(limits = c(0, 20))+#facet_grid(type ~ method) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.ticks.x = element_blank(), # Remove x-axis ticks
    plot.margin = margin(t = 10, r = 10, b = -20, l = 10)  # Add space below the plot
  )

chr_labels <- ggplot() +
  geom_text(
    data = odd_ticks,
    aes(x = center, y = 0.2, label = chr),
    size = 3
  ) +
  geom_text(
    data = even_ticks,
    aes(x = center, y = 0.1, label = chr),
    size = 3
  ) +
  scale_x_continuous(
    limits = c(0,max(data$pos_cum))
  ) +
  scale_y_continuous(limits = c(0.05, 0.25))+
  theme_void() +  # Completely empty background
  theme(
    plot.margin = margin(t = -10, r = 10, b = 5, l = 10)  # Adjust margin for alignment
  )

manhattan <- man / chr_labels + plot_layout(heights = c(8, 1))
ggsave(
  plot = manhattan,
  filename = paste0("manhattan_",biomarker,".png"),
  width = 10, height = 6
)