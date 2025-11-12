library(dplyr)
library(ggplot2)

# setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/")
# model10_data = readRDS("results3/model10.w1gfap_final.RDS")
# model11_data = readRDS("results3/model11.w1gfap_final.RDS")
biomarker = snakemake@params[["biomarker"]]
model10_data <- readRDS(snakemake@input[["model10_data"]])
model10_data = model10_data[,-7]
model11_data <- readRDS(snakemake@input[["model11_data"]])
model11_data = model11_data[,-7]

all.equal(model10_data$CpG,model11_data$CpG)
model10_data[model10_data$row_varcor <= 1e-8, c("beta", "SE", "p", "plate_varcor")] <- model11_data[model10_data$row_varcor <= 1e-8, c("beta", "SE", "p", "plate_varcor")]
model10_data = model10_data %>% filter(!is.nan(p))
print(nrow(model10_data))
model10_data$p.adj = p.adjust(as.numeric(model10_data$p), method = "fdr")

saveRDS(model10_data, file = snakemake@output[["fillin_nan_result"]])

qqplots = list()
pvals <- sort(model10_data$p)  # remove NAs and sort
n <- length(pvals)
qq_data <- data.frame(expected = -log10(ppoints(n)),
                      observed = -log10(pvals))
chisq <- qchisq(1 - pvals, df = 1)
lambda <- median(chisq) / qchisq(0.5, df = 1)
qqplots = ggplot(qq_data, aes(x = expected, y = observed)) +
  geom_point(size = 1, alpha = 0.6) +
  geom_abline(
    intercept = 0,
    slope = 1,
    color = "red",
    linetype = "dashed"
  ) +
  labs(
    title = paste0("QQ Plot of ", biomarker, " P-values (p=", n,")"),
    x = expression(Expected ~  ~ -log[10](italic(p))),
    y = expression(Observed ~  ~ -log[10](italic(p)))
  ) +
  theme_minimal() +
  annotate(
    "text",
    x = max(qq_data$expected) * 0.3,
    y = max(qq_data$observed) * 0.5,
    label = paste0("lambda == ", round(lambda, 3)),
    parse = TRUE
  )
ggsave(qqplots,
       file = snakemake@output[["qqplot"]],
       width = 7,
       height = 5)


sum(model10_data$p < 5e-8)
sum(model10_data$p.adj < 0.05)

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
cpg_gr <- GRanges(seqnames = model10_data$Chr,
                  ranges = IRanges(start = model10_data$pos, end = model10_data$pos))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
nearest_idx <- nearest(cpg_gr, genes_gr)
nearest_gene_id <- genes_gr$gene_id[nearest_idx]
nearest_gene_symbol <- mapIds(org.Hs.eg.db,
                              keys = nearest_gene_id,
                              column = "SYMBOL",
                              keytype = "ENTREZID",
                              multiVals = "first")
model10_data$nearest_gene <- nearest_gene_symbol

library(knitr)
model10_data %>%
  arrange(p.adj) %>%
  mutate(beta = round(beta,2), SE= round(SE,2), 
         p = formatC(p, format = "e", digits = 2), 
         p.adj = formatC(p.adj, format = "e", digits = 2)) %>%
  dplyr::select(CpG, chr, pos, beta, SE, p, p.adj, nearest_gene) %>%
  head(20) %>%
  kable(format = "latex", booktabs = TRUE)