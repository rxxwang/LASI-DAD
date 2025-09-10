library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

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

feature = c("biomarker", "age", "gender", "smoke")
qqplots = list()
for(i in 1:4){
  variable = feature[i]
  pvals <- sort(as.numeric(combined_result[,4*i]))  # remove NAs and sort
  n <- length(pvals)
  qq_data <- data.frame(
    expected = -log10(ppoints(n)),
    observed = -log10(pvals)
  )
  chisq <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq) / qchisq(0.5, df = 1)
  qqplots[[i]] = ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_point(size = 1, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(
      title = paste0("QQ Plot of ",variable," P-values"),
      x = expression(Expected~~-log[10](italic(p))),
      y = expression(Observed~~-log[10](italic(p)))
    ) +
    theme_minimal() + 
    annotate(
      "text",
      x = max(qq_data$expected) * 0.3,
      y = max(qq_data$observed) * 0.5,
      label = paste0("lambda == ", round(lambda, 3)),
      parse = TRUE
    )
  ggsave(qqplots[[i]], file = snakemake@output[["qqplot"]][i], width = 7, height = 5)
}

# adjust p
data = data.frame(CpG = combined_result[,19], 
                  chr = as.numeric(str_sub(combined_result[,20], 4)), 
                  pos = as.numeric(combined_result[,21]),
                  beta = as.numeric(combined_result[,1]),
                  SE = as.numeric(combined_result[,2]),
                  p = as.numeric(combined_result[,4]),
                  p.adj = p.adjust(as.numeric(combined_result[,4]), method = "fdr"),
                  Chr = combined_result[,20]) %>% filter(!is.nan(p)) %>%
  mutate(CpG = sub("_.*", "", CpG)) %>%
  filter(chr != 0)
sum(data$p < 5e-8)
sum(data$p.adj < 0.05)
saveRDS(data, file = snakemake@output[["combined_result"]])

library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
cpg_gr <- GRanges(seqnames = data$Chr,
                  ranges = IRanges(start = data$pos, end = data$pos))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
nearest_idx <- nearest(cpg_gr, genes_gr)
nearest_gene_id <- genes_gr$gene_id[nearest_idx]
nearest_gene_symbol <- mapIds(org.Hs.eg.db,
                              keys = nearest_gene_id,
                              column = "SYMBOL",
                              keytype = "ENTREZID",
                              multiVals = "first")
data$nearest_gene <- nearest_gene_symbol

# library(ChIPseeker)
# peakAnno <- annotatePeak(cpg_gr, 
#                          TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
#                          tssRegion = c(-2000, 2000))
# as.data.frame(peakAnno)[, c("seqnames","start","annotation","geneId")]
# exons <- exons(txdb)
# introns <- intronsByTranscript(txdb, use.names=TRUE)
# promoters <- promoters(txdb, upstream=2000, downstream=200)
# overlaps_exon <- findOverlaps(cpg_gr, exons)
# overlaps_promoter <- findOverlaps(cpg_gr, promoters)
# overlaps_intron <- findOverlaps(cpg_gr, introns)
# 
# data$region <- "other"
# data$region[queryHits(overlaps_intron)] <- "intron"
# data$region[queryHits(overlaps_promoter)] <- "promoters"
# data$region[queryHits(overlaps_exon)] <- "exon"

library(knitr)
data %>%
  arrange(p.adj) %>%
  mutate(beta = round(beta,2), SE= round(SE,2), 
         p = formatC(p, format = "e", digits = 2), 
         p.adj = formatC(p.adj, format = "e", digits = 2)) %>%
  dplyr::select(CpG, chr, pos, beta, SE, p, p.adj, nearest_gene) %>%
  head(20) %>%
  kable(format = "latex", booktabs = TRUE)
