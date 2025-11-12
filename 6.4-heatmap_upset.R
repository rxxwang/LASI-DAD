library(pheatmap)
library(ComplexUpset)
library(ggplot2)
library(dplyr)
library(data.table)

thresh <- 1e-5
adbio = c("w1gfap_final", "w1nfl_final", "w1ptau_final", "w1totaltau_final", "abeta_ratio")
setwd("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250822_EWAS4/")

data = list()
data_all = list()
for(i in 1:5){
  #data[[i]] = readRDS(snakemake@input[["data"]][i])  %>% dplyr::select(CpG, p)
  data[[i]] = readRDS(paste0("results3_fillin_nan/model2.fillin.nan.",adbio[i],".RDS")) %>% dplyr::select(CpG, p)
  data_all[[i]] = readRDS(paste0("results3_fillin_nan/model2.fillin.nan.",adbio[i],".RDS"))
}
dt_list <- lapply(data, as.data.table)
for (i in seq_along(dt_list)) {
  setnames(dt_list[[i]], old = "p", new = adbio[i])
}
merged <- Reduce(function(x, y) merge(x, y, by = "CpG", all = FALSE), dt_list)
cpg = merged$CpG
merged = as.matrix(merged[,-1])
rownames(merged) = cpg
pmat_filtered <- merged[rowSums(merged < thresh, na.rm = TRUE) > 0, ]
shared_cpg = rownames(pmat_filtered)
pmat_log <- -log10(pmat_filtered)

zscore = list()
for(i in 1:5){
  #data[[i]] = readRDS(snakemake@input[["data"]][i])  %>% dplyr::select(CpG, p)
  zscore[[i]] = readRDS(paste0("results3_fillin_nan/model2.fillin.nan.",adbio[i],".RDS")) %>% 
    mutate(z = beta/SE) %>% dplyr::select(CpG, z) %>% filter(CpG %in% shared_cpg)
}
zscore_list <- lapply(zscore, as.data.table)
for (i in seq_along(dt_list)) {
  setnames(zscore_list[[i]], old = "z", new = adbio[i])
}
zscore_merged <- Reduce(function(x, y) merge(x, y, by = "CpG", all = FALSE), zscore_list)
zmat_filtered <- as.matrix(zscore_merged[,-1])
rownames(zmat_filtered) = shared_cpg

pdf("plots/heatmap_model2.pdf")
# pdf(snakemake@output[["heatmap"]])
pheatmap(
  zmat_filtered,
  scale = "none",
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE,   # hide CpG names if many
  main = paste0("z-score heatmap (CpGs with at least one p < ",thresh,")")
)
dev.off()

sig_mat <- merged < thresh
sig_mat <- sig_mat[rowSums(sig_mat) > 0, ]
sig_df <- as.data.frame(sig_mat)
colnames(sig_df) <- adbio
sig_df$CpG <- rownames(sig_mat)

pdf(snakemake@output[["upset"]])
upset(
  sig_df,
  intersect = adbio,
  name = "Significant CpGs (p < 1e-5)",
  min_size = 1
)
dev.off()

combination = list(
  V1 = c("w1gfap_final", "w1nfl_final"),
  V2 = c("w1gfap_final", "w1totaltau_final"),
  V3 = c("w1totaltau_final", "abeta_ratio"),
  V4 = c("w1gfap_final", "w1ptau_final"),
  v5 = c("w1nfl_final", "w1ptau_final"),
  V6 = c("w1ptau_final", "abeta_ratio"),
  V7 = c("w1nfl_final", "w1totaltau_final")
)
manifest = readRDS("manifest3.RDS")
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(knitr)
for(i in 1:7){
  sig_cpg = names(which(rowSums(pmat_log[,combination[[i]]] > 5) == 2))
  sig_cpg = manifest[manifest$CpG %in% sig_cpg,]
  cpg_gr <- GRanges(seqnames = sig_cpg$Chr,
                    ranges = IRanges(start = sig_cpg$Position, end = sig_cpg$Position))
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  genes_gr <- genes(txdb)
  nearest_idx <- nearest(cpg_gr, genes_gr)
  nearest_gene_id <- genes_gr$gene_id[nearest_idx]
  nearest_gene_symbol <- mapIds(org.Hs.eg.db,
                                keys = nearest_gene_id,
                                column = "SYMBOL",
                                keytype = "ENTREZID",
                                multiVals = "first")
  sig_cpg = data.frame(CpG = sig_cpg$CpG,
                       Chr = sig_cpg$Chr,
                       pos = sig_cpg$Position,
                       nearest_gene = nearest_gene_symbol,
                       p1 = pmat_filtered[rowSums(pmat_log[,combination[[i]]] > 5) == 2,combination[[i]][1]],
                       p2 = pmat_filtered[rowSums(pmat_log[,combination[[i]]] > 5) == 2,combination[[i]][2]]) %>% 
    left_join(data_all[[which(adbio == combination[[i]][1])]] %>% dplyr::select(CpG, beta, SE), by = join_by(CpG)) %>%
    left_join(data_all[[which(adbio == combination[[i]][2])]] %>% dplyr::select(CpG, beta, SE), by = join_by(CpG))
  

  sig_cpg %>%
    arrange(p1) %>%
    mutate(beta.x = round(beta.x,3), SE.x= round(SE.x,3),
           beta.y = round(beta.y,3), SE.y= round(SE.y,3),
           p1 = formatC(p1, format = "e", digits = 2), 
           p2 = formatC(p2, format = "e", digits = 2)) %>%
    dplyr::select(CpG, Chr, pos, beta.x, SE.x, p1, beta.y, SE.y, p2, nearest_gene) %>%
    kable(format = "latex", booktabs = TRUE)
}

sig_cpg = names(which(rowSums(pmat_log[,c(1,2,5)] > 5) >= 3))
sig_cpg = manifest[manifest$CpG %in% sig_cpg,]
cpg_gr <- GRanges(seqnames = sig_cpg$Chr,
                  ranges = IRanges(start = sig_cpg$Position, end = sig_cpg$Position))
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)
nearest_idx <- nearest(cpg_gr, genes_gr)
nearest_gene_id <- genes_gr$gene_id[nearest_idx]
nearest_gene_symbol <- mapIds(org.Hs.eg.db,
                              keys = nearest_gene_id,
                              column = "SYMBOL",
                              keytype = "ENTREZID",
                              multiVals = "first")
sig_cpg = data.frame(CpG = sig_cpg$CpG,
                     Chr = sig_cpg$Chr,
                     pos = sig_cpg$Position,
                     nearest_gene = nearest_gene_symbol,
                     p1 = pmat_filtered[rowSums(pmat_log[,c(1,2,5)] > 5) == 3 ,1],
                     p2 = pmat_filtered[rowSums(pmat_log[,c(1,2,5)] > 5) == 3,2],
                     p3 = pmat_filtered[rowSums(pmat_log[,c(1,2,5)] > 5) == 3,5]) %>% 
  left_join(data_all[[1]] %>% dplyr::select(CpG, beta, SE), by = join_by(CpG)) %>%
  left_join(data_all[[2]] %>% dplyr::select(CpG, beta, SE), by = join_by(CpG)) %>%
  left_join(data_all[[5]] %>% dplyr::select(CpG, beta, SE), by = join_by(CpG))


sig_cpg %>%
  arrange(p1) %>%
  mutate(beta.x = round(beta.x,3), SE.x= round(SE.x,3),
         beta.y = round(beta.y,3), SE.y= round(SE.y,3),
         beta.z = round(beta.y,3), SE.z= round(SE.y,3),
         p1 = formatC(p1, format = "e", digits = 2), 
         p2 = formatC(p2, format = "e", digits = 2),
         p3 = formatC(p3, format = "e", digits = 2)) %>%
  dplyr::select(CpG, Chr, pos, beta.x, SE.x, p1, beta.y, SE.y, p2, beta.z, SE.z, p3, nearest_gene) %>%
  kable(format = "latex", booktabs = TRUE)

