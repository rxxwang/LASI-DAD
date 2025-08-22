
load("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/all_data.Rda")
cpg_all = cpg
lasi_all = lasi
lasi_unmeth_all = lasi_unmeth
lasi_meth_all = lasi_meth
for(i in 1:100){
  if(i != 100){
    seq = ((i-1)*9000+1):(i*9000)
  }else{
    seq = ((i-1)*9000+1):(length(cpg_all))
  }
  cpg = cpg_all[seq]
  lasi = lasi_all[,c(seq+1, (length(cpg_all)+2):ncol(lasi_all))]
  lasi_meth = lasi_meth_all[,c(seq+1, (length(cpg_all)+2):ncol(lasi_meth_all))]
  lasi_unmeth = lasi_unmeth_all[,c(seq+1, (length(cpg_all)+2):ncol(lasi_unmeth_all))]
  save(cpg, lasi, lasi_meth, lasi_unmeth, manifest_excludeX, file = paste0("/net/orion/skardia_lab/clubhouse/research/projects/LASI/morrison_lab/20250813_EWAS3/all_data/all_data_",i,".Rda"))
}