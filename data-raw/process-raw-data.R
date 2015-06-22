# Biase, F. H., Cao, X. & Zhong, S. Cell fate inclination within 2-cell and
# 4-cell mouse embryos revealed by single-cell RNA sequencing. Genome Res. 24,
# 1787–1796 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57249
zhong <- read.table("inst/extdata/GSE57249_fpkm.txt", header = T)
zhong <- as.matrix(zhong[ , 2:50])
colnames(zhong) <- c(rep(1, 9), rep(2, 20), rep(3, 20))
save(zhong, file = "data/zhong.rda")
