# Biase, F. H., Cao, X. & Zhong, S. Cell fate inclination within 2-cell and
# 4-cell mouse embryos revealed by single-cell RNA sequencing. Genome Res. 24,
# 1787–1796 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57249
zhong <- read.table("inst/extdata/GSE57249_fpkm.txt", header = T)
zhong <- as.matrix(zhong[ , 2:50])
colnames(zhong) <- c(rep(1, 9), rep(2, 20), rep(3, 20))
save(zhong, file = "data/zhong.rda")

# Klein, A. M. et al. Droplet Barcoding for Single-Cell Transcriptomics Applied
# to Embryonic Stem Cells. Cell 161, 1187–1201 (2015).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65525
files <- list.files("inst/extdata")
files <- files[grepl("GSM15994", files)]

kirschner <- read.csv(paste0("inst/extdata/", files[1]))
kirschner <- as.matrix(kirschner[2:dim(kirschner)[2]])
colnames(kirschner) <- rep(1, dim(kirschner)[2])

i <- 2
for(f in files[2:length(files)]) {
    d <- read.csv(paste0("inst/extdata/", f))
    d <- as.matrix(d[2:dim(d)[2]])
    colnames(d) <- rep(i, dim(d)[2])
    i <- i + 1
    kirschner <- cbind(kirschner, d)
}

save(kirschner, file = "data/kirschner.rda")
