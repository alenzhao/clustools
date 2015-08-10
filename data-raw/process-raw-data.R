# quake
#
# Treutlein, B. et al. Reconstructing lineage hierarchies of the distal lung
# epithelium using single-cell RNA-seq. Nature 509, 371–375 (2014).
#
# One files required to reproduce paper results:
# 1. File containing clusters identified by the paper's authors:
# http://www.nature.com/nature/journal/v509/n7500/extref/nature13173-s4.txt

# dowload and process original data file
system("sh data-raw/process-quake.sh")
# import matrix
quake_all_fpkm <- read.table("inst/extdata/quake.txt")
labs <- as.character(read.table("inst/extdata/quake.txt", nrows = 1))
system("rm inst/extdata/quake.txt")
# process matrix
colnames(quake_all_fpkm) <- labs
quake_all_fpkm <- as.matrix(quake_all_fpkm)
quake_all_fpkm <- quake_all_fpkm[ , 1:(dim(quake_all_fpkm)[2] - 2)]
# convert to FPKM values (log3 values in the original file)
quake_all_fpkm <- 3^quake_all_fpkm - 1
save(quake_all_fpkm, file = "data/quake_all_fpkm.rda")

# Sandberg
#
# Deng, Q., Ramsköld, D., Reinius, B. & Sandberg, R. Single-cell RNA-seq reveals
# dynamic, random monoallelic gene expression in mammalian cells. Science 343,
# 193–196 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45719
system("sh data-raw/process-sandberg.sh")
sandberg_new_rpkm <- read.table("inst/extdata/GSE45719_RAW/all-rpkms.txt", header = F)
sandberg_new_rpkm <- as.matrix(sandberg_new_rpkm)
colnames(sandberg_new_rpkm) <- sandberg_new_rpkm[1, ]
sandberg_new_rpkm <- sandberg_new_rpkm[2:dim(sandberg_new_rpkm)[1], ]
save(sandberg_new_rpkm, file = "data/sandberg_new_rpkm.rda")

sandberg_new_read <- read.table("inst/extdata/GSE45719_RAW/all-reads.txt", header = F)
sandberg_new_read <- as.matrix(sandberg_new_read)
colnames(sandberg_new_read) <- sandberg_new_read[1, ]
sandberg_new_read <- sandberg_new_read[2:dim(sandberg_new_read)[1], ]
save(sandberg_new_read, file = "data/sandberg_new_read.rda")

# Zhong
# Biase, F. H., Cao, X. & Zhong, S. Cell fate inclination within 2-cell and
# 4-cell mouse embryos revealed by single-cell RNA sequencing. Genome Res. 24,
# 1787–1796 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57249
zhong <- read.table("inst/extdata/GSE57249_fpkm.txt", header = T)
zhong <- as.matrix(zhong[ , 2:50])
colnames(zhong) <- c(rep(1, 9), rep(2, 20), rep(3, 20))
save(zhong, file = "data/zhong.rda")

# Kirschner
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

# Linnarsson
# Zeisel, A. et al. Brain structure. Cell types in the mouse cortex and
# hippocampus revealed by single-cell RNA-seq. Science 347, 1138–1142 (2015).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361

# Bernstein
# Patel, A. P. et al. Single-cell RNA-seq highlights intratumoral heterogeneity
# in primary glioblastoma. Science 344, 1396–1401 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57872
