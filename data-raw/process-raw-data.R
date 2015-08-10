# quake
#
# Treutlein, B. et al. Reconstructing lineage hierarchies of the distal lung
# epithelium using single-cell RNA-seq. Nature 509, 371–375 (2014).
#
# One files required to reproduce paper results:
# 1. File containing clusters identified by the paper's authors:
# http://www.nature.com/nature/journal/v509/n7500/extref/nature13173-s4.txt

# dowload and process original data file
system("sh data-raw/quake.sh")
# import matrix
quake_all_fpkm <- read.table("inst/extdata/quake.txt")
labs <- as.character(read.table("inst/extdata/quake.txt", nrows = 1))
# process matrix
colnames(quake_all_fpkm) <- labs
quake_all_fpkm <- as.matrix(quake_all_fpkm)
quake_all_fpkm <- quake_all_fpkm[ , 1:(dim(quake_all_fpkm)[2] - 2)]
# convert to FPKM values (log3 values in the original file)
quake_all_fpkm <- 3^quake_all_fpkm - 1
save(quake_all_fpkm, file = "data/quake_all_fpkm.rda")
system("rm inst/extdata/quake.txt")

# Sandberg
#
# Deng, Q., Ramsköld, D., Reinius, B. & Sandberg, R. Single-cell RNA-seq reveals
# dynamic, random monoallelic gene expression in mammalian cells. Science 343,
# 193–196 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45719
system("sh data-raw/sandberg.sh")
sandberg_all_rpkm <- read.table("inst/extdata/sandberg/all-rpkms.txt", header = F)
sandberg_all_read <- read.table("inst/extdata/sandberg/all-reads.txt", header = F)

sandberg_all_rpkm <- as.matrix(sandberg_all_rpkm)
sandberg_all_read <- as.matrix(sandberg_all_read)

colnames(sandberg_all_rpkm) <- sandberg_all_rpkm[1, ]
colnames(sandberg_all_read) <- sandberg_all_read[1, ]

sandberg_all_rpkm <- sandberg_all_rpkm[2:dim(sandberg_all_rpkm)[1], ]
sandberg_all_read <- sandberg_all_read[2:dim(sandberg_all_read)[1], ]

save(sandberg_all_rpkm, file = "data/sandberg_all_rpkm.rda")
save(sandberg_all_read, file = "data/sandberg_all_read.rda")

system("rm -r inst/extdata/sandberg")
system("rm inst/extdata/GSE45719_RAW.tar")

# Linnarsson
#
# Zeisel, A. et al. Brain structure. Cell types in the mouse cortex and
# hippocampus revealed by single-cell RNA-seq. Science 347, 1138–1142 (2015).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361
system("sh data-raw/linnarsson.sh")
linnarsson <- read.table("inst/extdata/linnarsson.txt", sep = "\t", header = T)
labs <- read.table("inst/extdata/linnarsson.txt", nrows = 1, stringsAsFactors = F)
labs <- as.character(labs[2:length(labs)])
rownames(linnarsson) <- linnarsson[ , 1]
linnarsson <- linnarsson[ , 3:dim(linnarsson)[2]]
linnarsson <- as.matrix(linnarsson)
colnames(linnarsson) <- labs
save(linnarsson, file = "data/linnarsson.rda")
system("rm inst/extdata/linnarsson.txt")

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

# Bernstein
# Patel, A. P. et al. Single-cell RNA-seq highlights intratumoral heterogeneity
# in primary glioblastoma. Science 344, 1396–1401 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57872
