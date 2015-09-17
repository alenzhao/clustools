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

quake_all_read <- read.table("inst/extdata/quake_counts.txt", header = T)
quake_labs <- read.table("inst/extdata/quake-labs-processed.txt", stringsAsFactors = F)
quake_labs <- quake_labs[, 1]
quake_all_read <- quake_all_read[ , colnames(quake_all_read) %in% quake_labs]
save(quake_all_read, file = "data/quake_all_read.rda")


# Sandberg
#
# Deng, Q., Ramsköld, D., Reinius, B. & Sandberg, R. Single-cell RNA-seq reveals
# dynamic, random monoallelic gene expression in mammalian cells. Science 343,
# 193–196 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45719
system("sh data-raw/sandberg.sh")
sandberg_all_rpkm <- read.table("inst/extdata/sandberg/sandberg-all-rpkms.txt", check.names = F, header = T)
sandberg_all_read <- read.table("inst/extdata/sandberg/sandberg-all-reads.txt", check.names = F, header = T)

genes <- sandberg_all_rpkm[ , 1]
labs <- colnames(sandberg_all_rpkm)[2:dim(sandberg_all_rpkm)[2]]

sandberg_all_rpkm <- as.matrix(sandberg_all_rpkm[ , 2:dim(sandberg_all_rpkm)[2]])
colnames(sandberg_all_rpkm) <- labs
rownames(sandberg_all_rpkm) <- genes

sandberg_all_read <- as.matrix(sandberg_all_read[ , 2:dim(sandberg_all_read)[2]])
colnames(sandberg_all_read) <- labs
rownames(sandberg_all_read) <- genes

save(sandberg_all_rpkm, file = "data/sandberg_all_rpkm.rda")
save(sandberg_all_read, file = "data/sandberg_all_read.rda")
system("rm -r inst/extdata/sandberg")
system("rm inst/extdata/GSE45719_RAW.tar")

# Pollen
#
# I got this table directly from Alex Pollen (see Gmail communication)
# it contains linear TPM values derived from HiSeq data for all genes and cells
pollen <- read.table("inst/extdata/NBT_hiseq_linear_tpm_values.txt")
labs <- colnames(pollen)
labs[grepl("Hi_2338", labs)] <- 1
labs[grepl("Hi_2339", labs)] <- 2
labs[grepl("Hi_K562", labs)] <- 3
labs[grepl("Hi_BJ", labs)] <- 4
labs[grepl("Hi_HL60", labs)] <- 5
labs[grepl("Hi_iPS", labs)] <- 6
labs[grepl("Hi_Kera", labs)] <- 7
labs[grepl("Hi_GW21.2", labs)] <- 8
labs[grepl("Hi_GW21", labs)] <- 9
labs[grepl("Hi_NPC", labs)] <- 10
labs[grepl("Hi_GW16", labs)] <- 11
labs <- as.numeric(labs)
pollen <- as.matrix(pollen)
colnames(pollen) <- labs
save(pollen, file = "data/pollen.rda")

# Usoskin
#
d <- read.csv("inst/extdata/usoskin.csv", check.names = F)
d <- as.matrix(d)
d <- d[,!grepl("Empty well", colnames(d))]
d <- d[,!grepl("NF outlier", colnames(d))]
d <- d[,!grepl("TH outlier", colnames(d))]
d <- d[,!grepl("NoN outlier", colnames(d))]
d <- d[,!grepl("NoN", colnames(d))]
d <- d[,!grepl("Central, unsolved", colnames(d))]
d <- d[,!grepl(">1 cell", colnames(d))]
d <- d[,!grepl("Medium", colnames(d))]

colnames(d)[colnames(d) == "NP1"] <- 1
colnames(d)[colnames(d) == "NP2"] <- 2
colnames(d)[colnames(d) == "NP3"] <- 3
colnames(d)[colnames(d) == "PEP1"] <- 4
colnames(d)[colnames(d) == "PEP2"] <- 5
colnames(d)[colnames(d) == "NF1"] <- 6
colnames(d)[colnames(d) == "NF2"] <- 7
colnames(d)[colnames(d) == "NF3"] <- 8
colnames(d)[colnames(d) == "NF4"] <- 9
colnames(d)[colnames(d) == "NF5"] <- 10
colnames(d)[colnames(d) == "TH"] <- 11

usoskin <- d
save(usoskin, file = "data/usoskin.rda")

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

# Bernstein
# Patel, A. P. et al. Single-cell RNA-seq highlights intratumoral heterogeneity
# in primary glioblastoma. Science 344, 1396–1401 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57872


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
#
# Klein, A. M. et al. Droplet Barcoding for Single-Cell Transcriptomics Applied
# to Embryonic Stem Cells. Cell 161, 1187–1201 (2015).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65525
system("sh data-raw/kirschner.sh")
kirschner <- read.table("inst/extdata/kirschner/reads.txt", sep = ",", row.names = 1)
kirschner <- as.matrix(kirschner)
d0 <- read.table("inst/extdata/kirschner/d0.txt")
d2 <- read.table("inst/extdata/kirschner/d2.txt")
d4 <- read.table("inst/extdata/kirschner/d4.txt")
d7 <- read.table("inst/extdata/kirschner/d7.txt")
colnames(kirschner) <- c(rep(1, d0 - 1), rep(2, d2 - 1), rep(3, d4 - 1), rep(4, d7 - 1))
save(kirschner, file = "data/kirschner.rda")
system("rm -r inst/extdata/kirschner")
system("rm inst/extdata/GSE65525.tar")
