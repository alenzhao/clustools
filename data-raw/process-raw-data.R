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
system("sh data-raw/bernstein.sh")
bernstein <- read.table("inst/extdata/bernstein.txt")

bernstein <- bernstein[,grepl("MGH26_", colnames(bernstein)) |
                           grepl("MGH264_", colnames(bernstein)) |
                           grepl("MGH28_", colnames(bernstein)) |
                           grepl("MGH29_", colnames(bernstein)) |
                           grepl("MGH30_", colnames(bernstein)) |
                           grepl("MGH31_", colnames(bernstein))]
# number of cells from 5 patients:
# MGH26 118
# MGH28 94
# MGH29 75
# MGH30 73
# MGH31 70
colnames(bernstein)[grepl("MGH26_", colnames(bernstein))] <- 1
colnames(bernstein)[grepl("MGH264_", colnames(bernstein))] <- 1
colnames(bernstein)[grepl("MGH28_", colnames(bernstein))] <- 2
colnames(bernstein)[grepl("MGH29_", colnames(bernstein))] <- 3
colnames(bernstein)[grepl("MGH30_", colnames(bernstein))] <- 4
colnames(bernstein)[grepl("MGH31_", colnames(bernstein))] <- 5
bernstein <- as.matrix(bernstein)
save(bernstein, file = "data/bernstein.rda")
system("rm inst/extdata/bernstein.txt")

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


# Kim
#
# Kim, D. H. et al. Single-cell transcriptome analysis reveals dynamic changes
# in lncRNA expression during reprogramming. Cell Stem Cell 16, 88–101 (2015).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE55291
system("sh data-raw/kim.sh")
kim_fpkm <- read.table("inst/extdata/kim/kim-all-fpkm.txt", check.names = F, header = T)

genes <- kim_fpkm[ , 1]
labs <- colnames(kim_fpkm)[2:dim(kim_fpkm)[2]]

kim_fpkm <- as.matrix(kim_fpkm[ , 2:dim(kim_fpkm)[2]])
colnames(kim_fpkm) <- labs
rownames(kim_fpkm) <- genes

save(kim_fpkm, file = "data/kim_fpkm.rda")
system("rm -r inst/extdata/kim")
system("rm inst/extdata/GSE55291_RAW.tar")


# Ting
#
# Ting, D. T. et al. Single-cell RNA sequencing identifies extracellular matrix
# gene expression by pancreatic circulating tumor cells. Cell Rep. 8, 1905–1918 (2014).
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51372
system("sh data-raw/ting.sh")
ting <- read.table("inst/extdata/GSE51372_readCounts.txt", check.names = F, header = T, sep = "\t", quote = "")
genes <- ting[ , 4]
labs <- colnames(ting)[7:dim(ting)[2]]

ting <- as.matrix(ting[ , 7:dim(ting)[2]])
colnames(ting) <- labs
rownames(ting) <- genes

# this was done manually!
ting <- ting[ , !grepl("TuGMP3", colnames(ting)) & !grepl("GMP1", colnames(ting)) & !grepl("GMP2", colnames(ting))]
clust1 <- c("MP4-28", "MP4-29", "MP4-31", "MP4-32", "MP6-10", "MP6-9", "MP6-17", "MP6-20", "MP6-2", "MP6-7", "MP6-3", "MP6-4", "MP6-5", "MP4-20", "MP4-22", "MP4-7", "MP4-8", "MP4-1", "MP4-6", "MP4-3", "MP4-4", "MP4-13", "MP4-14", "MP4-17")
clust2 <- c("WBC-12", "WBC-11", "WBC-9", "WBC-10", "WBC-8", "WBC-6", "WBC-7", "WBC-5", "WBC-3", "WBC-4", "WBC-1", "WBC-2")
clust3 <- c("MP2-1", "MP2-2", "MP4-24", "MP6-18", "MP6-6", "MP6-11", "MP6-15", "MP6-16", "MP7-21", "MP7-1", "MP7-3", "MP7-8", "MP7-9", "MP7-4", "MP7-7", "MP7-18", "MP7-20", "MP7-16", "MP7-12", "MP7-13", "MP7-31", "MP7-42", "MP7-30", "MP7-25", "MP7-29", "MP7-33", "MP7-34", "MP7-41", "MP7-37", "MP7-40", "MP2-36", "MP2-30", "MP2-32", "MP2-24", "MP2-26", "MP2-11", "MP2-4", "MP2-17", "MP2-18", "MP2-20", "MP2-21")
clust7 <- c("MP6-21", "MP3-2", "nb508-1", "MP6-19", "MP3-3", "MP3-5", "MP3-17", "MP3-21", "MP3-8", "MP3-15", "MP3-9")

labels <- colnames(ting)

labels[labels %in% clust1] <- 1
labels[labels %in% clust2] <- 2
labels[labels %in% clust3] <- 3
labels[labels %in% clust7] <- 7
labels[grepl("TuMP", labels)] <- 4
labels[grepl("nb508", labels)] <- 6
labels[grepl("MEF", labels)] <- 5

colnames(ting) <- labels

save(ting, file = "data/ting.rda")
system("rm inst/extdata/GSE51372_readCounts.txt")
