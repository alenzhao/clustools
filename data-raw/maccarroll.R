# first, download the following files from:
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63473
# number of rows in each file
# 20479 GSM1626793_P14Retina_1.digital_expression.txt
# 20620 GSM1626794_P14Retina_2.digital_expression.txt
# 20170 GSM1626795_P14Retina_3.digital_expression.txt
# 20240 GSM1626796_P14Retina_4.digital_expression.txt
# 19720 GSM1626797_P14Retina_5.digital_expression.txt
# 20648 GSM1626798_P14Retina_6.digital_expression.txt
# 20106 GSM1626799_P14Retina_7.digital_expression.txt

tab5rows <-
    read.table("inst/extdata/dropseq/GSM1626793_P14Retina_1.digital_expression.txt",
               header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
d1 <- read.table("inst/extdata/dropseq/GSM1626793_P14Retina_1.digital_expression.txt",
                header = TRUE, colClasses = classes, nrow = 20479)
labs1 <- d1[,1]

tab5rows <-
    read.table("inst/extdata/dropseq/GSM1626794_P14Retina_2.digital_expression.txt",
               header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
d2 <- read.table("inst/extdata/dropseq/GSM1626794_P14Retina_2.digital_expression.txt",
                 header = TRUE, colClasses = classes, nrow = 20620)
labs2 <- d2[,1]

tab5rows <-
    read.table("inst/extdata/dropseq/GSM1626795_P14Retina_3.digital_expression.txt",
               header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
d3 <- read.table("inst/extdata/dropseq/GSM1626795_P14Retina_3.digital_expression.txt",
                 header = TRUE, colClasses = classes, nrow = 20170)
labs3 <- d3[,1]

tab5rows <-
    read.table("inst/extdata/dropseq/GSM1626796_P14Retina_4.digital_expression.txt",
               header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
d4 <- read.table("inst/extdata/dropseq/GSM1626796_P14Retina_4.digital_expression.txt",
                 header = TRUE, colClasses = classes, nrow = 20240)
labs4 <- d4[,1]

tab5rows <-
    read.table("inst/extdata/dropseq/GSM1626797_P14Retina_5.digital_expression.txt",
               header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
d5 <- read.table("inst/extdata/dropseq/GSM1626797_P14Retina_5.digital_expression.txt",
                 header = TRUE, colClasses = classes, nrow = 19720)
labs5 <- d5[,1]

tab5rows <-
    read.table("inst/extdata/dropseq/GSM1626798_P14Retina_6.digital_expression.txt",
               header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
d6 <- read.table("inst/extdata/dropseq/GSM1626798_P14Retina_6.digital_expression.txt",
                 header = TRUE, colClasses = classes, nrow = 20648)
labs6 <- d6[,1]

tab5rows <-
    read.table("inst/extdata/dropseq/GSM1626799_P14Retina_7.digital_expression.txt",
               header = TRUE, nrows = 5)
classes <- sapply(tab5rows, class)
d7 <- read.table("inst/extdata/dropseq/GSM1626799_P14Retina_7.digital_expression.txt",
                 header = TRUE, colClasses = classes, nrow = 20106)
labs7 <- d7[,1]

# number of genes in each file is different - find gene common to all files
labs.all <- labs1[labs1 %in% labs2 &
                      labs1 %in% labs3 &
                      labs1 %in% labs4 &
                      labs1 %in% labs5 &
                      labs1 %in% labs6 &
                      labs1 %in% labs7]

# select only genes from labs.all
d1 <- d1[d1[ , 1] %in% labs.all, ]
d2 <- d2[d2[ , 1] %in% labs.all, ]
d3 <- d3[d3[ , 1] %in% labs.all, ]
d4 <- d4[d4[ , 1] %in% labs.all, ]
d5 <- d5[d5[ , 1] %in% labs.all, ]
d6 <- d6[d6[ , 1] %in% labs.all, ]
d7 <- d7[d7[ , 1] %in% labs.all, ]

# make matrices
d1 <- d1[ , 2:dim(d1)[2]]
rownames(d1) <- labs.all
d1 <- as.matrix(d1)

d2 <- d2[ , 2:dim(d2)[2]]
rownames(d2) <- labs.all
d2 <- as.matrix(d2)

d3 <- d3[ , 2:dim(d3)[2]]
rownames(d3) <- labs.all
d3 <- as.matrix(d3)

d4 <- d4[ , 2:dim(d4)[2]]
rownames(d4) <- labs.all
d4 <- as.matrix(d4)

d5 <- d5[ , 2:dim(d5)[2]]
rownames(d5) <- labs.all
d5 <- as.matrix(d5)

d6 <- d6[ , 2:dim(d6)[2]]
rownames(d6) <- labs.all
d6 <- as.matrix(d6)

d7 <- d7[ , 2:dim(d7)[2]]
rownames(d7) <- labs.all
d7 <- as.matrix(d7)

cols1 <- colnames(d1)
cols1 <- paste("r1", cols1, sep = "_")
cols2 <- colnames(d2)
cols2 <- paste("r2", cols2, sep = "_")
cols3 <- colnames(d3)
cols3 <- paste("r3", cols3, sep = "_")
cols4 <- colnames(d4)
cols4 <- paste("r4", cols4, sep = "_")
cols5 <- colnames(d5)
cols5 <- paste("r5", cols5, sep = "_")
cols6 <- colnames(d6)
cols6 <- paste("r6", cols6, sep = "_")
cols7 <- colnames(d7)
cols7 <- paste("p1", cols7, sep = "_")

colnames(d1) <- cols1
colnames(d2) <- cols2
colnames(d3) <- cols3
colnames(d4) <- cols4
colnames(d5) <- cols5
colnames(d6) <- cols6
colnames(d7) <- cols7

# merge data files
t <- cbind(d1, d2)
t <- cbind(t, d3)
t <- cbind(t, d4)
t <- cbind(t, d5)
t <- cbind(t, d6)
t <- cbind(t, d7)

maccarroll <- t

# cell labels data file is downloaded from:
# http://mccarrolllab.com/dropseq/
cells <- read.table("inst/extdata/dropseq/retina_clusteridentities.txt", stringsAsFactors = F)

cols <- colnames(maccarroll)
cols <- cols[cols %in% cells[,1]]

maccarroll <- maccarroll[ , colnames(maccarroll) %in% cols]

cells <- cells[order(cells[,1]),]
cols <- colnames(maccarroll)
maccarroll <- maccarroll[ , order(cols)]
colnames(maccarroll) <- cells[ , 2]

# save one big table
save(maccarroll, file = "data/maccarroll.rda")
