run_tsne <- function(dat) {
    dataset <- get(dat)
    if(dat != "bernstein") {
        filter1.params <- filter1_params(dataset)
        min.cells <- filter1.params$min.cells
        max.cells <- filter1.params$max.cells
        min.reads <- filter1.params$min.reads
        d <- gene_filter1(dataset, min.cells, max.cells, min.reads)
        d <- log2(1 + d)
    } else {
        d <- dataset
    }
    k <- length(unique(colnames(d)))
    if(dat == "quake_all_fpkm") {
        tsne_out <- Rtsne(t(d), perplexity = 0.5) # Run TSNE
    } else {
        tsne_out <- Rtsne(t(d)) # Run TSNE
    }
    t <- kmeans(tsne_out$Y, k, iter.max = 1e9, nstart = 1000)$clust
    return(adjustedRandIndex(t, colnames(d)))
}
