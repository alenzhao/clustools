consensus_prove_concept <- function(filename, k, cell.filter, n.starts) {
    dataset <- get(filename)

    # hard cell filter
    # more than 2000 genes have to be expressed in each cell
    filt <- ""
    if(cell.filter) {
        filt <- "filt"
        
        if(filename == "quake_all_fpkm") {
            read.th <- 2000
        }

        dataset <- dataset[ , colSums(dataset > 1e-2) > read.th]
    
        if(dim(dataset)[2] == 0) {
            cat("Your dataset did not pass cell filter (more than 2000 genes have to be expressed in each cell)! Stopping now...")
            return()
        }
    }
    
    labs.known <- colnames(dataset)
    
    distances <- c("euclidean", "pearson", "spearman")
    dimensionality.reductions <- c("pca", "spectral")
    
    cat("1. Preliminary gene filtering...\n")
    filter1.params <- filter1_params(dataset)
    min.cells <- filter1.params$min.cells
    max.cells <- filter1.params$max.cells
    min.reads <- filter1.params$min.reads
    dataset <- gene_filter1(dataset, min.cells, max.cells, min.reads)
    
    if(dim(dataset)[1] == 0) {
        cat("Your dataset did not pass gene filter! Stopping now...")
        return()
    }
    
    cat("2. Log2-transforming data...\n")
    if(filename != "bernstein") {
        dataset <- log2(1 + dataset)
    }

    n.cells <- dim(dataset)[2]
    n.dim <- floor(0.05 * n.cells) : ceiling(0.08 * n.cells)
    
    if(length(n.dim) > 15) {
        n.dim <- sample(n.dim, 15)
    }
    
    hash.table <- expand.grid(distan = distances,
                              dim.red = dimensionality.reductions,
                              n.dim = n.dim, stringsAsFactors = F)
    
    # register local cluster
    cl <- makeCluster(detectCores() - 1, outfile="")
    registerDoParallel(cl, cores = detectCores() - 1)
    
    cat("3. Calculating distance matrices...\n")
    dists = foreach(i = distances, .packages = "clustools") %dopar% {
        try({
            calculate_distance(dataset, i)
        })
    }
    names(dists) <- distances
    
    pb <- txtProgressBar(min = 1, max = dim(hash.table)[1], style = 3)
    
    cat("4. Performing dimensionality reduction and kmeans clusterings...\n")
    labs = foreach(i = 1:dim(hash.table)[1], .packages = "clustools",
                   .combine = rbind) %dopar% {
                       try({
                           t <- transformation(get(hash.table[i, 1], dists), hash.table[i, 2])[[1]]
                           s <- kmeans(t[, 1:hash.table[i, 3]],
                                             k,
                                             iter.max = 1e+09,
                                             nstart = n.starts)$cluster
                           setTxtProgressBar(pb, i)
                           return(s)
                       })
                   }
    
    close(pb)
    
    res <- cbind(hash.table, apply(labs, 1, paste, collapse = " "))
    colnames(res)[4] <- "labs"
    res$labs <- as.character(res$labs)
    rownames(res) <- NULL
    
    res1 <- cbind(hash.table, apply(labs, 1, adjustedRandIndex, labs.known))
    colnames(res1)[4] <- "ARI"
    rownames(res1) <- NULL
    
    cat("5. Computing consensus matrix and labels...\n")
    all.combinations <- NULL
    for(i in 1:length(distances)) {
        for(j in 1:length(dimensionality.reductions)) {
            dist.combs <- combn(distances, i)
            dim.red.combs <- combn(dimensionality.reductions, j)
            for(m in 1:dim(dist.combs)[2]) {
                for(n in 1:dim(dim.red.combs)[2]) {
                    all.combinations <- rbind(
                        all.combinations,
                        cbind(paste(dist.combs[, m], collapse = " "),
                              paste(dim.red.combs[, n], collapse = " ")))
                }
            }
        }
    }

    cons = foreach(i = 1:dim(all.combinations)[1], .packages = "clustools", .combine = "rbind") %dopar% {
        try({
            d <- res[res$distan %in% strsplit(all.combinations[i, 1], " ")[[1]] &
                         res$dim.red %in% strsplit(all.combinations[i, 2], " ")[[1]], ]
            
            dat <- consensus_clustering(d$labs)
            
            diss <- dist(dat)
            clust <- hclust(diss)
            clusts <- cutree(clust, k = k)
            return(adjustedRandIndex(clusts, labs.known))
        })
    }
    
    
    
    # stop local cluster
    stopCluster(cl)
    res2 <- cbind(all.combinations, cons)
    rownames(res2) <- NULL
    colnames(res2) <- c("distan", "dim.red", "ARI")
    res2 <- as.data.frame(res2)
    res2[ , 3] <- as.numeric(as.character(res2[ , 3]))
    
    # how consensus over n changes the ARI
    t <- merge(res1, res2, by = c("distan", "dim.red"))
    colnames(t)[4:5] <- c("Individual", "Consensus")
    t <- as.data.table(t)
    t <- unique(melt(t[, c(1,2,4,5), with = FALSE]))
    p <- ggplot(t, aes(variable, value)) +
        geom_boxplot() +
        facet_grid(distan ~ dim.red) +
        labs(x = "", y = "ARI") +
        theme_bw()
    ggsave(paste0(filename, "-", filt, "-cons1.pdf"))
    
    t1 <- data.frame(res1[,4], rep("Individual", length(res1[,4])), stringsAsFactors = F)
    colnames(t1) <- c("ARI", "Method")
    t2 <- data.frame(res2[,3], rep("Consensus", length(res2[,3])), stringsAsFactors = F)
    colnames(t2) <- c("ARI", "Method")
    res.all <- rbind(t1, t2)
    
    write.table(res.all, paste0(filename, "-", filt, "-cons2.txt"))
    
    p <- ggplot(res.all, aes(Method, as.numeric(ARI))) +
        geom_boxplot() +
        labs(x = "", y = "ARI") +
        theme_bw()
    ggsave(paste0(filename, "-", filt, "-cons2.pdf"), w = 4, h = 3)
}
