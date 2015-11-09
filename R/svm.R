svm_real <- function(dataset, j, kern) {
    if(!grepl("bernstein", deparse(substitute(dataset)))) {
        filter1.params <- filter1_params(dataset)
        min.cells <- filter1.params$min.cells
        max.cells <- filter1.params$max.cells
        min.reads <- filter1.params$min.reads
        dataset <- gene_filter1(dataset, min.cells, max.cells, min.reads)
        dataset <- log2(1 + dataset)
    }
    
    k <- length(unique(colnames(dataset)))
    svm.num.cells <- dim(dataset)[2]*j
    
    teach <- sample(1:dim(dataset)[2], svm.num.cells)
    study <- setdiff(1:dim(dataset)[2], teach)
    study <- dataset[ , study]
    dataset <- dataset[, teach]
    labs1 <- colnames(dataset)
    colnames(dataset) <- NULL
    labs2 <- colnames(study)
    colnames(study) <- NULL
    
    distances <- c("euclidean", "pearson", "spearman")
    dimensionality.reductions <- c("pca", "spectral")
    
    n.cells <- dim(dataset)[2]
    n.dim <- floor(0.04 * n.cells) : ceiling(0.07 * n.cells)
    
    if(length(n.dim) > 15) {
        n.dim <- sample(n.dim, 15)
    }
    
    hash.table <- expand.grid(distan = distances,
                              dim.red = dimensionality.reductions,
                              n.dim = n.dim, stringsAsFactors = F)
    
    # register local cluster
    cl <- makeCluster(detectCores() - 1, outfile="")
    registerDoParallel(cl, cores = detectCores() - 1)
    
    dists = foreach(i = distances, .packages = "clustools") %dopar% {
        try({
            calculate_distance(dataset, i)
        })
    }
    names(dists) <- distances

    labs = foreach(i = 1:dim(hash.table)[1], .packages = "clustools",
                   .combine = rbind) %dopar% {
                       try({
                           t <- transformation(get(hash.table[i, 1], dists), hash.table[i, 2])[[1]]
                           s <- kmeans(t[, 1:hash.table[i, 3]],
                                       k,
                                       iter.max = 1e+09,
                                       nstart = 1000)$cluster
                           return(s)
                       })
                   }
    # stop local cluster
    stopCluster(cl)
    
    res <- cbind(hash.table, apply(labs, 1, paste, collapse = " "))
    colnames(res)[4] <- "labs"
    res$labs <- as.character(res$labs)
    rownames(res) <- NULL

    dat <- consensus_clustering(res$labs)
    
    diss <- dist(dat)
    clust <- hclust(diss)
    clusts <- cutree(clust, k = k)

    dataset <- t(dataset)
    rownames(dataset) <- NULL

    model <- tryCatch(svm(dataset, factor(as.character(clusts)), kernel = kern),
                      error = function(cond) return(NA))
    if(!is.na(model)) {
        pred <- predict(model, t(study))
        return(adjustedRandIndex(c(clusts, pred), c(labs1, labs2)))
    } else {
        return(NA)
    }
}

svm_real_bernstein <- function(dataset, j, kern) {
    k <- length(unique(colnames(dataset)))
    svm.num.cells <- dim(dataset)[2]*j
    
    teach <- sample(1:dim(dataset)[2], svm.num.cells)
    study <- setdiff(1:dim(dataset)[2], teach)
    study <- dataset[ , study]
    dataset <- dataset[, teach]
    labs1 <- colnames(dataset)
    colnames(dataset) <- NULL
    labs2 <- colnames(study)
    colnames(study) <- NULL
    
    distances <- c("euclidean", "pearson", "spearman")
    dimensionality.reductions <- c("pca", "spectral")
    
    n.cells <- dim(dataset)[2]
    n.dim <- floor(0.04 * n.cells) : ceiling(0.07 * n.cells)
    
    if(length(n.dim) > 15) {
        n.dim <- sample(n.dim, 15)
    }
    
    hash.table <- expand.grid(distan = distances,
                              dim.red = dimensionality.reductions,
                              n.dim = n.dim, stringsAsFactors = F)
    
    # register local cluster
    cl <- makeCluster(detectCores() - 1, outfile="")
    registerDoParallel(cl, cores = detectCores() - 1)
    
    dists = foreach(i = distances, .packages = "clustools") %dopar% {
        try({
            calculate_distance(dataset, i)
        })
    }
    names(dists) <- distances
    
    labs = foreach(i = 1:dim(hash.table)[1], .packages = "clustools",
                   .combine = rbind) %dopar% {
                       try({
                           t <- transformation(get(hash.table[i, 1], dists), hash.table[i, 2])[[1]]
                           s <- kmeans(t[, 1:hash.table[i, 3]],
                                       k,
                                       iter.max = 1e+09,
                                       nstart = 1000)$cluster
                           return(s)
                       })
                   }
    # stop local cluster
    stopCluster(cl)
    
    res <- cbind(hash.table, apply(labs, 1, paste, collapse = " "))
    colnames(res)[4] <- "labs"
    res$labs <- as.character(res$labs)
    rownames(res) <- NULL
    
    dat <- consensus_clustering(res$labs)
    
    diss <- dist(dat)
    clust <- hclust(diss)
    clusts <- cutree(clust, k = k)
    
    dataset <- t(dataset)
    rownames(dataset) <- NULL
    
    model <- tryCatch(svm(dataset, factor(as.character(clusts)), kernel = kern),
                      error = function(cond) return(NA))
    if(!is.na(model)) {
        pred <- predict(model, t(study))
        return(adjustedRandIndex(c(clusts, pred), c(labs1, labs2)))
    } else {
        return(NA)
    }
}

