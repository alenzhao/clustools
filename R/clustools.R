#' First filtering step
#' 
#' Filter genes that would not contribute to clustering, because they are either
#' expreseed or not expressed in almost all cells.
#' 
#' @param d Expression matrix with rows as genes and columns as cells
#' @param min.cells Minimum number of cells in which a given gene is expressed
#' @param max.cells Maximum number of cells in which a given gene is expressed
#' @param min.reads Minimum number of reads per gene per cell
#' @return Filtered expression matrix in which only genes that are expressed in
#' more than \code{min.cells} and in less than [total number of cells - \code{max.cells}].
#' @examples
#' gene_filter1(quake, 3, 3, 2)
gene_filter1 <- function(d, min.cells, max.cells, min.reads) {
    d <- d[rowSums(d > min.reads) >= min.cells & rowSums(d > 0) <= dim(d)[2] - max.cells, ]
    d <- unique(d)
    return(d)
}

#' Second filtering step
#' 
#' Filter genes that would not contribute to clustering, because statistically they
#' have a very low contribution
#' 
#' @param d Expression matrix with rows as genes and columns as cells
#' @param method Filtering method. There are four options: "correlation",
#' "variance", "variance_weight", "shannon_weight" or "none". If "correlation" is used
#' a Pearson correlation value is used with a threshold of 0.8.
#' @return Filtered expression matrix where statistically weak genes are excluded
#' @examples
#' gene_filter2(quake, "correlation")
gene_filter2 <- function(d, method) {
    if (method == "correlation") {
        rhos <- cor(t(scale(d, scale = T, center = T))) - diag(dim(d)[1])
        inds.cor <- which(rhos >= 0.8, arr.ind = T)
        inds.cor <- unique(inds.cor[, 2])
        inds <- setdiff(c(1:dim(d)[1]), inds.cor)
        d[inds, ]
    } else if (method == "variance") {
        v <- diag(var(t(d)))
        # extra coefficient found somewhere by Martin
        alpha_0 <- 1
        alpha <- alpha_0 * sqrt(log(dim(d)[1])/dim(d)[1])
        inds <- which(v > (1 + alpha) * median(v), arr.ind = T)
        d[inds, ]
    } else if (method == "variance_weight") {
        # add variance weight to genes
        v <- diag(var(t(d)))
        # extra coefficient found somewhere by Martin -- need to find the reference
        ws <- sqrt(1 + v)
        ws * d
    } else if (method == "shannon_weight") {
        # add Shannon weight to genes binarise the matrix first, then find number of bits in each row (normalised
        # by number of cells)
        p <- rowSums(d > 0)/dim(d)[2]
        # calculate Shannon entropy
        h <- -p * log2(p) - (1 - p) * log2(1 - p)
        # check this - at the moment the preferred genes get smaller weight
        ws <- 1/(2 - h)
        ws * d
    } else if (method == "none") {
        d
    }
}

calculate_distance <- function(d, method) {
    return(if (method == "spearman") {
        # there is no spearman distance method in 'proxy' package - have to define manually
        as.matrix(1 - cor(d, method = "spearman"))
        # the output is bit different from Martin's results - need to figure out
    } else if (method == "pearson") {
        as.matrix(1 - cor(d, method = "pearson"))
    } else {
        as.matrix(dist(t(d), method = method))
    })
}

calculate_correlation <- function(d, method) {
    return(cor(d, method = method))
}

transformation <- function(dists, method) {
    if (method == "pca") {
        prcomp(dists, center = TRUE, scale. = TRUE)$rotation
    } else if (method == "spectral") {
        L <- norm_laplacian(exp(-dists/max(dists)), 0)
        # here need to sort eigenvectors by their eigenvalues in increasing order!
        list(eigen(L)$vectors[, order(eigen(L)$values)], eigen(L)$values[order(eigen(L)$values)])
    } else if (method == "spectral_reg") {
        L <- norm_laplacian(exp(-dists/max(dists)), 1000)
        # here need to sort eigenvectors by their eigenvalues in increasing order!
        list(eigen(L)$vectors[, order(eigen(L)$values)], eigen(L)$values[order(eigen(L)$values)])
    } else if (method == "mds") {
        cmdscale(dists, k = length(labs.known) - 1)
    }
}

nearest_neighbor <- function(cors, k) {
    # the distance matrix has been calculated, but most of methods want to calculate the distance again
    # (usually Euclidean), so here I ignore these methods and do just simple ordering of distances and take
    # the first k neighbors the code is taken from http://stackoverflow.com/a/23449860/1365915
    n <- nrow(cors)
    for (i in 1:n) {
        cors[i, setdiff(1:n, order(cors[i, ], decreasing = T)[1:(k + 1)])] <- 0
    }
    return(cors)
}

check_kmeans_clustering <- function(w, n.dim, n.clusters, labs.known) {
    ids <- kmeans(w[, 1:n.dim], n.clusters, iter.max = 1e+09, nstart = 1000)
    
    ari <- adjustedRandIndex(ids$cluster, labs.known)
    rand <- extCriteria(ids$cluster, as.integer(labs.known), "rand")[[1]]
    jaccard <- extCriteria(ids$cluster, as.integer(labs.known), "jaccard")[[1]]
    dunn <- intCriteria(w, ids$cluster, "dunn")[[1]]
    davies_bouldin <- intCriteria(w, ids$cluster, "davies_bouldin")[[1]]
    silhouette <- intCriteria(w, ids$cluster, "silhouette")[[1]]
    
    return(list(ari = ari, rand = rand, jaccard = jaccard, dunn = dunn, davies_bouldin = davies_bouldin, silhouette = silhouette, 
        labs = ids$cluster, labs.known = as.numeric(labs.known)))
} 
