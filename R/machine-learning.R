#' Compare kmeans clustering with ground truth clusters
#' 
#' Perform kmeans clustering and then calculates several external and internal
#' clustering indecies. Parameters of kmeans clustering: iter.max = 1e+09; nstart = 1000.
#' 
#' @param w Transformed distance matrix (transformed before by transformation() function)
#' @param n.dim Number of dimensions of \code{w} to use for kmeans clustering
#' @param n.clusters Expected number of clusters required by kmeans
#' @param labs.known Ground truth clustering labels (known from the experiment)
#' @return List of clustering indecies and clutering labels:
#' \describe{
#'   \item{ari}{Adjusted Rand Index}
#'   \item{rand}{Rand Index}
#'   \item{jaccard}{Jaccard Index}
#'   \item{dunn}{Dunn Index}
#'   \item{davies_bouldin}{Davies Bouldin Index}
#'   \item{silhouette}{Silhouette Index}
#'   \item{labs}{Clustering labels of the cells}
#'   \item{labs.known}{Known (from experiment) clustering labels of the cells}
#' }
#' @examples
#' check_kmeans_clustering(w, 4, 5, labs.known)
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

#' Check the first filtering step.
#' 
#' Evaluate Adjusted Rand index based on gene_filter1() parameters. NOTE that 
#' gene_filter2 is set to "none", distance is set to "spearman" and transformation 
#' is set to "spectral". These parameters were chosen because they provide the 
#' best clustering index.
#' 
#' @param d Expression matrix with rows as genes and columns as cells. Column 
#' names of d should represent the ground truth clustering indecies.
#' @param min.cells Minimum number of cells in which a given gene is expressed.
#' @param max.cells Maximum number of cells in which a given gene is expressed.
#' @param min.reads Minimum number of reads per gene per cell.
#' @param n.dim Number of dimension of the transformed distance matrix which is used
#' in kmeans clustering.
#' @return Adjusted Rand index of the clustering.
#' @examples
#' check_gene_filter1(quake, 3, 3, 2, 4)
check_gene_filter1 <- function(d, min.cells, max.cells, min.reads, n.dim) {
    labs.known <- as.numeric(colnames(d))
    n.clusters <- length(unique(labs.known))
    cat("Performing filtering1...\n")
    d <- gene_filter1(d, min.cells, max.cells, min.reads)
    cat("Log-trasforming data...\n")
    d <- log2(1 + d)
    cat("Performing filtering2...\n")
    d <- gene_filter2(d, "none")
    cat("Computing distance matrix...\n")
    dists <- calculate_distance(d, "spearman")
    cat("Performing data transformation...\n")
    w <- transformation(dists, "spectral")
    cat("Performing kmeans clustering...\n")
    res <- check_kmeans_clustering(w[[1]], n.dim, n.clusters, labs.known)
    return(res$ari)
}
