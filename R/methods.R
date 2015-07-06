#' Laplacian computation
#' 
#' @param x Adjacency/distance matrix
#' @param tau Regularization term
#' @return Laplacian of the adjacency/distance matrix
#' @examples
#' L <- norm_laplacian(dists, 0)
norm_laplacian <- function(x, tau) {
    x <- x + tau * matrix(1, dim(x)[1], dim(x)[2])
    D <- diag(colSums(x))
    D1 <- D^(-0.5)
    D1[D1 == Inf] <- 0
    return(diag(dim(D)[1]) - D1 %*% x %*% D1)
} 

#' Match clustering solutions
#' 
#' You can relate this type of question to weighted bipartite graphs and subsets
#' of them. In a bipartite graph a matching is a subset of the edges so that no
#' two edges in the subset share a common vertex. It is called a minimum weighted
#' bipartite matching when the graph is a weighted bipartite graph and the sum of
#' all edges in the subset is minimal.
#' 
#' This could be represented as a distance matrix having the dimension of the
#' number of clusters where the value between two instances depicts the agreement
#' between these two partitions (one constraint for this approach is that there
#' is the same number of partitions in both clusterings). One clustering is
#' represented by columns and the other one by row or vice versa.
#' 
#' The agreement can be calculated as follows: Calculate the number of elements
#' in the intersection of the two partitions and subtract it twice from the sum
#' of the number of elements in both clusters. The notion behind this computation
#' is that if all elements are in the intersection, the value is zero and hence
#' it is very likely that these two partitions are mapped on each other. The
#' higher the value the more different are the two partitions.
#' 
#' Adapted from: \url{http://things-about-r.tumblr.com/post/36087795708/matching-clustering-solutions-using-the-hungarian}
#' 
#' To make the final match based on the created distance matrix solve_LSAP
#' function from 'clue' package is used. It implements Hungarian method of
#' matching with minimization/maximization.
#' 
#' @param clusteringA Cluster labels from clustering A
#' @param clusteringB Cluster labels from clustering B
#' @return Match between labels of clustering A and labels of clustering B
#' @examples
#' res <- minWeightBipartiteMatching(labs1, labs2)
minWeightBipartiteMatching <- function(clusteringA, clusteringB) {
    idsA <- unique(clusteringA)  # distinct cluster ids in a
    idsB <- unique(clusteringB)  # distinct cluster ids in b
    nA <- length(clusteringA)  # number of instances in a
    nB <- length(clusteringB)  # number of instances in b
    if (length(idsA) != length(idsB) || nA != nB) {
        stop("number of cluster or number of instances do not match")
    }
    
    nC <- length(idsA)
    tupel <- c(1:nA)
    
    # computing the distance matrix
    assignmentMatrix <- matrix(rep(-1, nC * nC), nrow = nC)
    for (i in 1:nC) {
        tupelClusterI <- tupel[clusteringA == i]
        solRowI <- sapply(1:nC,
                          function(i, clusterIDsB, tupelA_I) {
                              nA_I <- length(tupelA_I)  # number of elements in cluster I
                              tupelB_I <- tupel[clusterIDsB == i]
                              nB_I <- length(tupelB_I)
                              nTupelIntersect <- length(intersect(tupelA_I, tupelB_I))
                              return((nA_I - nTupelIntersect) + (nB_I - nTupelIntersect))
                          },
                          clusteringB, tupelClusterI)
        assignmentMatrix[i, ] <- solRowI
    }
    
    # optimization
    result <- solve_LSAP(assignmentMatrix, maximum = FALSE)
    attr(result, "assignmentMatrix") <- assignmentMatrix
    return(result)
}

#' Confusion calculation
#' 
#' First test labels are matched to gold standard labels by using
#' minWeightBipartiteMatching function from this package. Then when labels of
#' labs and labs.known are of the same order, confusions between them are
#' calculated. Confusion metrics introduced by myself:
#' 
#' 1 - number-of-matching-labels/total-number-of-labels
#' 
#' where number-of-matching-labels - number of labels in each cluster group of
#' labs that match to the labels of labs.known for the same cluster group; and
#' total-number-of-labels - total number of labels of labs.known from the same
#' cluster group.
#' 
#' @param labs Cluster labels from test clustering
#' @param labs.known Cluster labels from gold standard clustering
#' @return Vector containing confusions between labs and labs.known for each
#' cluster group
#' @examples
#' res <- confusion(labs, labs.known)
confusion <- function(labs, labs.known) {
    matching <- minWeightBipartiteMatching(labs, labs.known)
    labs.new <- labs
    for(i in 1:length(matching)) {
        labs.new[which(labs == i)] <- matching[i]
    }
    
    res <- unique(labs.known)
    
    for(i in unique(labs.known)) {
        res[i] <- 1 - sum(labs.new[labs.known == i] == i)/length(labs.new[labs.known == i])
    }
    return(res)
}

# consensus clustering analysis Cluster-based similarity partitioning algorithm
consensus_clustering <- function(clusts) {
    n.cells <- length(unlist(strsplit(clusts[1], " ")))
    res <- matrix(0, nrow = n.cells, ncol = n.cells)
    for (i in 1:length(clusts)) {
        t <- clusts[i]
        t <- as.numeric(unlist(strsplit(t, " ")))
        t <- as.matrix(dist(t))
        t[t != 0] <- -1
        t[t == 0] <- 1
        t[t == -1] <- 0
        res <- res + t
    }
    res <- res/i
    return(res)
} 

