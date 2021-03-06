% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/machine-learning.R
\name{check_kmeans_clustering}
\alias{check_kmeans_clustering}
\title{Compare kmeans clustering with ground truth clusters}
\usage{
check_kmeans_clustering(w, n.dim, n.clusters, labs.known)
}
\arguments{
\item{w}{Transformed distance matrix (transformed before by transformation() function)}

\item{n.dim}{Number of dimensions of \code{w} to use for kmeans clustering}

\item{n.clusters}{Expected number of clusters required by kmeans}

\item{labs.known}{Ground truth clustering labels (known from the experiment)}
}
\value{
List of clustering indecies and clutering labels:
\describe{
  \item{ari}{Adjusted Rand Index}
  \item{rand}{Rand Index}
  \item{jaccard}{Jaccard Index}
  \item{dunn}{Dunn Index}
  \item{davies_bouldin}{Davies Bouldin Index}
  \item{silhouette}{Silhouette Index}
  \item{labs}{Clustering labels of the cells}
  \item{labs.known}{Known (from experiment) clustering labels of the cells}
}
}
\description{
Perform kmeans clustering and then calculates several external and internal
clustering indecies. Parameters of kmeans clustering: iter.max = 1e+09; nstart = 1000.
}
\examples{
labs.known <- as.numeric(colnames(quake))
w <- transformation(as.matrix(1 - cor(quake, method = "spearman")), "spectral")
check_kmeans_clustering(w, 4, 5, labs.known)
}

