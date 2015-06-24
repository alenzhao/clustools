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
#' labs.known <- as.numeric(colnames(quake))
#' w <- transformation(as.matrix(1 - cor(quake, method = "spearman")), "spectral")
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

#' Default parameters for the first filtering step
#' 
#' Returns default min.cells, max.cells, min.reads parameters used in
#' gene_filter1()
#' 
#' @param dataset Name of the toy dataset (either "quake", "sandberg",
#' "linnarsson" or "bernstein")
#' @return min.cells, max.cells, min.reads required to be able to filter genes
#' using filter1.
#' @examples
#' filter1_params("quake")
filter1_params <- function(dataset) {
    cat("Performing filtering1...\n")
    if (dataset == "quake") {
        min.cells <- 3
        max.cells <- 3
        min.reads <- 2
    } else if (dataset == "sandberg") {
        min.cells <- 12
        max.cells <- 12
        min.reads <- 2
    } else if (dataset == "linnarsson") {
        min.cells <- 180
        max.cells <- 180
        min.reads <- 2
    } else if (dataset == "bernstein") {
        min.cells <- 0
        max.cells <- 0
        min.reads <- 0
    } else if (dataset == "zhong") {
        min.cells <- 3
        max.cells <- 3
        min.reads <- 2
    } else if (dataset == "kirschner") {
        min.cells <- 160
        max.cells <- 160
        min.reads <- 2
    }
    return(list(min.cells = min.cells, max.cells = max.cells, min.reads = min.reads))
}

create_distance_matrix <- function(dataset, sel, distan) {
    filter1.params <- filter1_params(dataset)
    min.cells <- filter1.params$min.cells
    max.cells <- filter1.params$max.cells
    min.reads <- filter1.params$min.reads
    
    d <- gene_filter1(get(dataset), min.cells, max.cells, min.reads)
    cat("Log-trasforming data...\n")
    if (dataset != "bernstein") {
        d <- log2(1 + d)
    }
    cat("Performing filtering2...\n")
    d <- gene_filter2(d, sel)
    cat("Computing distance matrix...\n")
    dists <- calculate_distance(d, distan)
    return(dists)
}

#' Learning pipeline
#' 
#' Performs and evaluates kmeans clustering with a given combination of gene
#' filtering 2, distance metrics, transformation and number of dimensions used
#' in clustering.
#' 
#' @param dataset Name of the dataset. Either "quake", "sandberg", "bernstein" or
#' "linnarsson"
#' @param sel Selection method used by gene_filter2() function (either 
#' "none", "correlation", "variance", "variance_weight", "shannon_weight")
#' @param distan Distance metrics for calculating a distance matrix (either 
#' "pearson", "spearman", "euclidean", "manhattan" or "minkowski").
#' @param clust Distance matrix transformation method (either "pca", "spectral",
#' "spectral_reg" or "mds")
#' @param n.dim Number of dimension of the transformed distance matrix which is used
#' in kmeans clustering.
#' @return Two small files (depending on the transformation method) containing
#' a set of clustering evaluation indecies (*-inds.txt) together with known and
#' calculated clustering labels of the cells (*-labs.txt). Evaluation indecies are
#' the following:
#' \describe{
#'   \item{ari}{Adjusted Rand Index}
#'   \item{rand}{Rand Index}
#'   \item{jaccard}{Jaccard Index}
#'   \item{dunn}{Dunn Index}
#'   \item{davies_bouldin}{Davies Bouldin Index}
#'   \item{silhouette}{Silhouette Index}
#' }
#' 
#' @examples
#' machine_learning_pipeline("quake", "none", "spearman", "spectral", 4)
machine_learning_pipeline <- function(dataset, sel, distan, clust, n.dim) {
    labs.known <- as.numeric(colnames(get(dataset)))
    n.clusters <- length(unique(labs.known))
    dists <- create_distance_matrix(dataset, sel, distan)
    cat("Performing data transformation...\n")
    w <- transformation(dists, clust)
    cat("Performing kmeans clustering...\n")
    res <- check_kmeans_clustering(w[[1]], n.dim, n.clusters, labs.known)
    sink(paste0(clust, "-inds.txt"))
    cat(c(unlist(res)[1:6], dataset, sel, distan, clust, n.dim))
    cat("\n")
    sink()
    sink(paste0(clust, "-labs.txt"))
    cat(res$labs)
    cat("\n")
    cat(res$labs.known)
    cat("\n")
    sink()
}

#' Nearest neighbour pipeline
#' 
#' Performs and evaluates kmeans clustering with a given combination of gene
#' filtering 2, distance metrics, transformation, number of dimensions used in
#' clustering and number of nearest neighbours.
#' 
#' @param dataset Name of the dataset. Either "quake", "sandberg", "bernstein" or
#' "linnarsson"
#' @param sel Selection method used by gene_filter2() function (either 
#' "none", "correlation", "variance", "variance_weight", "shannon_weight")
#' @param distan Distance metrics for calculating a distance matrix (either 
#' "pearson", "spearman", "euclidean", "manhattan" or "minkowski").
#' @param clust Distance matrix transformation method (either "pca", "spectral",
#' "spectral_reg" or "mds")
#' @param n.dim Number of dimension of the transformed distance matrix which is used
#' in kmeans clustering.
#' @param nn Number of nearest neighbours used
#' @return Two small files (depending on the transformation method) containing
#' a set of clustering evaluation indecies (*-inds.txt) together with known and
#' calculated clustering labels of the cells (*-labs.txt). Evaluation indecies are
#' the following:
#' \describe{
#'   \item{ari}{Adjusted Rand Index}
#'   \item{rand}{Rand Index}
#'   \item{jaccard}{Jaccard Index}
#'   \item{dunn}{Dunn Index}
#'   \item{davies_bouldin}{Davies Bouldin Index}
#'   \item{silhouette}{Silhouette Index}
#' }
#' 
#' @examples
#' nearest_neighbour_pipeline("quake", "none", "spearman", "spectral", 4, 3)
nearest_neighbour_pipeline <- function(dataset, sel, distan, clust, n.dim, nn) {
    labs.known <- as.numeric(colnames(get(dataset)))
    n.clusters <- length(unique(labs.known))
    dists <- create_distance_matrix(dataset, sel, distan)
    graph <- exp(-dists/max(dists))
    cat("Identifying nearest neighbours...\n")
    graph.nn <- nearest_neighbor(graph, nn)
    graph.nn <- (graph.nn + t(graph.nn))/2
    cat("Performing data transformation...\n")
    w <- transformation(graph.nn, clust)
    cat("Performing kmeans clustering...\n")
    res <- check_kmeans_clustering(w[[1]], n.dim, n.clusters, labs.known)
    sink(paste0(clust, "-inds.txt"))
    cat(c(unlist(res)[1:6], dataset, sel, distan, clust, n.dim, nn))
    cat("\n")
    sink()
    sink(paste0(clust, "-labs.txt"))
    cat(res$labs)
    cat("\n")
    cat(res$labs.known)
    cat("\n")
    sink()
}

#' Support vector machine (SVM) learning on toy datasets
#' 
#' Perform SVM learning on either quake, sandberg or linnarsson data sets.
#' 
#' @param dataset Name of the toy dataset
#' @param teach.proportion Proportion of cells used for teaching, the rest of the
#' cells are used for studying
#' @param kern The kernel used in training and predicting: linear, polynomal, radial,
#' sigmoid. See ?svm from "clue" package for more information.
#' @return 
#' #' \describe{
#'   \item{pred}{Predicted cell labels for cells used in study}
#'   \item{ari}{Adjusted Rand Index - comparison of pred wiht gold standard}
#'   \item{model}{model parameters used in SVM}
#'   \item{training}{Training dataset}
#' }
#' @examples
#' res <- support_vector_machines("sandberg", 0.7)
support_vector_machines <- function(dataset, teach.proportion, kern) {
    filter1.params <- filter1_params(dataset)
    min.cells <- filter1.params$min.cells
    max.cells <- filter1.params$max.cells
    min.reads <- filter1.params$min.reads
    
    d <- get(dataset)
    dat <- gene_filter1(d, min.cells, max.cells, min.reads)
    samp <- sample(dim(dat)[2], round(teach.proportion*dim(dat)[2]))
    teach <- dat[ , samp]
    cat("Dimensions of teacher:\n")
    cat(dim(teach))
    cat("\n")
    study <- dat[ , setdiff(1:dim(dat)[2], samp)]
    cat("Dimensions of study:\n")
    cat(dim(study))
    cat("\n")
    
    teach <- t(teach)
    labs <- factor(rownames(teach))
    rownames(teach) <- NULL
    # length(unique(colnames(teach)))
    cat("Performing svm...\n")
    model <- tryCatch(svm(teach, labs, kernel = kern),
                      error = function(cond) return(NA))
    cat("Performing prediction...\n")
    if(!is.na(model)) {
        pred <- predict(model, t(study))
        cat(paste0("ARI: ",
                   adjustedRandIndex(as.numeric(names(pred)), as.numeric(pred)), "\n"))
        return(list(pred = pred,
                    ari = adjustedRandIndex(as.numeric(names(pred)), as.numeric(pred)),
                    model = model,
                    training = teach))
    } else {
        return(NA)
    }
}

#' Calclulate and plot confusions
#' 
#' Only works on the output of machine_learning_pipeline function from this
#' package.
#' 
#' Calculate and plot confusions.
#' 
#' @param dataset Name of the toy dataset.
#' @param n.dim Number of dimensions used for clustering of the toy dataset
#' @return Two files with values of confusions and plot of the confussions
#' facetted by different methods
#' @examples
#' confusion_pipeline("quake", 4)
confusion_pipeline <- function(dataset, n.dim) {
    sink(paste0(strsplit(dataset, "\\/")[[1]][3], "-", n.dim, ".txt"))
    files1 <- list.files(paste0(dataset, "/", n.dim, "/"))
    for(f1 in files1) {
        files2 <- list.files(paste0(dataset, "/", n.dim, "/", f1))
        for(f2 in files2) {
            for(clust in c("pca", "mds", "spectral", "spectral_reg")) {
                d <- readLines(paste0(dataset, "/", n.dim, "/", f1, "/", f2, "/",
                                      clust, "-labs.txt"))
                labs <- as.numeric(unlist(strsplit(d[1], " ")))
                labs.known <- as.numeric(unlist(strsplit(d[2], " ")))
                out <- confusion(labs, labs.known)
                for(i in 1:length(out)) {
                    cat(paste(f1, f2, clust, i, out[i], sep = "\t"))
                    cat("\n")
                }
            }
        }
    }
    sink()
    
    d <- read.table(paste0(strsplit(dataset, "\\/")[[1]][3], "-", n.dim, ".txt"),
                    sep = "\t")
    colnames(d) <- c("filter2", "distan", "transformation", "clust.ind", "confusion")
    
    p <- ggplot(d, aes(as.factor(clust.ind), confusion, group = transformation,
                       fill = transformation)) +
        geom_bar(stat = "identity", position = "dodge") +
        # geom_bar(position = "dodge", stat = "identity", width = 0.7) +
        facet_grid(filter2 ~ distan) +
        labs(x = "Cluster index", y = "Confusion measure") +
        theme_bw()
    
    ggsave(paste0(strsplit(dataset, "\\/")[[1]][3], "-", n.dim, ".pdf"), w = 15,
           h = 10)
}
