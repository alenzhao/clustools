#' Find structure of the data (number of clusters) based on Markov stability
#' method (http://github.com/michaelschaub/PartitionStability)
#' 
#' WARNING: Matlab has to be installed on your computer.
#' 
#' Diffusion times are defined in the interval [10^1:10^5].
#' 
#' @param dataset Name of the dataset. Either "quake", "sandberg", "bernstein" or
#' "linnarsson"
#' @param sel Selection method used by gene_filter2() function (either 
#' "none", "correlation", "variance", "variance_weight", "shannon_weight")
#' @param distan Distance metrics for calculating a distance matrix (either 
#' "pearson", "spearman", "euclidean", "manhattan" or "minkowski").
#' @param nn Number of nearest neighbours used
#' @param matlab Path to matlab source file (e.g. "/Applications/MATLAB_R2015a.app/bin/matlab")
#' @return File dataset-sel-distan-nn-N.csv containing number of clusters at
#' different diffusion times. File script.m contains a sequence of commands used
#' in Matlab to run Markov stability algorithm.
#' 
#' @examples
#' markov_stability("quake", "none", "spearman", 3, "/Applications/MATLAB_R2015a.app/bin/matlab")
markov_stability <- function(dataset, sel, distan, nn, matlab) {
    d <- get(dataset)
    dists <- create_distance_matrix(dataset, d, sel, distan)
    graph <- exp(-dists/max(dists))
    cat("Identifying nearest neighbours...\n")
    if(nn != nrow(graph)) {
        graph.nn <- nearest_neighbor(graph, nn)
        graph.nn <- (graph.nn + t(graph.nn))/2
    }
    write.table(graph.nn, file = paste0(dataset, "-", nn, "-graph.csv"),
                sep = ",", row.names = F, col.names = F)

    matlab.lines <- c(
        "T = 10.^[1:0.1:5];",
        paste0("A = csvread('", dataset, "-", nn, "-graph.csv", "');"),
        "[S, N, VI, C] = stability(A,T,'v');",
        paste0("csvwrite('", dataset, "-", sel, "-", distan, "-", nn, "-N.csv",
               "', N);"))
    writeLines(matlab.lines, con = "script.m")
    system(paste0(matlab, " -nodisplay -r \"run('script.m'); exit\""))
    
    system(paste0("rm ", dataset, "-", nn, "-graph.csv"))
    # system(paste0("rm script.m"))
}
