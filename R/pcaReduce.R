pca_reduce <- function(dataset) {
    if(deparse(substitute(dataset)) != "bernstein") {
        Input <- t(log2(dataset + 1)) # data matrix, cells in rows, genes in columns
    } else {
        print("bernstein")
        Input <- t(dataset) # data matrix, cells in rows, genes in columns
    }
    true_labs <- rownames(Input)
    Output_S <- PCAreduce(Input, nbt=100, q=30, method='S')
    N <- length(Output_S)
    M <- dim(Output_S[[1]])[2]
    K <- c()
    for (n in 1:N){
        cls_cell <- c()
        labels <- c()
        for (m in 1:M){
            cls_cell <- c(cls_cell, adjustedRandIndex(Output_S[[n]][,m], true_labs))
            labels <- c(labels, length(unique(Output_S[[n]][,m])))
        }
        K <- cbind(K, cls_cell)
    }
    return(K[which(labels == length(unique(rownames(Input)))), ])
}
