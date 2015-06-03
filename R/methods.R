norm_laplacian <- function(x, tau) {
    # Assume that x is the adjacency matrix and that tau is a regularization term. Default is to have tau = 0.
    x <- x + tau * matrix(1, dim(x)[1], dim(x)[2])
    D <- diag(colSums(x))
    D1 <- D^(-0.5)
    D1[D1 == Inf] <- 0
    return(diag(dim(D)[1]) - D1 %*% x %*% D1)
} 
