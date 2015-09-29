find_outliers <- function(dataset, clusts, cell.filt = F) {
    if(cell.filt) {
        dataset <- dataset[ , colSums(dataset > 1e-2) > 2000]
    }
    
    # gene filter
    filter1.params <- filter1_params(dataset)
    min.cells <- filter1.params$min.cells
    max.cells <- filter1.params$max.cells
    min.reads <- filter1.params$min.reads
    dataset <- gene_filter1(dataset, min.cells, max.cells, min.reads)
    
    dataset <- log2(1 + dataset)

    outl.res <- list()
    
    for(i in clusts) {
        print(i)
        # reduce p dimensions by using
        t <- tryCatch({
            PcaHubert(dataset[ , colnames(dataset) == i])
        }, warning = function(cond) {
            message(cond)
        }, error = function(cond) {
            message(paste0("No outliers detected in cluster ", i, ". Distribution of gene expression in cells is too skewed towards 0."))
            return(NULL)
        })
        if(class(t) != "NULL") {
            # degrees of freedom used in mcd and chisquare distribution
            if(dim(t@loadings)[1] <= 6) {
                message(paste0("No outliers detected in cluster ", i, ". Small number of cells in the cluster."))
                out <- rep(0, dim(dataset[ , colnames(dataset) == i])[2])
                names(out) <- rep(i, dim(dataset[ , colnames(dataset) == i])[2])
                outl.res[[i]] <- out
            } else {
                df <- ifelse(dim(t@loadings)[2] > 3, 3, dim(t@loadings)[2])
                mcd <- tryCatch({
                    covMcd(t@loadings[ , 1:df])
                }, warning = function(cond) {
                    message(cond)
                }, error = function(cond) {
                    message("No outliers detected in the cluster. Error in MCD.")
                    return(NULL)
                })
                if(class(mcd) != "NULL") {
                    # sqrt(mcd$mah) - sqrt of robust distance
                    # sqrt(qchisq(.95, df = length(mcd$best))) - sqrt of 97.5% quantile of a
                    # chi-squared distribution with p degrees of freedom
                    outliers <- sqrt(mcd$mah) - sqrt(qchisq(.9999, df = df))
                    outliers[which(outliers < 0)] <- 0
                    outl.res[[i]] <- outliers
                } else {
                    out <- rep(0, dim(dataset[ , colnames(dataset) == i])[2])
                    names(out) <- rep(i, dim(dataset[ , colnames(dataset) == i])[2])
                    outl.res[[i]] <- out
                }
            }
        } else {
            out <- rep(0, dim(dataset[ , colnames(dataset) == i])[2])
            names(out) <- rep(i, dim(dataset[ , colnames(dataset) == i])[2])
            outl.res[[i]] <- out
        }
    }
    return(outl.res)
}
