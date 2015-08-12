all_clusterings <- function(filename, ks) {
    if(filename %in% c("quake", "quake_all_fpkm", "quake_all_read", "sandberg", "sandberg_all_read", "sandberg_all_rpkm",
                       "bernstein", "linnarsson", "zhong", "kirschner")) {
        dataset <- get(filename)
    } else {
        dataset <- read.table(filename, header = T)
        if(filename == "../clustering/tallulah/data/Bergiers_Exp1_FeatureCounts_BaseGenomeAnnotations_No_Multimapping.out" |
           filename == "../clustering/tallulah/data/Bergiers_Exp2_FeatureCounts_BaseGenomeAnnotations_No_Multimapping.out") {
            dataset <- dataset[ , 2:dim(dataset)[2]]
        }
    }
    
    distances <- c("euclidean", "pearson", "spearman")
    dimensionality.reductions <- c("pca", "spectral")
    
    cat("Preliminary gene filtering...\n")
    filter1.params <- filter1_params(dataset)
    min.cells <- filter1.params$min.cells
    max.cells <- filter1.params$max.cells
    min.reads <- filter1.params$min.reads
    dataset <- gene_filter1(dataset, min.cells, max.cells, min.reads)
    
    if(filename == "linnarsson" | filename == "kirschner") {
        dataset <- dataset[, sample(1:dim(dataset)[2], 500)]
    }
    
    n.cells <- dim(dataset)[2]
    n.dim <- floor(0.05 * n.cells) : ceiling(0.08 * n.cells)
    
    if(length(n.dim) > 15) {
        n.dim <- sample(n.dim, 15)
    }
    
    hash.table <- expand.grid(distan = distances,
                              dim.red = dimensionality.reductions,
                              k = c(min(ks) - 1, ks),
                              n.dim = n.dim, stringsAsFactors = F)
    
    cat("Log2-transforming data...\n")
    if(filename != "bernstein") {
        dataset <- log2(1 + dataset)
    }
    
    # register local cluster
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl, cores = detectCores() - 1)
    
    cat("Calculating distance matrices...\n")
    dists = foreach(i = distances, .packages = "clustools") %dopar% {
        try({
            calculate_distance(dataset, i)
        })
    }
    names(dists) <- distances
    
    cat("Performing kmeans clusterings...\n")
    labs = foreach(i = 1:dim(hash.table)[1], .packages = "clustools", .combine = rbind) %dopar% {
        try({
            t <- transformation(get(hash.table[i, 1], dists), hash.table[i, 2])[[1]]
            paste(kmeans(t[, 1:hash.table[i, 4]],
                         hash.table[i, 3],
                         iter.max = 1e+09,
                         nstart = 1000)$cluster,
                  collapse = " ")
        })
    }
    
    res <- cbind(hash.table, labs)
    res$labs <- as.character(res$labs)
    rownames(res) <- NULL
    
    cat("Computing consensus clusterings and labels...\n")
    all.combinations <- NULL
    for(k in c(min(ks) - 1, ks)) {
        for(i in 1:length(distances)) {
            for(j in 1:length(dimensionality.reductions)) {
                dist.combs <- combn(distances, i)
                dim.red.combs <- combn(dimensionality.reductions, j)
                for(m in 1:dim(dist.combs)[2]) {
                    for(n in 1:dim(dim.red.combs)[2]) {
                        all.combinations <- rbind(
                            all.combinations,
                            cbind(paste(dist.combs[, m], collapse = " "),
                                  paste(dim.red.combs[, n], collapse = " "),
                                  as.numeric(k)))
                    }
                }
            }
        }
    }
    
    cons = foreach(i = 1:dim(all.combinations)[1], .packages = "clustools") %dopar% {
        try({
            d <- res[res$distan %in% strsplit(all.combinations[i, 1], " ")[[1]] &
                         res$dim.red %in% strsplit(all.combinations[i, 2], " ")[[1]] &
                         res$k == as.numeric(all.combinations[i, 3]), ]
            
            dat <- consensus_clustering(d$labs)
            
            diss <- dist(dat)
            clust <- hclust(diss)
            clusts <- cutree(clust, k = as.numeric(all.combinations[i, 3]))
            
            labs <- NULL
            for(j in 1:as.numeric(all.combinations[i, 3])) {
                labs <- rbind(labs, paste(names(clusts[clusts == j]), collapse = " "))
            }
            labs <- as.data.frame(labs)
            colnames(labs) <- "Labels"
            # rownames(labs) <- paste("Cluster", 1:as.numeric(all.combinations[i, 3]))
            
            return(list(dat, labs))
        })
    }
    
    
    
    # stop local cluster
    stopCluster(cl)
    
    show_consensus(filename, distances, dimensionality.reductions, cbind(all.combinations, cons), dataset)
}

show_consensus <- function(filename, distances, dimensionality.reductions, cons.table, dataset) {
    
    dist.opts <- strsplit(unlist(cons.table[,1]), " ")
    dim.red.opts <- strsplit(unlist(cons.table[,2]), " ")
    
    distances <- as.list(distances)
    names(distances) <- distances
    
    dimensionality.reductions <- as.list(dimensionality.reductions)
    names(dimensionality.reductions) <- dimensionality.reductions
    
    shinyApp(
        ui = pageWithSidebar(
            headerPanel(
                HTML("SC<sup>3</sup> - Single-Cell Consensus Clustering")
            ),
            sidebarPanel(
                
                sliderInput("clusters", label = "k",
                            min = min(as.numeric(unlist(cons.table[,3]))) + 1,
                            max = max(as.numeric(unlist(cons.table[,3]))),
                            value = min(as.numeric(unlist(cons.table[,3])))),
                
                checkboxGroupInput("distance", label = "Distance metrics",
                                   choices = distances,
                                   selected = distances[1]),
                
                checkboxGroupInput("dimRed", label = "Dimensionality reduction",
                                   choices = dimensionality.reductions,
                                   selected = dimensionality.reductions[1]),
                
                div("\n"),
                downloadLink('datalink', label = "Download Labels")
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Consensus Matrix", plotOutput('plot')),
                    tabPanel("Expression matrix", htmlOutput('matrix')),
                    tabPanel("Cell Labels", div(htmlOutput('labels'), style = "font-size:80%"))
                )
            )
        ),
        server = function(input, output) {
            get_consensus <- reactive({
                res <- cons.table[unlist(lapply(dist.opts, function(x){setequal(x, input$distance)})) &
                                      unlist(lapply(dim.red.opts, function(x){setequal(x, input$dimRed)})) &
                                      as.numeric(cons.table[ , 3]) == input$clusters, 4]
                return(res[[1]])
            })
            get_consensus_1 <- reactive({
                res <- cons.table[unlist(lapply(dist.opts, function(x){setequal(x, input$distance)})) &
                                      unlist(lapply(dim.red.opts, function(x){setequal(x, input$dimRed)})) &
                                      as.numeric(cons.table[ , 3]) == (input$clusters - 1), 4]
                return(res[[1]])
            })
            output$plot <- renderPlot({
                d <- get_consensus()[[1]]
                show_labs <- TRUE
                if(dim(d)[1] > 80) {
                    show_labs <- FALSE
                }
                pheatmap(d,
                         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(median(as.numeric(unlist(cons.table[,3])))),
                         cutree_rows = input$clusters, cutree_cols = input$clusters,
                         show_rownames = show_labs, show_colnames = show_labs)
            }, height = 600, width = 600)
            output$matrix <- renderPlot({
                d <- get_consensus()[[2]]
                labs <- NULL
                for(i in 1:input$clusters) {
                    ind <- as.numeric(unlist(strsplit(as.character(d[i, ]), " ")))
                    labs <- c(labs, ind)
                }
                pheatmap(dataset[ , labs],
                         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(median(as.numeric(unlist(cons.table[,3])))),
                         cutree_rows = input$clusters, cutree_cols = input$clusters,
                         show_rownames = F, show_colnames = F)
            }, height = 600, width = 600)
            output$labels <- renderUI({
                d <- get_consensus()[[2]]
                d1 <- get_consensus_1()[[2]]
                
                labs1 <- list()
                cols <- brewer.pal(input$clusters - 1, "Paired")
                for(i in 1:(input$clusters - 1)) {
                    col <- cols[i]
                    ind <- unlist(strsplit(as.character(d1[i, ]), " "))
                    for(j in ind) {
                        labs1[[j]] <- paste0("<font color=\"", col, "\">", j, "</font>")
                    }
                }
                labs <- "<br/>"
                for(i in 1:input$clusters) {
                    ind <- unlist(strsplit(as.character(d[i, ]), " "))
                    for(j in ind) {
                        labs <- c(labs, labs1[[j]])
                    }
                    labs <- c(labs, c("<br/>", "<hr>"))
                }
                
                HTML(paste0(labs))
            })
            output$datalink <- downloadHandler(
                filename <- function() {
                    paste0("Exp", input$replicate, "-", input$norm, "-k=", input$clusters, "-labels.csv")
                },
                content <- function(file) {
                    diss <- dist(data())
                    clust <- hclust(diss)
                    clusts <- cutree(clust, k = input$clusters)
                    d <- rep(0, length(clusts))
                    for(i in 2:input$clusters) {
                        inds <- as.numeric(names(clusts[clusts == i]))
                        d[inds] <- i
                    }
                    
                    write.table(d, file = file)
                }
            )
        }
    )
}
