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
    
    if(dim(dataset)[2] > 300) {
        original_dataset <- dataset
        dataset <- dataset[, sample(1:dim(dataset)[2], 300)]
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
    labs = foreach(i = 1:dim(hash.table)[1],
                   .packages = "clustools",
                   .combine = rbind) %dopar% {
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

            return(list(dat, labs))
        })
    }
    
    
    
    # stop local cluster
    stopCluster(cl)
    
    show_consensus(filename, distances, dimensionality.reductions,
                   cbind(all.combinations, cons),
                   dataset, original_dataset)
}

show_consensus <- function(filename, distances, dimensionality.reductions, cons.table, dataset, original_dataset) {
    
    dist.opts <- strsplit(unlist(cons.table[,1]), " ")
    dim.red.opts <- strsplit(unlist(cons.table[,2]), " ")
    
    distances <- as.list(distances)
    names(distances) <- distances
    
    dimensionality.reductions <- as.list(dimensionality.reductions)
    names(dimensionality.reductions) <- dimensionality.reductions
    
    colour.pallete <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(median(as.numeric(unlist(cons.table[,3]))))
    plot.width <- 600
    plot.height <- 600
    
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
                div("\n"), div("\n"),
                actionButton("get_de_genes", label = "Get DE genes"),
                div("\n"), div("\n"),
                actionButton("get_mark_genes", label = "Get Marker genes"),
                div("\n"), div("\n"),
                actionButton("svm", label = "Run SVM"),
                
                div("\n"), div("\n"),
                downloadLink('datalink', label = "Download Labels"),
                downloadLink('svm', label = "Download Predicted Labels")
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel("Consensus Matrix", plotOutput('plot')),
                    tabPanel("Expression matrix", plotOutput('matrix')),
                    tabPanel("Cell Labels", div(htmlOutput('labels'), style = "font-size:80%")),
                    tabPanel("DE genes", plotOutput('de_genes')),
                    tabPanel("Marker genes", plotOutput('mark_genes'))
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
            
            sort_cells_by_clusters <- isolate(function(labs.table) {
                labs <- NULL
                gaps_col <- NULL
                gap <- 0
                for(i in 1:input$clusters) {
                    ind <- as.numeric(unlist(strsplit(as.character(labs.table[i, ]), " ")))
                    gap <- gap + length(ind)
                    labs <- c(labs, ind)
                    if(i != input$clusters) {
                        gaps_col <- c(gaps_col, gap)
                    }
                }
                return(list(labs, gaps_col))
            })
            
            output$plot <- renderPlot({
                d <- get_consensus()[[1]]
                show_labs <- TRUE
                if(dim(d)[1] > 80) {
                    show_labs <- FALSE
                }
                pheatmap(d,
                         color = colour.pallete,
                         cutree_rows = input$clusters, cutree_cols = input$clusters,
                         show_rownames = show_labs, show_colnames = show_labs)
            }, height = plot.height, width = plot.width)
            
            output$matrix <- renderPlot({
                d <- get_consensus()[[2]]
                d <- sort_cells_by_clusters(d)
                labs <- d[[1]]
                gaps_col <- d[[2]]
                pheatmap(dataset[ , labs], cluster_cols = F, gaps_col = gaps_col,
                         color = colour.pallete,
                         kmeans_k = 10, show_rownames = F, show_colnames = F)
            }, height = plot.height, width = plot.width)
            
            get_de_genes <- eventReactive(input$get_de_genes, {
                d <- get_consensus()[[2]]
                labs <- rep(1, dim(dataset)[2])
                for(i in 2:input$clusters) {
                    ind <- as.numeric(unlist(strsplit(as.character(d[i, ]), " ")))
                    labs[ind] <- i
                }
                res <- kruskal_statistics(dataset, labs)
                res <- head(res, 70)
                d <- sort_cells_by_clusters(d)
                labs <- d[[1]]
                gaps_col <- d[[2]]
                d <- dataset[rownames(dataset) %in% names(res), labs]
                d <- d[names(res), ]
                
                p.value.ann <- split(res, ceiling(seq_along(res)/17))
                p.value.ranges <- as.vector(unlist(lapply(p.value.ann, function(x){rep(max(x), length(x))})))
                p.value.ranges <- format(p.value.ranges, scientific = T, digits = 2)
                p.value.ranges <- paste("<", p.value.ranges, sep=" ")
                
                p.value.ann <- data.frame(p.value = factor(p.value.ranges, levels = unique(p.value.ranges)))
                rownames(p.value.ann) <- names(res)
                
                pheatmap(d,
                         color = colour.pallete, show_colnames = F,
                         cluster_cols = F, cluster_rows = F,
                         gaps_col = gaps_col, annotation_row = p.value.ann)
            })
            
            get_mark_genes <- eventReactive(input$get_mark_genes, {
                d <- get_consensus()[[2]]
                labs <- rep(1, dim(dataset)[2])
                for(i in 2:input$clusters) {
                    ind <- as.numeric(unlist(strsplit(as.character(d[i, ]), " ")))
                    labs[ind] <- i
                }
                res <- get_marker_genes(dataset, labs)
                res1 <- NULL
                for(i in 1:input$clusters) {
                    tmp <- res[res[,2] == i, ]
                    if(dim(tmp)[1] > 10) {
                        tmp <- tmp[1:10, ]
                    }
                    res1 <- rbind(res1, tmp)
                }
                
                d <- sort_cells_by_clusters(d)
                labs <- d[[1]]
                gaps_col <- d[[2]]
                d <- dataset[rownames(dataset) %in% rownames(res1), labs]
                d <- d[rownames(res1), ]
                
                pheatmap(d,
                         color = colour.pallete, show_colnames = F,
                         cluster_cols = F, cluster_rows = F,
                         gaps_col = gaps_col)
            })
            
            get_svm <- eventReactive(input$svm, {
                prediction <- support_vector_machines()
            })
            
            output$de_genes <- renderPlot({
                get_de_genes()
            }, height = plot.height, width = plot.width)
            
            output$mark_genes <- renderPlot({
               get_mark_genes()
            }, height = plot.height, width = plot.width)
            
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
            output$svmlink <- downloadHandler(
                filename <- function() {
                    paste0("k=", input$clusters, "-svm-labels.csv")
                },
                content <- function(file) {
                    d <- get_svm()
                    write.table(d, file = file)
                }
            )
        }
    )
}
