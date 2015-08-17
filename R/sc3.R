all_clusterings <- function(filename, ks) {
    if(filename %in% c("quake", "quake_all_fpkm", "quake_all_read", "sandberg", "sandberg_all_read", "sandberg_all_rpkm",
                       "bernstein", "linnarsson", "zhong", "kirschner")) {
        dataset <- get(filename)
    } else {
        dataset <- read.table(filename, header = T)
        if(filename == "../clustering/tallulah/data/Bergiers_Exp1_FeatureCounts_BaseGenomeAnnotations_No_Multimapping.out" |
           filename == "../clustering/tallulah/data/Bergiers_Exp2_FeatureCounts_BaseGenomeAnnotations_No_Multimapping.out" |
           filename == "../clustering/tallulah/data/Bergiers_Exp1_FeatureCounts_BaseGenomeAnnotations_No_Multimapping_RUVnorm.txt" |
           filename == "../clustering/tallulah/data/Bergiers_Exp2_FeatureCounts_BaseGenomeAnnotations_No_Multimapping_RUVnorm.txt" |
           filename == "../clustering/tallulah/data/Bergiers_Exp1_FeatureCounts_BaseGenomeAnnotations_No_Multimapping_SFnorm.txt" |
           filename == "../clustering/tallulah/data/Bergiers_Exp2_FeatureCounts_BaseGenomeAnnotations_No_Multimapping_SFnorm.txt") {
            rownames(dataset) <- dataset[, 1]
            dataset <- dataset[ , 2:dim(dataset)[2]]
        }
    }
    
    svm.num.cells <- 500
    distances <- c("euclidean", "pearson", "spearman")
    dimensionality.reductions <- c("pca", "spectral")
    
    cat("1. Preliminary gene filtering...\n")
    filter1.params <- filter1_params(dataset)
    min.cells <- filter1.params$min.cells
    max.cells <- filter1.params$max.cells
    min.reads <- filter1.params$min.reads
    dataset <- gene_filter1(dataset, min.cells, max.cells, min.reads)
    
    if(dim(dataset)[2] > svm.num.cells) {
        cat("\n")
        cat("Your dataset contains more than 1000 cells, therefore clustering
            wil be performed on a random sample of 1000 cells, the rest of the
            cells will be predicted using SVM.")
        cat("\n")
        cat("\n")
        working.sample <- sample(1:dim(dataset)[2], svm.num.cells)
        study.dataset <- dataset[ , setdiff(1:dim(dataset)[2], working.sample)]
        dataset <- dataset[, working.sample]
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
    
    cat("2. Log2-transforming data...\n")
    if(filename != "bernstein") {
        dataset <- log2(1 + dataset)
    }
    
    # register local cluster
    cl <- makeCluster(detectCores() - 1)
    registerDoParallel(cl, cores = detectCores() - 1)
    
    cat("3. Calculating distance matrices...\n")
    dists = foreach(i = distances, .packages = "clustools") %dopar% {
        try({
            calculate_distance(dataset, i)
        })
    }
    names(dists) <- distances
    
    cat("4. Performing dimensionality reduction and kmeans clusterings...\n")
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
    
    cat("5. Computing consensus matrix and labels...\n")
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
            for(j in unique(clusts[clust$order])) {
                labs <- rbind(labs, paste(names(clusts[clusts == j]), collapse = " "))
            }

            labs <- as.data.frame(labs)
            colnames(labs) <- "Labels"

            return(list(dat, labs, clust))
        })
    }
    
    
    
    # stop local cluster
    stopCluster(cl)
    
    show_consensus(filename, distances, dimensionality.reductions,
                   cbind(all.combinations, cons),
                   dataset, study.dataset)
}

show_consensus <- function(filename, distances, dimensionality.reductions, cons.table, dataset, study.dataset) {
    
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
        ui = fluidPage(
            headerPanel(
                HTML("SC<sup>3</sup> - Single-Cell Consensus Clustering")
            ),
            sidebarPanel(
                h4("1. Clustering"),
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
                h4("2. Gene identificatiion"),
                p("\n\n"),
                actionButton("get_de_genes", label = "Get DE genes"),
                p("\n\n"),
                actionButton("get_mark_genes", label = "Get Marker genes"),

                h4("3. Save results"),
                p("\n\n"),
                downloadLink('datalink', label = "Save cell labels")
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
            
            output$plot <- renderPlot({
                d <- get_consensus()
                hc <- d[[3]]
                d <- d[[1]]
                show_labs <- TRUE
                if(dim(d)[1] > 80) {
                    show_labs <- FALSE
                }
                withProgress(message = 'Plotting...', value = 0, {
                    pheatmap(d,
                             color = colour.pallete,
                             cluster_rows = hc,
                             cutree_rows = input$clusters, cutree_cols = input$clusters,
                             show_rownames = show_labs, show_colnames = show_labs)
                })
            }, height = plot.height, width = plot.width)
            
            output$matrix <- renderPlot({
                hc <- get_consensus()[[3]]
                withProgress(message = 'Plotting...', value = 0, {
                    pheatmap(dataset, cluster_cols = hc,
                             cutree_cols = input$clusters,
                             color = colour.pallete,
                             kmeans_k = 10, show_rownames = F, show_colnames = F)
                })
            }, height = plot.height, width = plot.width)
            
            get_de_genes <- eventReactive(input$get_de_genes, {
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                hc <- get_consensus()[[3]]
                clusts <- cutree(hc, input$clusters)
                res <- kruskal_statistics(dataset, clusts)
                res <- head(res, 70)
                d <- dataset[rownames(dataset) %in% names(res), ]
                d <- d[names(res), ]
                
                p.value.ann <- split(res, ceiling(seq_along(res)/17))
                p.value.ranges <- as.vector(unlist(lapply(p.value.ann, function(x){rep(max(x), length(x))})))
                p.value.ranges <- format(p.value.ranges, scientific = T, digits = 2)
                p.value.ranges <- paste("<", p.value.ranges, sep=" ")
                
                p.value.ann <- data.frame(p.value = factor(p.value.ranges, levels = unique(p.value.ranges)))
                rownames(p.value.ann) <- names(res)
                
                pheatmap(d,
                         color = colour.pallete, show_colnames = F,
                         cluster_cols = hc,
                         cutree_cols = input$clusters, cluster_rows = F,
                         annotation_row = p.value.ann)
            })
            
            get_mark_genes <- eventReactive(input$get_mark_genes, {
                validate(
                    need(try(!is.null(rownames(dataset))), "\nNo gene names provided in the input expression matrix!")
                )
                hc <- get_consensus()[[3]]
                clusts <- cutree(hc, input$clusters)
                res <- get_marker_genes(dataset, clusts)
                res1 <- NULL
                for(i in 1:input$clusters) {
                    tmp <- res[res[,2] == i, ]
                    if(dim(tmp)[1] > 10) {
                        tmp <- tmp[1:10, ]
                    }
                    res1 <- rbind(res1, tmp)
                }
                
                d <- dataset[rownames(dataset) %in% rownames(res1), ]
                d <- d[rownames(res1), ]
                
                pheatmap(d,
                         color = colour.pallete, show_colnames = F,
                         cluster_cols = hc,
                         cutree_cols = input$clusters, cluster_rows = F)
            })
            
            get_svm <- reactive({
                withProgress(message = 'Running SVM...', value = 0, {
                    d <- get_consensus()[[2]]
                    labs <- rep(1, dim(dataset)[2])
                    for(i in 2:input$clusters) {
                        ind <- as.numeric(unlist(strsplit(as.character(d[i, ]), " ")))
                        labs[ind] <- i
                    }
                    colnames(dataset) <- labs
                    cat("svm prediction started")
                    prediction <- support_vector_machines1(dataset, study.dataset, "linear")
                    cat("svm prediction finished")
                })
            })
            
            output$de_genes <- renderPlot({
                withProgress(message = 'Calculating...', value = 0, {
                    get_de_genes()
                })
            }, height = plot.height, width = plot.width)
            
            output$mark_genes <- renderPlot({
                withProgress(message = 'Calculating...', value = 0, {
                    get_mark_genes()
                })
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
                filename = function() {
                    paste0("k=", input$clusters, "-labels.csv")
                },
                content = function(file) {
                    if(dim(dataset)[2] > 20) {
                        write.table(get_svm(), file = file)
                    } else {
                        write.table(get_consensus()[[2]][[1]], file = file)
                    }
                }
            )
        }
    )
}
