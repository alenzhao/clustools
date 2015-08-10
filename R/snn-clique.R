
#########################################################
# This program is part of the SNN-Cliq method           #
# Contact Chen Xu at UNC-Charlotte for more information.# 
#########################################################
#----- example of use------#
#data<-read.table(infile, header=TRUE, sep="\t", row.names=1);
#data<-log2(data+1)
#source('SNN.R')
#SNN(data, edge_file, k=3, distance='euclidean')
#--------------------------#

SNN<-function(data, outfile, k, distance){
    
    if(missing(data)){
        stop(paste("Input data missing.",help,sep="\n"))
    }
    if(missing(outfile)){
        stop(paste("Output file name missing.",help,sep="\n"))
    }
    if(missing(k)){
        k=3
    }
    if(missing(distance)){
        distance<-"euclidean"  # other distance options refer to dist() in R
    }
    m<-as.data.frame(data)
    numSpl<-dim(data)[1]
    m<-dist(data, distance, diag=TRUE, upper=TRUE)
    x<-as.matrix(m)
    IDX<-t(apply(x,1,order)[1:k,]) # knn list
    
    edges<-list()              # SNN graph
    for (i in 1:numSpl){
        j<-i
        while (j<numSpl){
            j<-j+1
            shared<-intersect(IDX[i,], IDX[j,])
            if(length(shared)>0){			
                s<-k-0.5*(match(shared, IDX[i,])+match(shared, IDX[j,]))
                strength<-max(s)
                if (strength>0)
                    edges<-rbind(edges, c(i,j,strength))
            }				
        }
    }
    write.table(edges, outfile, quote=FALSE, sep='\t',col.names=FALSE,row.names=FALSE)
}

run_snn_cliq <- function(path.to.Cliq.py, dataset, name, log.tr, par.k, par.r, par.m) {
#     filter1.params <- filter1_params(dataset)
#     min.cells <- filter1.params$min.cells
#     max.cells <- filter1.params$max.cells
#     min.reads <- filter1.params$min.reads
#     d <- gene_filter1(dataset, min.cells, max.cells, min.reads)
#     
#     if (log.tr) {
#         d <- log2(1 + d)
#     }
    
    d <- t(d)
    
    SNN(d, paste0("snn-clique-", name, ".txt"), k=par.k, distance='euclidean')
    system(paste0("python ", path.to.Cliq.py,"Cliq.py -i ", "snn-clique-",
                  name, ".txt", "  -o ", "res-snn-clique-", name, ".txt -r ",
                  par.r, " -m ", par.m))
    clusts <- read.table(paste0("res-snn-clique-", name, ".txt"))
    
    system(paste0("rm snn-clique-", name, ".txt"))
    system(paste0("rm res-snn-clique-", name, ".txt"))
    
    return(list(c(t(clusts)), adjustedRandIndex(rownames(d), c(t(clusts)))))
}
