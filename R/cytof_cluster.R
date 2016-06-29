#' Subset detection by clustering
#' 
#' Apply clustering algorithms to detect cell subsets. \code{DensVM} and \code{ClusterX} 
#' clustering is based on the transformend ydata and use xdata to train the model. 
#' While \code{Rphenograph} directly works on the high dimemnional xdata. \code{FlowSOM} is 
#' integrated from FlowSOM pacakge (https://bioconductor.org/packages/release/bioc/html/FlowSOM.html).
#' 
#' @param ydata A matrix of the dimension reduced data.
#' @param xdata A matrix of the expression data.
#' @param method Cluster method including \code{DensVM}, \code{densityClustX}, \code{Rphenograph} and \code{FlowSOM}.
#' @param FlowSOM_k Number of clusters for meta clustering in FlowSOM.
#' 
#' @return a vector of the clusters assigned for each row of the ydata
#' @export
#' @examples
#' d<-system.file('extdata', package='cytofkit')
#' fcsFile <- list.files(d, pattern='.fcs$', full=TRUE)
#' parameters <- list.files(d, pattern='.txt$', full=TRUE)
#' markers <- as.character(read.table(parameters, sep = "\t", header = TRUE)[, 1])
#' xdata <- cytof_exprsMerge(fcsFile, markers = markers, mergeMethod = 'fixed', fixedNum = 100)
#' ydata <- cytof_dimReduction(xdata, method = "tsne")
#' clusters <- cytof_cluster(ydata, xdata, method = "ClusterX")
cytof_cluster <- function(ydata = NULL, 
                          xdata = NULL, 
                          method = c("Rphenograph", "ClusterX", "DensVM", "FlowSOM", "NULL"),
                          FlowSOM_k = 40){
    
    method = match.arg(method)
    if(method == "NULL"){
        return(NULL)
    }
    switch(method, 
           Rphenograph = {
               cat("  Runing PhenoGraph...")
               clusters <- as.numeric(membership(Rphenograph(xdata, k=30)))
           },
           ClusterX = {
               cat("  Runing ClusterX...")
               clusters <- ClusterX(ydata, gaussian=TRUE, alpha = 0.001, detectHalos = FALSE)$cluster
           },
           DensVM = {
               cat("  Runing DensVM...")
               clusters <- DensVM(ydata, xdata)$cluster$cluster
           },
           FlowSOM = {
               cat("  Runing FlowSOM...")
               clusters <- FlowSOM_integrate2cytofkit(xdata, FlowSOM_k)
           })
    
    if( length(clusters) != ifelse(is.null(ydata), nrow(xdata), nrow(ydata)) ){
        message("Cluster is not complete, cluster failed, try other cluster method!")
        return(NULL)
    }else{
        if(!is.null(xdata) && !is.null(row.names(xdata))){
            names(clusters) <- row.names(xdata)
        }else if(!is.null(ydata) && !is.null(row.names(ydata))){
            names(clusters) <- row.names(ydata)
        }
        cat(" DONE!\n")
        return(clusters)
    }
}


#' FlowSOM algorithm
#' 
#' @param xdata Input data matrix.
#' @param k Number of clusters.
#' @param ... Other parameters passed to SOM.
#' 
#' @noRd
#' @importFrom FlowSOM SOM metaClustering_consensus
FlowSOM_integrate2cytofkit <- function(xdata, k, ...){
    cat("    Building SOM...\n")
    xdata <- as.matrix(xdata)
    
    ord <- tryCatch({
        map <- SOM(xdata, silent = TRUE, ...)
        cat("    Meta clustering to", k, "clusters...\n")
        metaClusters <- suppressMessages(metaClustering_consensus(map$codes, k = k))
        cluster <- metaClusters[map$mapping[,1]]
    }, error=function(cond) {
        message("Run FlowSOM failed")
        message("Here's the error message:")
        message(cond)
        return(NULL)
    }) 
    
    if(is.null(ord)){
        cluster <- NULL
    }else{
        if(length(ord) != nrow(xdata)){
            message("Run FlowSOM failed!")
            return(NULL)
        }
        cluster <- ord
    }
    
    return(cluster)
}
