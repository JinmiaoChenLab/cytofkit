#' Progression estimation of cytof expression data 
#' 
#' Infer the progression based on the relationship of cell subsets estimated 
#' using ISOMAP or Diffusion map.
#' 
#' @param data Expression data matrix.
#' @param cluster A vector of cluster results for the data.
#' @param method Method for estimation of cell progression, isomap or diffusionmap.
#' @param distMethod Method for distance calculation, default is "euclidean", other choices like "manhattan", "cosine", "rankcor".
#' @param out_dim Number of transformed dimensions choosen for output.
#' @param clusterSampleMethod Cluster sampling method including \code{ceil}, \code{all}, \code{min}, \code{fixed}. The default option is 
#' \code{ceil}, up to a fixed number (specified by \code{fixedNum}) of cells are sampled without replacement from each cluster and combined for analysis.
#' \code{all}: all cells from each cluster are combined for analysis. 
#' \code{min}: The minimum number of cells among all clusters are sampled from cluster and combined for analysis. 
#' \code{fixed}: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each cluster and combined for analysis.
#' @param clusterSampleSize The number of cells to be sampled from each cluster.
#' @param sampleSeed The seed for random down sample of the clusters.
#' 
#' @return a list. Includes: sampleData, sampleCluster and progressionData.
#' 
#' @export
#' @examples
#' d<-system.file('extdata', package='cytofkit')
#' fcsFile <- list.files(d, pattern='.fcs$', full=TRUE)
#' parameters <- list.files(d, pattern='.txt$', full=TRUE)
#' markers <- as.character(read.table(parameters, header = TRUE)[, 1])
#' xdata <- cytof_exprsMerge(fcsFile, mergeMethod = 'fixed', fixedNum = 2000)
#' clusters <- cytof_cluster(xdata = xdata, method = "Rphenograph")
#' prog <- cytof_progression(data = xdata, cluster = clusters, clusterSampleSize = 100)
#' d <- as.data.frame(cbind(prog$progressionData, cluster = factor(prog$sampleCluster)))
#' cytof_clusterPlot(data =d, xlab = "diffusionmap_1", ylab="diffusionmap_2",
#'                   cluster = "cluster", sampleLabel = FALSE)
cytof_progression <- function(data, cluster, 
                              method=c("diffusionmap", "isomap", "NULL"), 
                              distMethod = "euclidean",
                              out_dim = 2,
                              clusterSampleMethod = c("ceil", "all", "fixed", "min"),
                              clusterSampleSize = 500, 
                              sampleSeed = 123){
    
    data <- as.matrix(data)
    method <- match.arg(method)
    if(method == "NULL"){
        return(NULL)
    }
    clusterSampleMethod <- match.arg(clusterSampleMethod)
    if(missing(cluster)){
        message("No cluster vector provided, take all data for estimation.")
        clusterSampleMethod <- "all"
        cluster <- rep(1, length.out = nrow(data))
    }
    if(is.numeric(sampleSeed))
        set.seed(sampleSeed) # Set a seed if you want reproducible results
    
    cellClusterList <- split(1:nrow(data), cluster)
    switch(clusterSampleMethod,
           ceil = {
               mergeFunc <- function(x) {
                   if (length(x) < clusterSampleSize) {
                       x
                   } else {
                       sample(x, size = clusterSampleSize, replace = FALSE)
                   }
               }
               sampleCellID <- do.call(base::c, lapply(cellClusterList, mergeFunc))
           },
           all = {
               sampleCellID <- do.call(base::c, cellClusterList)
           },
           fixed = {
               mergeFunc <- function(x) {
                   sample(x, size = clusterSampleSize, replace = ifelse(length(x) < clusterSampleSize, TRUE, FALSE))
               }
               sampleCellID <- do.call(base::c, lapply(cellClusterList, mergeFunc))
           },
           min = {
               minSize <- min(sapply(cellClusterList, length))
               mergeFunc <- function(x) {
                   sample(x, size = minSize, replace = FALSE)
               }
               sampleCellID <- do.call(base::c, lapply(cellClusterList, mergeFunc))
           })
    
    sampleData <- data[sampleCellID, ,drop=FALSE]
    nCluster <- cluster[sampleCellID]
    progressionData <- cytof_dimReduction(sampleData, method = method, 
                                          distMethod = distMethod, out_dim = out_dim)
    
    progressRes <- list(sampleData = sampleData, 
                        sampleCluster = nCluster, 
                        progressionData = progressionData)
    return(progressRes)
}