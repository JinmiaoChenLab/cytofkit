#' Dimension reduction for high dimensional data 
#' 
#' Apply dimension reduction on the cytof expression data, 
#' with method \code{pca}, \code{tsne}, \code{diffusionmap} or \code{isomap}. 
#' 
#' @param data Input expression data matrix.
#' @param method Method chosed for dimensition reduction, must be one of \code{isomap}, \code{pca} , \code{diffusionmap} or \code{tsne}. 
#' @param out_dim The dimensionality of the output.
#' @param tsneSeed Set a seed if you want reproducible t-SNE results.
#' @param distMethod Method for distance calcualtion, default is "euclidean", other choices like "manhattan", "cosine", "rankcor"....
#' @param isomap_k Number of shortest dissimilarities retained for a point, parameter for \code{isomap} method.
#' @param isomap_ndim Number of axes in metric scaling, parameter for \code{isomap} method.
#' @param isomapFragmentOK What to do if dissimilarity matrix is fragmented, parameter for \code{isomap} method.
#' @param ... Other parameters passed to the method, check \code{\link{Rtsne}}, \code{\link{DiffusionMap}}, \code{\link{isomap}}.
#' @return a matrix of the dimension reducted data, with colnames method_ID, and rownames same as the input data.
#' 
#' @importFrom vegan vegdist spantree isomap
#' @importFrom Rtsne Rtsne
#' @importFrom destiny DiffusionMap
#' @import stats
#' @export
#' @examples
#' data(iris)
#' in_data <- iris[, 1:4]
#' out_data <- cytof_dimReduction(in_data, method = "tsne")
cytof_dimReduction <- function(data, 
                               method = c("tsne", "pca", "isomap", "diffusionmap", "NULL"), 
                               distMethod = "euclidean", 
                               out_dim = 2,
                               tsneSeed = 42,
                               isomap_k = 5, 
                               isomap_ndim = NULL, 
                               isomapFragmentOK = TRUE,
                               ...) {
    
    data <- as.matrix(data)
    rnames <- row.names(data)
    
    method <- match.arg(method)
    if(method == "NULL"){
        return(NULL)
    }
    
    switch(method,
           tsne={
               cat("  Runing t-SNE...with seed", tsneSeed)
               if(is.numeric(tsneSeed))
                   set.seed(tsneSeed) # Set a seed if you want reproducible results
               tsne_out <- Rtsne(data, initial_dims = ncol(data), 
                                 dims = 2, 
                                 check_duplicates = FALSE, 
                                 pca = TRUE, ...)
               mapped <- tsne_out$Y
           },
           pca={
               cat("  Runing PCA...")
               mapped <- prcomp(data, scale = TRUE)$x
           },
           diffusionmap={
               cat("  Runing Diffusion Map...")
               ord <- tryCatch({
                   DiffusionMap(data, distance = distMethod, ...)
                   }, error=function(cond) {
                   message("Run Diffusion Map failed")
                   message("Here's the error message:")
                   message(cond)
                   return(NULL)
                   }) 
               
               if(is.null(ord)){
                   mapped <- NULL
               }else{
                   if(nrow(ord@eigenvectors) != nrow(data) || any(!complete.cases(ord@eigenvectors))){
                       message("Run Diffusion Map failed!")
                       return(NULL)
                   }
                   mapped <- ord@eigenvectors
                   mapped <- apply(mapped, 2, function(x) {
                       ## replace inf value to max finite value
                       x[is.infinite(x)] <- max(x[is.finite(x)])
                       x
                   })
               }
           },
           isomap={
               cat("  Runing ISOMAP...")
               if (is.null(isomap_ndim))
                   isomap_ndim <- ncol(data)
               
               ord <- tryCatch({
                       dis <- vegdist(data, method = distMethod)
                       isomap(dis, ndim = isomap_ndim, k = isomap_k, fragmentedOK = isomapFragmentOK, ...)
                       }, error=function(cond) {
                       message("Run isomap failed")
                       message("Here's the error message:")
                       message(cond)
                       return(NULL)
                       })    
               
               if(is.null(ord)){
                   mapped <- NULL
               }else{
                   if(nrow(ord$points) != nrow(data) || any(!complete.cases(ord$points))){
                       message("Run ISOMAP failed!")
                       return(NULL)
                   }
                   mapped <- ord$points
               }
           })
    
    ## extract out_dim dimensions, organize output
    if(!is.null(mapped)){
        if(ncol(mapped) < out_dim){
            out_dim <- ncol(mapped)
            message("Run ",method," for dimensional reduction, out dimension coerced to ",out_dim)
        }
        mapped <- mapped[ ,c(1:out_dim)]
        colnames(mapped) <- paste(method, c(1:out_dim), sep = "_")
        rownames(mapped) <- rnames
    }
    cat("  DONE\n")
    return(mapped)
} 
