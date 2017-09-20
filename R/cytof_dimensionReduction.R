#' Dimension reduction for high dimension data 
#' 
#' Apply dimension reduction on the cytof expression data, 
#' with method \code{pca}, \code{tsne}, \code{diffusionmap} or \code{isomap}. 
#' 
#' @param data Input expression data matrix.
#' @param markers Selected markers for dimension reduction, either marker names/descriptions or marker IDs.
#' @param method Method chosen for dimensition reduction, must be one of \code{isomap}, \code{pca} , \code{diffusionmap} or \code{tsne}. 
#' @param out_dim The dimensionality of the output.
#' @param tsneSeed Set a seed if you want reproducible t-SNE results.
#' @param distMethod Method for distance calcualtion, default is "euclidean", other choices like "manhattan", "cosine", "rankcor"....
#' @param isomap_k Number of shortest dissimilarities retained for a point, parameter for \code{isomap} method.
#' @param isomap_ndim Number of axes in metric scaling, parameter for \code{isomap} method.
#' @param isomapFragmentOK What to do if dissimilarity matrix is fragmented, parameter for \code{isomap} method.
#' @param ... Other parameters passed to the method, check \code{\link{Rtsne}}, \code{\link{DiffusionMap}}, \code{\link{isomap}}.
#' @return A matrix of the dimension reduced data, with colnames method_ID, and rownames same as the input data.
#' 
#' @importFrom vegan vegdist spantree isomap
#' @importFrom Rtsne Rtsne
#' @importFrom destiny DiffusionMap
#' @importFrom utils compareVersion packageVersion
#' @import stats
#' @export
#' @examples
#' data(iris)
#' in_data <- iris[, 1:4]
#' markers <- colnames(in_data[, 1:4])
#' out_data <- cytof_dimReduction(in_data, markers = markers, method = "tsne")
#' @note Currently, \code{diffusionmap} will not work with R 3.4.0, due to an issue with the latest CRAN release of its dependency \code{\link{igraph}} 
#' If this is the case, consider manually updating \code{\link{igraph}} using;
#' \code{install.packages("https://github.com/igraph/rigraph/releases/download/v1.1.0/igraph_1.1.0.zip", repos=NULL, method="libcurl")}
cytof_dimReduction <- function(data,
                               markers = NULL,
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
    
    ##markers
    if (!(is.null(markers))) {
      if(is.character(markers)){
        right_marker <- markers %in% colnames(data)
        if(!all(right_marker)){
          stop("\n Selected marker(s) is/are not in the input fcs files \n please check your selected marker(s)! \n")
        }else{
          marker_id <- markers
        }
      }else{
        stop("Sorry, input markers cannot be recognized!")
      }
    }else{
      ## NULL default to all
      marker_id <- colnames(data)
    }
    
    marker_filtered_data <- data[, marker_id]
    
    method <- match.arg(method)
    if(method == "NULL"){
        return(NULL)
    }
    
    switch(method,
           tsne={
               cat("  Running t-SNE...with seed", tsneSeed)
               if(is.numeric(tsneSeed))
                   set.seed(tsneSeed) # Set a seed if you want reproducible results
               tsne_out <- Rtsne(marker_filtered_data, initial_dims = ncol(marker_filtered_data), 
                                 dims = 2, 
                                 check_duplicates = FALSE, 
                                 pca = TRUE, ...)
               mapped <- tsne_out$Y
           },
           pca={
               cat("  Running PCA...")
               mapped <- prcomp(marker_filtered_data, scale = TRUE)$x
           },
           diffusionmap={
               cat("  Running Diffusion Map...\n")
               versiontest <- compareVersion(as.character(packageVersion("igraph")), "1.1.0")
               if(versiontest == 0 || versiontest == 1){
                 message("igraph up to date!")
               }else{
                 stop("igraph not at least version 1.1.0! Stopping...")
               }
               ord <- tryCatch({
                   DiffusionMap(marker_filtered_data, distance = distMethod, ...)
                   }, error=function(cond) {
                   message("Run Diffusion Map failed")
                   message("Here's the error message:")
                   message(cond)
                   return(NULL)
                   }) 
               
               if(is.null(ord)){
                   mapped <- NULL
               }else{
                   if(nrow(ord@eigenvectors) != nrow(marker_filtered_data) || any(!complete.cases(ord@eigenvectors))){
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
               cat("  Running ISOMAP...")
               if (is.null(isomap_ndim))
                   isomap_ndim <- ncol(marker_filtered_data)
               
               ord <- tryCatch({
                       dis <- vegdist(marker_filtered_data, method = distMethod)
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
                   if(nrow(ord$points) != nrow(marker_filtered_data) || any(!complete.cases(ord$points))){
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
