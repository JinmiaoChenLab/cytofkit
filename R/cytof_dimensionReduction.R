#' Dimension reduction of CyTof expression data 
#' 
#' Apply dimension reduction on the CyTof expression data, 
#' with method \code{isomap}, \code{pca}, \code{lle}, or \code{tsne}. 
#' 
#' @param data An expression data matrix
#' @param method Method chosed for dimensition reduction, must be one of \code{isomap}, \code{pca}, \code{lle} or \code{tsne} 
#' @param lle_k Number of neighbours, parameter for \code{lle} method
#' @param lle_m Intrinsic dimension of the data, parameter for \code{lle} method
#' @param distMethod Method for distance calcualtion
#' @param isomap_k Number of shortest dissimilarities retained for a point, parameter for \code{isomap} method
#' @param isomap_ndim Number of axes in metric scaling, parameter for \code{isomap} method
#' @param isomapFragmentOK What to do if dissimilarity matrix is fragmented, parameter for \code{isomap} method
#' @return a matrix of the dimension reducted data, with colnames and rownames(if have, same as the input)
#' @author Chen Jinmiao
#' @importFrom lle lle
#' @importFrom vegan vegdist spantree isomap
#' @importFrom Rtsne Rtsne
#' @import stats
#' @export
#' @examples
#' data(iris)
#' in_data <- iris[, 1:4]
#' out_data <- cytof_dimReduction(in_data)
#' 
cytof_dimReduction <- function(data, method = "tsne", lle_k = 12, 
    lle_m = 2, distMethod = "euclidean", isomap_k = 5, isomap_ndim = NULL, 
    isomapFragmentOK = TRUE) {
    
    if (method == "pca") {
        res <- prcomp(as.matrix(data), scale = TRUE)
        colnames(res$x) <- paste("PCA_dim", c(1:ncol(res$x)), 
            sep = "")
        rownames(res$x) <- row.names(data)
        return(res$x)
    }
    if (method == "lle") {
        res <- lle(X = as.matrix(data), m = lle_m, k = lle_k, 
            reg = 2, ss = FALSE, id = TRUE, v = 0.9)
        colnames(res$Y) <- paste("LLE_dim", c(1:ncol(res$Y)), 
            sep = "")
        rownames(res$Y) <- row.names(data)
        return(res$Y)
    }
    if (method == "isomap") {
        if (is.null(isomap_ndim)) {
            isomap_ndim <- dim(data)[2]
        }
        dis <- vegdist(as.matrix(data), method = distMethod)
        ord <- isomap(dis, ndim = isomap_ndim, k = isomap_k, 
            fragmentedOK = isomapFragmentOK)
        colnames(ord$points) <- paste("ISOMAP_dim", c(1:ncol(ord$points)), 
            sep = "")
        rownames(ord$points) <- row.names(data)
        return(ord$points)
    }
    if (method == "tsne") {
        tsne_out <- Rtsne(as.matrix(data), initial_dims = dim(as.matrix(data))[2], 
            dims = 2, perplexity = 30, theta = 0.5, check_duplicates = FALSE, 
            pca = TRUE)
        
        mapped <- tsne_out$Y
        colnames(mapped) <- paste("tsne", c(1:ncol(mapped)), 
            sep = "_")
        rownames(mapped) <- row.names(data)
        return(mapped)
    }
} 
