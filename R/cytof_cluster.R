#' Density-based local maxima cluster with SVM 
#' 
#' Density-based local maxima peak finding, subpopulation assigning with the power of SVM
#' 
#' @param ydata a matrix of the dimension reduced(transformed) data
#' @param xdata a matrix of the expression data
#' @return a list contains a matrix \code{peakdata} of the peak numbers with different kernel bandwidth, and a matrix \code{clusters} of the cluster results 
#' @export
#' @examples
#' d<-system.file('extdata',package='cytofkit')
#' fcsFile <- list.files(d,pattern='.fcs$',full=TRUE)
#' xdata <- fcs_lgcl_merge(fcsFile, mergeMethod = 'fixed', fixedNum = 100)
#' ydata <- cytof_dimReduction(xdata)
#' clusters <- densVM_cluster(ydata, xdata)
densVM_cluster <- function(ydata, xdata) {
    
    y_range_x <- max(ydata[, 1]) - min(ydata[, 1])
    y_range_y <- max(ydata[, 2]) - min(ydata[, 2])
    y_range <- min(c(y_range_x, y_range_y))
    
    ## Search density peaks of different kernel bandwidth(gamma or
    ## sig_opt)
    num_band <- 20
    kernel_min <- y_range/100
    kernel_max <- y_range/10
    sig_tol_range <- seq(kernel_min, kernel_max, length.out = num_band)
    
    cat("Testing kernel bandwidth for ", num_band, "points in the range min=", 
        kernel_min, " to max=", kernel_max, "\n")
    
    cat("This will take a while...\n")
    peak_list <- mapply(peakFind, sig_tol = sig_tol_range, MoreArgs = list(ydata = ydata), 
        SIMPLIFY = FALSE)
    
    Npeaks <- sapply(peak_list, function(x) {
        x$Npeaks
    })
    peakdata <- data.frame(numpeaks = Npeaks, sig_range = sig_tol_range)
    
    ## Locate plateau(knee)
    flag <- 0
    sig_id <- 0
    for (i in 1:length(sig_tol_range)) {
        if (i > 1) {
            if ((Npeaks[i - 1] - Npeaks[i]) <= 1) {
                sig_opt <- sig_tol_range[i - 1]
                Nsub <- Npeaks[i - 1]
                flag <- 1
                sig_id <- i - 1
                break
            }
        }
    }
    
    ## determine the optimal gamma
    if (flag == 0) {
        cat("Could not locate plateau in the Npeaks vs. sigma graph!\n")
        cat("Consider changing the search space. Increase 'num_band' or change 'kernel_min'/'kernel_max' \n")
        input <- readline("Select one (y) - If you would like to proceed with a user-specified bandwidth,\n (n) - Exit and change bandwidth search parameters : ")
        if (input == "y") {
            sig_opt <- readline("Please specify the desired kernel bandwidth (Warning : This will directly determine # of subpopulations):")
            sig_opt <- as.numeric(sig_opt)
        }
    }
    
    ## Compute optimal density map
    if (sig_id > 0) {
        opt_densityIMG <- peak_list[[sig_id]]
    } else {
        opt_densityIMG <- peakFind(ydata, sig_opt)
    }
    
    clusters <- assign_pop(opt_densityIMG, ydata, xdata)
    
    return(list(peakdata = peakdata, clusters= clusters))
}

#' assign cells into clusters 
#'
#' @importFrom e1071 svm 
#' @noRd
assign_pop <- function(opt_densityIMG, ydata, xdata) {
    
    density_output <- opt_densityIMG$denImg
    p <- opt_densityIMG$peakCoords
    p_num <- opt_densityIMG$Npeaks
    
    ## assign subpopulations
    Subpop <- array(list(), p_num)
    ind_vector <- vector()
    cluster_vector <- vector()
    
    for (i in 1:p_num) {
        # Find closest subpopulation and determine appropriate radius
        # of 2-d space to draw sample from
        dists <- (kronecker(matrix(1, p_num - 1, 1), matrix(p[i, 
            ], 1, 2)) - p[-i, ])^2
        dists <- matrix(apply(dists, 1, sum), p_num - 1, 1)
        min_dist <- sqrt(min(dists))
        dist2 <- ydata - kronecker(matrix(1, dim(ydata)[1], 1), 
            matrix(p[i, ], 1, 2))
        dist2 <- matrix(apply(dist2^2, 1, sum), dim(ydata)[1], 
            1)
        ind <- which(dist2 < (min_dist/2)^2)
        ind_vector <- c(ind_vector, ind)
        cluster_vector <- c(cluster_vector, rep(i, length(ind)))
    }
    
    ## SVM prediction using the density cluster results
    train_data <- xdata[ind_vector, ]
    train_class <- cluster_vector
    svm.obj <- svm(train_data, train_class, type = "C-classification")
    pred <- predict(svm.obj, xdata)
    cluster <- data.frame(cluster = pred)
    
    ## check and merge the cluster results
    if (nrow(xdata) != nrow(cluster)) 
        print("Cluster Error!.\n")
    clusterResults <- data.frame(ydata, cluster)
    
    return(clusterResults)
}

peakFind <- function(ydata, sig_tol) {
    cat("Computing number of peaks for kernel bandwidth = ", 
        sig_tol, "\n")
    img <- densityIMG(ydata, sig_tol)
    
    d <- img[[1]]
    edg <- 6
    threshold <- max(c(min(apply(d, 1, max)), min(apply(t(d), 
        1, max))))
    
    # Apply threshold
    indicator <- (d > threshold) * 1  # Indicator matrix
    d <- d * indicator
    
    if (sum(d) != 0) {
        # If the image is still non-zero Peak Find - Using the local
        # maxima approach Skip the edge pixels
        
        rows_cols <- which(d[edg:(dim(d)[1] - edg), edg:(dim(d)[2] - 
            edg)] > 0, arr.ind = TRUE)
        x <- rows_cols[, 1]
        y <- rows_cols[, 2]
        cent <- c(0, 0)  #Initialize
        
        # Initialize outputs
        x <- x + edg - 1
        y <- y + edg - 1
        
        for (j in 1:length(y)) {
            
            if (d[x[j], y[j]] >= d[x[j] - 1, y[j] - 1] && d[x[j], 
                y[j]] > d[x[j] - 1, y[j]] && d[x[j], y[j]] >= 
                d[x[j] - 1, y[j] + 1] && d[x[j], y[j]] > d[x[j], 
                y[j] - 1] && d[x[j], y[j]] > d[x[j], y[j] + 1] && 
                d[x[j], y[j]] >= d[x[j] + 1, y[j] - 1] && d[x[j], 
                y[j]] > d[x[j] + 1, y[j]] && d[x[j], y[j]] >= 
                d[x[j] + 1, y[j] + 1]) {
                cent <- rbind(cent, c(img$x[x[j], y[j]], img$y[x[j], 
                  y[j]]))
            }
        }
    }
    
    ## get peak nums and return peaks with out the first row which
    ## is a dummy
    peakNums <- dim(cent)[1] - 1
    cent <- cent[-1, ]
    return(list(Npeaks = peakNums, peakCoords = cent, denImg = img))
}

densityIMG <- function(ydata, sig_tol) {
    ## Partitions the y-space into a 400 x 400 pixel grid and
    ## estimate local density
    density_output <- kernalDensity(ydata, 400, 400, sig_tol)
    ## normize the density by the sum
    density_output[[1]] <- density_output[[1]]/sum(density_output[[1]])
    ## remove potential noise
    density_output[[1]][density_output[[1]] < 5e-08] <- 0
    
    return(density_output)
}

kernalDensity <- function(ydata, width, height, sig_tol) {
    
    # Returns a density image of the data ydata - M x 2 vector
    # containing the y1-y2 coordinates of map-points width,
    # height - width and height of the map area in pixels sig_tol
    # - kernel density function variance (equal for x and y)
    
    limits <- matrix(0, 1, 4)
    limits[1, 1] <- min(ydata[, 1]) - 5
    limits[1, 2] <- max(ydata[, 1]) + 5
    limits[1, 3] <- min(ydata[, 2]) - 5
    limits[1, 4] <- max(ydata[, 2]) + 5
    
    deltax <- (limits[2] - limits[1])/width
    deltay <- (limits[4] - limits[3])/height
    
    dmap <- matrix(0, height, width)
    
    ## kernal density transform
    for (i in 0:(height - 1)) {
        yi <- limits[3] + i * deltay + deltay/2
        for (j in 0:(width - 1)) {
            xi <- limits[1] + j * deltax + deltax/2
            dist2 <- (ydata[, 1] - xi)^2 + (ydata[, 2] - yi)^2
            dd <- sum(exp(-dist2/(2 * sig_tol^2)))
            dmap[i + 1, j + 1] <- (1/sqrt(2 * pi * sig_tol^2)) * 
                dd
        }
    }
    
    y2_range <- matrix(0, height, 1)
    for (i in 0:(height - 1)) {
        y2_range[i + 1] <- limits[3] + i * deltay + deltay/2
    }
    y2_range <- kronecker(matrix(1, 1, dim(dmap)[2]), y2_range)
    
    y1_range <- matrix(0, 1, width)
    for (i in 0:(width - 1)) {
        y1_range[i + 1] <- limits[1] + i * deltax + deltax/2
    }
    y1_range <- kronecker(matrix(1, dim(dmap)[1], 1), y1_range)
    
    return(list(z = dmap, x = y1_range, y = y2_range))
} 
