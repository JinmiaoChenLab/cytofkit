#' Save the cytofkit analysis results
#' 
#' Scatter dot plot and heatmap of the cluster results, and all intermediate files will be 
#' generated and saved in the \code{resDir}
#'
#' @param analysis_results result data from output of \code{\link{densVM_cluster}}
#' @param vizMethods visualization methods for clustering results, including \code{tsne}, \code{pca} and \code{isomap}.
#' @param baseName a prefix that will be added to the names of result files.
#' @param rawFCSdir the directory that contains fcs files to be analysed.
#' @param resDir the directory where result files will be generated.
#' @return save all results in the \code{resDir}
#' @importFrom gplots heatmap.2 bluered    
#' @importFrom ggplot2 ggplot ggsave aes_string facet_wrap geom_point geom_rug theme_bw theme xlab ylab ggtitle coord_fixed guides guide_legend scale_shape_manual scale_colour_manual
#' @importFrom reshape cast melt.data.frame
#' @importFrom ggplot2 ggplot ggsave aes_string geom_line geom_point xlab ylab ggtitle theme_bw
#' @importFrom flowCore write.FCS flowFrame inverseLogicleTransform
#' @export
#' @seealso \code{\link{cytof_tsne_densvm}}, \code{\link{cytofkit}}
#' @examples
#' dir <- system.file('extdata',package='cytofkit')
#' f <- list.files(dir, pattern='.fcs$', full=TRUE)
#' p <- list.files(dir, pattern='.txt$', full=TRUE)
#' #tr <- cytof_tsne_densvm(fcsFile=f,paraFile=p,baseName='t',writeResults=FALSE)
#' #cytof_write_results(tr,baseName = 'test',rawFCSdir=dir)

cytof_write_results <- function(analysis_results, vizMethods, baseName = "cytofkit_analysis", 
    rawFCSdir = getwd(), resDir = getwd()) {
    setwd(resDir)
    ## exprs
    exprs <- analysis_results$lgclMergedExprs
    write.table(exprs, paste(baseName, "_lgcl_merged_markerFiltered_exprsData.txt", 
        sep = ""), sep = "\t", col.names = NA)
    ## transformed
    transformed <- analysis_results$transData
    write.table(transformed, paste(baseName, "_dimension_transformed.txt", 
        sep = ""), sep = "\t", col.names = NA)
    
    ## clusters
    if (!(is.null(analysis_results$clustersRes))) {
        peakData <- analysis_results$clustersRes[[1]]
        peaksGamma_graph <- peaksGamma_plot(peakData)
        ggsave(filename = paste(baseName, "NumOfpeaks_Vs_kernalBandwith.pdf", 
            sep = "_"), plot = peaksGamma_graph)
        
        clusters <- analysis_results$clustersRes[[2]]
        ifMultiFCS <- length(unique(sub("_[0-9]*$", "", row.names(clusters))))>1
        
        ## cluster heatmap
        exprs_cluster <- data.frame(exprs, cluster = clusters[, 3][match(rownames(clusters), rownames(exprs))])
        clust_statData <- clust_state(exprs_cluster, stat = "mean")
        write.table(clust_statData[[1]], paste(baseName, "clusterMean.txt", 
                                               sep = "_"), sep = "\t", col.names = NA)
        write.table(clust_statData[[2]], paste(baseName, "_clusterCellCount.txt", 
                                               sep = ""), sep = "\t", row.names = FALSE)
        
        pdf(paste(baseName, "clusterMeanHeatmap.pdf", sep = "_"))
        par(cex.main = 0.8)
        clust_mean_heatmap(clust_statData[[1]], baseName)
        dev.off()
        
        if (ifMultiFCS) {
                pdf(paste(baseName, "clusterPercHeatmap.pdf", sep = "_"))
                par(mar = rep(2, 4), cex.main = 0.7)
                clust_percentage_heatmap(clust_statData[[2]], baseName)
                dev.off()
        }
        
        ## data to add to fcsfiles
        to_add <- clusters
        trans_col_names <- setdiff(colnames(to_add), "cluster")
        
        ## cluster scatter plot
        for(x in sort(vizMethods, decreasing = TRUE)){
                if (x == "tsne"){
                        clusters <- analysis_results$clustersRes[[2]]
                } else{
                        transformed <- cytof_dimReduction(exprs, method = x)
                        clusters <- merge(transformed[ ,c(1,2)], subset(clusters, select = get("cluster")), by = "row.names")
                        rownames(clusters) <- clusters[ ,1]
                        clusters <- clusters[ ,-1]  
                        ## update to_add data
                        if(ncol(transformed) > 2){
                                s_cols <- c(1,2,3)  ## maxi 3 columns will be saved for each method
                        }else{
                                s_cols <- c(1 : ncol(transformed))
                        }
                        to_add <- merge(to_add, transformed[ ,s_cols], by = "row.names")
                        rownames(to_add) <- to_add[ ,1]
                        to_add <- to_add[ ,-1]  
                        trans_col_names <- setdiff(colnames(to_add), "cluster")
                }
                
                ## visualize clusters and save resutls
                write.table(clusters, paste(baseName, x, "cluster.txt", 
                                            sep = "_"), sep = "\t", col.names = NA)
                
                cluster_plot <- cluster_plot(clusters)
                ggsave(filename = paste(baseName, x, "clusterData_plot.pdf", 
                                        sep = "_"), cluster_plot, width = 13, height = 13)
                if (ifMultiFCS) {
                        samples <- sub("_[0-9]*$", "", row.names(clusters))
                        sample_num <- length(unique(samples))
                        grid_row_num <- ceiling(sample_num/4)
                        if (sample_num >= 4) {
                                grid_col_num <- 4
                        } else {
                                grid_col_num <- sample_num
                        }
                        grid_size <- sqrt(225/(grid_row_num * grid_col_num))
                        grid_width <- grid_size * grid_col_num
                        grid_height <- grid_size * grid_row_num
                        cluster_grid_plot <- cluster_gridPlot(clusters)
                        ggsave(filename = paste(baseName, x, "clusterData_grid_plot.pdf", 
                                                sep = "_"), cluster_grid_plot, width = grid_width, 
                               height = grid_height)
                }      
        }
        
        suppressWarnings(add_col_to_fcs(data = to_add, rawFCSdir = rawFCSdir, 
                                        analyzedFCSdir = paste(baseName, "analyzedFCS", sep = "_"), 
                                        transformed_col = trans_col_names, cluster_col = c("cluster")))     
    } else {
        suppressWarnings(add_col_to_fcs(data = transformed, rawFCSdir = rawFCSdir, 
            analyzedFCSdir = "analyzedFCS", transformed_col = trans_col_names, 
            cluster_col = NULL))
    }
}

#' Plot varaition of peak nums with increasing gamma
#' 
#' @param peakdata a matrix of \code{peakdata} returned from \code{densVM_cluster}
#' @return a line graph of peak nums vs. increasing gamma
#' @export
#' @examples
#' x <- seq(0, 1, length.out = 20)
#' y <- c(20:6, 6, 6, 5:3)
#' peakdata <- data.frame(sig_range = x, numpeaks = y)
#' peaksGamma_plot(peakdata)
peaksGamma_plot <- function(peakdata) {
    ggplot(data = peakdata, aes_string(x = "sig_range", y = "numpeaks")) + 
        geom_line() + geom_point() + theme_bw() + xlab(expression(gamma)) + 
        ylab("N_peaks") + ggtitle("Number of peaks with increasing bandwith(gamma)")
}


#' Scatter plot of the cluster results
#' 
#' Dot plot visualization of the cluster results, with color indicating different clusters, 
#' and shape of different samples.
#' 
#' @param clusterData The matrix of cluster results, with rownames of their sample name and cell id.
#' @param title the title name of the plot.
#' @param point_size the size of the dot.
#' @return the scatter dot plot
#' @export
#' @examples
#' x <- c(rnorm(100, mean = 1), rnorm(100, mean = 3), rnorm(100, mean = 9))
#' y <- c(rnorm(100, mean = 2), rnorm(100, mean = 8), rnorm(100, mean = 5))
#' c <- c(rep(1,100), rep(2,100), rep(3,100))
#' rnames <- paste(paste('sample_', c('A','B','C'), sep = ''), rep(1:100,each = 3), sep='_') 
#' clusterData <- data.frame(dim1 = x, dim2 = y, cluster = c)
#' rownames(clusterData) <- rnames
#' cluster_plot(clusterData)
cluster_plot <- function(clusterData, title = "cluster", point_size = NULL) {
    
    # plot the cluster results with cluster labels
    clusterData$sample <- sub("_[0-9]*$", "", row.names(clusterData))
    clusterData$cluster <- as.factor(clusterData$cluster)
    lab1 <- colnames(clusterData)[1]
    lab2 <- colnames(clusterData)[2]
    cluster <- colnames(clusterData)[3]
    sample <- colnames(clusterData)[4]
    
    cluster_num <- length(unique(clusterData$cluster))
    sample_num <- length(unique(clusterData$sample))
    col_legend_row <- ceiling(cluster_num/15)
    size_legend_row <- ceiling(sample_num/4)
    if (sample_num >= 4) {
        grid_col_num <- 4
        shape_value <- LETTERS[1:sample_num]
    } else {
        grid_col_num <- sample_num
        shape_value <- c(1:sample_num) + 15
    }
    if (is.null(point_size)) {
        point_size <- ifelse(nrow(clusterData) > 10000, 1, 1.5)
    }
    
    ggplot(clusterData, aes_string(x = lab1, y = lab2, colour = cluster, 
        shape = sample)) + geom_point(size = point_size) + geom_rug(position = "jitter", 
        size = 0.2) + theme_bw() + scale_shape_manual(values = shape_value) + 
        scale_colour_manual(values = rainbow(cluster_num)) + 
        theme(legend.position = "bottom") + xlab(lab1) + ylab(lab2) + 
        ggtitle(paste(title, "scatter plot", sep = " ")) + coord_fixed() + 
        guides(colour = guide_legend(nrow = col_legend_row), 
            shape = guide_legend(nrow = size_legend_row))
}


#' Grid scatter plot of the cluster results with multiple samples
#' 
#' Grid dot plot visualization of the cluster results, with color indicating 
#' different clusters, and panels of different samples.
#' 
#' @param clusterData The matrix of cluster results, with rownames of their sample name and cell id.
#' @param title the title name of the plot.
#' @param point_size the size of the dot.
#' @return the grid scatter dot plot
#' @export
#' @examples
#' x <- c(rnorm(100, mean = 1), rnorm(100, mean = 3), rnorm(100, mean = 9))
#' y <- c(rnorm(100, mean = 2), rnorm(100, mean = 8), rnorm(100, mean = 5))
#' c <- c(rep(1,100), rep(2,100), rep(3,100))
#' rnames <- paste(paste('sample_', c('A','B','C'), sep = ''), rep(1:100,each = 3), sep='_') 
#' clusterData <- data.frame(dim1 = x, dim2 = y, cluster = c)
#' rownames(clusterData) <- rnames
#' cluster_gridPlot(clusterData)
cluster_gridPlot <- function(clusterData, title = "cluster", 
    point_size = NULL) {
    
    # plot the cluster results with cluster labels
    clusterData$sample <- sub("_[0-9]*$", "", row.names(clusterData))
    m <- regexpr("_[0-9]*$", row.names(clusterData))
    clusterData$sample_ID <- substring(regmatches(row.names(clusterData), 
        m), 2)
    clusterData$cluster <- as.factor(clusterData$cluster)
    
    lab1 <- colnames(clusterData)[1]
    lab2 <- colnames(clusterData)[2]
    cluster <- colnames(clusterData)[3]
    sample <- colnames(clusterData)[4]
    sample_ID <- colnames(clusterData)[5]
    
    cluster_num <- length(unique(clusterData$cluster))
    sample_num <- length(unique(clusterData$sample))
    col_legend_row <- ceiling(cluster_num/15)
    size_legend_row <- ceiling(sample_num/4)
    grid_row_num <- size_legend_row
    if (sample_num >= 4) {
        grid_col_num <- 4
        shape_value <- LETTERS[1:sample_num]
    } else {
        grid_col_num <- sample_num
        shape_value <- c(1:sample_num) + 15
    }
    grid_size <- sqrt(225/(grid_row_num * grid_col_num))
    grid_width <- grid_size * grid_col_num
    grid_height <- grid_size * grid_row_num
    if (is.null(point_size)) {
        point_size <- ifelse(nrow(clusterData) > 10000, 1, 1.5)
    }
    
    if (sample_num > 1) {
        ggplot(clusterData, aes_string(x = lab1, y = lab2, colour = cluster)) + 
            facet_wrap(~sample, ncol = 4, scales = "free") + 
            geom_point(size = point_size - 0.05 * sample_num) + 
            geom_rug(position = "jitter", size = 0.2) + theme_bw() + 
            scale_colour_manual(values = rainbow(cluster_num)) + 
            theme(legend.position = "bottom") + xlab(lab1) + 
            ylab(lab2) + ggtitle(paste(title, "grid plot", sep = " ")) + 
            coord_fixed() + guides(colour = guide_legend(nrow = col_legend_row), 
            shape = guide_legend(nrow = size_legend_row))
    } else {
        stop("No multiple samples in the data")
    }
}


#' Heatmap plot of cluster mean value results
#' 
#' @param baseName The name as a prefix in the title of the heatmap.
#' @param clust_mean cluster mean data from results of \code{clust_state}.
#' @param scaleMethod character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. The default is "none".
#' @return a heatmap object from \code{gplots}
#' @export
#' @examples
#' m1 <- c(rnorm(300, 10, 2), rnorm(400, 4, 2), rnorm(300, 7))
#' m2 <- c(rnorm(300, 4), rnorm(400, 16), rnorm(300, 10, 3))
#' m3 <- c(rnorm(300, 16), rnorm(400, 40, 3), rnorm(300, 10))
#' m4 <- c(rnorm(300, 7, 3), rnorm(400, 30, 2), rnorm(300, 10))
#' m5 <- c(rnorm(300, 27), rnorm(400, 40, 1),rnorm(300, 10))
#' c <- c(rep(1,300), rep(2,400), rep(3,300))
#' rnames <- paste(paste('sample_', c('A','B','C','D'), sep = ''), 
#' rep(1:250,each = 4), sep='_') 
#' exprs_cluster <- data.frame(cluster = c, m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5)
#' row.names(exprs_cluster) <- sample(rnames, 1000)
#' clust_statData <- clust_state(exprs_cluster)
#' clust_mean_heatmap(clust_statData[[1]])
clust_mean_heatmap <- function(clust_mean, baseName = "Cluster_mean", scaleMethod = "none") {
    
    rownames(clust_mean) <- paste("cluster_", clust_mean$cluster, 
        sep = "")
    clust_mean <- clust_mean[, -which(colnames(clust_mean) == 
        "cluster")]
    
    cex_row_label <- (9 - ceiling(nrow(clust_mean)/10))/10
    cex_col_label <- (9 - ceiling(ncol(clust_mean)/10))/10
    
    heatmap.2(as.matrix(clust_mean), col = bluered, trace = "none", 
        symbreaks = FALSE, scale = scaleMethod, cexRow = cex_row_label, 
        cexCol = cex_col_label, srtCol = 30, symkey = FALSE, 
        main = paste(baseName, "\n clusterMean heatmap", sep = " "))
}


#' Heatmap plot of cluster percentage results
#' 
#' @param baseName The name as a prefix in the title of the heatmap
#' @param clust_cellCount cluster count data from results of \code{clust_state}
#' @param scaleMethod character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. The default is "none".
#' @return a heatmap object from \code{gplots}
#' @export
#' @examples
#' m1 <- c(rnorm(300, 10, 2), rnorm(400, 4, 2), rnorm(300, 7))
#' m2 <- c(rnorm(300, 4), rnorm(400, 16), rnorm(300, 10, 3))
#' m3 <- c(rnorm(300, 16), rnorm(400, 40, 3), rnorm(300, 10))
#' m4 <- c(rnorm(300, 7, 3), rnorm(400, 30, 2), rnorm(300, 10))
#' m5 <- c(rnorm(300, 27), rnorm(400, 40, 1),rnorm(300, 10))
#' c <- c(rep(1,300), rep(2,400), rep(3,300))
#' rnames <- paste(paste('sample_', c('A','B','C','D'), sep = ''), 
#' rep(1:250,each = 4), sep='_') 
#' exprs_cluster <- data.frame(cluster = c, m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5)
#' row.names(exprs_cluster) <- sample(rnames, 1000)
#' clust_statData <- clust_state(exprs_cluster)
#' clust_percentage_heatmap(clust_statData[[2]])
clust_percentage_heatmap <- function(clust_cellCount, baseName = "Cluster percentage", scaleMethod = "none") {
    
    clust_percentage <- cast(clust_cellCount, sample ~ cluster, 
        value = "percentage")
    row.names(clust_percentage) <- clust_percentage$sample
    clust_percentage <- clust_percentage[, -which(colnames(clust_percentage) == 
        "sample")]
    colnames(clust_percentage) <- paste("cluster_", colnames(clust_percentage), 
        sep = "")
    
    cex_row_label <- (9 - ceiling(ncol(clust_percentage)/10))/10
    cex_col_label <- (8 - ceiling(nrow(clust_percentage)/10))/10
    
    heatmap.2(t(as.matrix(clust_percentage)), col = bluered, 
        trace = "none", scale = scaleMethod, labRow = colnames(clust_percentage), 
        labCol = rownames(clust_percentage), cexRow = cex_row_label, 
        cexCol = cex_col_label, srtCol = 30, symkey = FALSE, 
        symbreaks = FALSE, main = paste(baseName, "\n clusterPerc heatmap", 
            sep = " "))
}

#' Statistical analysis of the cluster results
#' 
#' Calculate the mean value of each markers in each cluster, If there are multiple samples, the percentage
#' of cells in each cluster in each sample will be calculated
#' 
#' @param exprs_cluster the expression matrix combined with the cluster results
#' @param stat the method used for statistical analysis, like mean, median...
#' @return a list contains a matrix of \code{clust_mean} and a matrix of \code{clust_cellCount}
#' @export
#' @examples
#' m1 <- c(rnorm(300, 10, 2), rnorm(400, 4, 2), rnorm(300, 7))
#' m2 <- c(rnorm(300, 4), rnorm(400, 16), rnorm(300, 10, 3))
#' m3 <- c(rnorm(300, 16), rnorm(400, 40, 3), rnorm(300, 10))
#' m4 <- c(rnorm(300, 7, 3), rnorm(400, 30, 2), rnorm(300, 10))
#' m5 <- c(rnorm(300, 27), rnorm(400, 40, 1),rnorm(300, 10))
#' c <- c(rep(1,300), rep(2,400), rep(3,300))
#' rnames <- paste(paste('sample_', c('A','B','C','D'), sep = ''), 
#' rep(1:250,each = 4), sep='_') 
#' exprs_cluster <- data.frame(cluster = c, m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5)
#' row.names(exprs_cluster) <- sample(rnames, 1000)
#' clust_statData <- clust_state(exprs_cluster)
clust_state <- function(exprs_cluster, stat = "mean") {
    
    # calculate cluster mean
    clust_mean <- aggregate(. ~ cluster, data = exprs_cluster, 
        stat)
    
    # calculate cluster percentage
    exprs_cluster$sample <- sub("_[0-9]*$", "", row.names(exprs_cluster))
    clust_cellCount <- as.data.frame(table(exprs_cluster[, c("sample", 
        "cluster")]))
    colnames(clust_cellCount)[3] <- "cellCount"
    sample_cellCount <- as.data.frame(table(exprs_cluster$sample))
    colnames(sample_cellCount) <- c("sample", "totalCellCount")
    clust_cellCount <- merge(clust_cellCount, sample_cellCount, 
        by = "sample")
    clust_cellCount$percentage <- clust_cellCount$cellCount/clust_cellCount$totalCellCount * 
        100
    
    return(list(clust_mean = clust_mean, clust_cellCount = clust_cellCount))
}

#' Add new columns to the fcs expression data
#' 
#' Store the new dimension transformed data and cluster data into the exprs 
#' matrix in new fcs files under \code{analyzedFCSdir}
#' 
#' @param data The new data matrix to be added in.
#' @param rawFCSdir The directory containing the original fcs files.
#' @param analyzedFCSdir The directory to store the new fcs files.
#' @param transformed_col the column name of the dimension transformend data in \code{data}.
#' @param cluster_col the column name of the cluster data in \code{data}.
#' @return new fcs files stored under \code{analyzedFCSdir}
#' @importMethodsFrom flowCore keyword
#' @importFrom Biobase  exprs exprs<- description description<- pData pData<- 
add_col_to_fcs <- function(data, rawFCSdir, analyzedFCSdir, 
    transformed_col = c("tsne_1", "tsne_2"), cluster_col = c("cluster")) {
    
    lgcl <- logicleTransform(w = 0.1, t = 4000, m = 4.5, a = 0)
    ilgcl <- inverseLogicleTransform(trans = lgcl)
    
    if (!file.exists(analyzedFCSdir)) {
        dir.create(analyzedFCSdir)
    }
    
    if (!is.null(transformed_col)) {
        transformed <- data[, transformed_col]
        row.has.na <- !(complete.cases(transformed))
        transformed <- transformed[!row.has.na, ]
#         N_transformed <- apply(transformed, 2, function(x) ((x-min(x))/(max(x)-min(x)))*4.4)
#         R_N_transformed <- apply(N_transformed,2,ilgcl)
        R_N_transformed <- apply(transformed, 2, function(x) (x-min(x)) + 0.1)
        row.names(R_N_transformed) <- row.names(transformed)
    }
    
    ## transform cluster column
    if (!is.null(cluster_col)) {
        if (row.has.na) {
                clust_cor_1 <- as.numeric(data[!row.has.na, cluster_col])%%10
                clust_cor_2 <- floor(as.numeric(data[!row.has.na, cluster_col])/10)
        } else {
            clust_cor_1 <- as.numeric(data[, cluster_col])%%10
            clust_cor_2 <- floor(as.numeric(data[, cluster_col])/10)
        }
        clust_cor_1 <- clust_cor_1 + runif(length(clust_cor_1), 
            0, 0.2)
        clust_cor_2 <- clust_cor_2 + runif(length(clust_cor_2), 
            0, 0.2)
        clust_cor <- cbind(clust_cor_1, clust_cor_2)
        N_clust_cor <- apply(clust_cor, 2, function(x) ((x - 
            min(x))/(max(x) - min(x))) * 4.4)
        R_N_clust_cor <- apply(N_clust_cor, 2, ilgcl)
        if (row.has.na) {
            row.names(R_N_clust_cor) <- row.names(data)[!row.has.na]
        } else {
            row.names(R_N_clust_cor) <- row.names(data)
        }
    }
    
    if ((!is.null(transformed_col)) && (!is.null(cluster_col))) {
        if (sum(row.names(R_N_transformed) != row.names(R_N_clust_cor)) == 
            0) {
            to_add <- cbind(R_N_transformed, R_N_clust_cor)
        } else {
            print("error:the row.names is not consistent between transformed coordinates and cluster.\n")
        }
    } else if (!is.null(transformed_col)) {
        to_add <- R_N_transformed
    } else if (!is.null(cluster_col)) {
        to_add <- R_N_clust_cor
    }
    addColNames <- colnames(to_add)
    
    sample <- unique(sub("_[0-9]*$", "", row.names(to_add)))
    for (i in 1:length(sample)) {
        
        fcs <- read.FCS(paste(rawFCSdir, "/", sample[i], ".fcs", 
            sep = ""))
        pattern <- paste(sample[i], "_", sep = "")
        to_add_i <- as.data.frame(to_add[grep(pattern, row.names(to_add), 
            fixed = TRUE), ])
        m <- regexpr("_[0-9]*$", row.names(to_add_i))
        cellNo_i <- as.integer(substring(regmatches(row.names(to_add_i), 
            m), 2))
        sub_exprs <- fcs@exprs[cellNo_i, ]
        params <- parameters(fcs)
        pd <- pData(params)
        keyval <- keyword(fcs)
        
        for (j in 1:length(addColNames)) {
            ## update the parameters
            if (addColNames[j] %in% colnames(fcs@exprs)) {
                addColNames[j] <- paste(addColNames[j], "_new", 
                  sep = "")
            }
            
            addColName <- addColNames[j]
            channel_number <- nrow(pd) + 1
            channel_id <- paste("$P", channel_number, sep = "")
            channel_name <- addColName
            minRange <- ceiling(min(to_add_i[[addColName]]))
            maxRange <- ceiling(max(to_add_i[[addColName]]))
            channel_range <- maxRange - minRange
            
            plist <- matrix(c(channel_name, "<NA>", channel_range, 
                minRange, maxRange))
            rownames(plist) <- c("name", "desc", "range", "minRange", 
                "maxRange")
            colnames(plist) <- c(channel_id)
            pd <- rbind(pd, t(plist))
            
            ## update the expression value
            out_col_names <- colnames(sub_exprs)
            sub_exprs <- cbind(sub_exprs, to_add_i[[addColName]])
            colnames(sub_exprs) <- c(out_col_names, addColName)
            
            ## update the description remove '\' in the keywords
            keyval <- lapply(keyval, function(x) {
                if (class(x) == "character") {
                  gsub("\\", "", x, fixed = TRUE)
                } else {
                  x
                }
            })
            keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"  # Number of bits
            keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)  # Range
            keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"  # Exponent
            keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name  # Name
            keyval[[paste("P", channel_number, "BS", sep = "")]] <- 0
            keyval[[paste("P", channel_number, "MS", sep = "")]] <- 0
            keyval[[paste("P", channel_number, "DISPLAY", sep = "")]] <- "LIN"  # data display
        }
        
        pData(params) <- pd
        out_frame <- flowFrame(exprs = sub_exprs, parameters = params, 
            description = keyval)
        
        write.FCS(out_frame, paste(analyzedFCSdir, "/", sample[i], 
            ".fcs", sep = ""))
    }
} 
