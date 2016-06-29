#' Save the cytofkit analysis results
#' 
#' Scatter dot plot and heatmap of the cluster results, and all intermediate files will be 
#' generated and saved in the \code{resultDir}
#'
#' @param analysis_results result data from output of \code{\link{cytofkit}}
#' @param projectName a prefix that will be added to the names of result files.
#' @param resultDir the directory where result files will be generated.
#' @param saveToFCS save the results back to the FCS files, new FCS files will be generated.
#' @param rawFCSdir the directory that contains fcs files to be analysed.
#' @return save all results in the \code{resultDir}
#' @importFrom ggplot2 ggplot ggsave aes_string facet_wrap geom_point geom_rug theme_bw theme xlab ylab ggtitle coord_fixed guides guide_legend scale_shape_manual scale_colour_manual
#' @importFrom reshape2 dcast
#' @importFrom ggplot2 ggplot ggsave aes_string geom_line geom_point xlab ylab ggtitle theme_bw
#' @importFrom flowCore write.FCS flowFrame inverseLogicleTransform
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom graphics par
#' @importFrom utils read.table write.csv
#' @export
#' @seealso \code{\link{cytofkit}}
#' @examples
#' d <- system.file('extdata',package='cytofkit')
#' f <- list.files(d, pattern='.fcs$', full=TRUE)
#' p <- list.files(d, pattern='.txt$', full=TRUE)
#' #tr <- cytofkit(fcsFile=f,markers=p,projectName='t',saveResults=FALSE)
#' #cytof_write_results(tr,projectName = 'test',resultDir=d,rawFCSdir =d)
cytof_writeResults <- function(analysis_results, 
                               projectName, 
                               resultDir, 
                               saveToFCS = TRUE, 
                               rawFCSdir = getwd()) {
    
    if(missing(projectName)){
        if(!is.null(analysis_results$projectName)){
            projectName <- analysis_results$projectName
        }else{
            projectName <- "cytofkit_"
        }
    }
    
    if(missing(resultDir)){
        if(!is.null(analysis_results$resultDir)){
            resultDir <- analysis_results$resultDir
        }else{
            resultDir <- getwd()
        }
    }
    
    if(missing(rawFCSdir)){
        if(!is.null(analysis_results$rawFCSdir)){
            rawFCSdir <- analysis_results$rawFCSdir
        }else{
            rawFCSdir <- getwd()
        }
    }
     
    curwd <- getwd()
    setwd(resultDir)
    
    ## save exprs
    exprs <- as.data.frame(analysis_results$expressionData)
    ifMultiFCS <- length(unique(sub("_[0-9]*$", "", row.names(exprs)))) > 1
    write.csv(exprs, paste0(projectName, "_markerFiltered_transformed_merged_exprssion_data.csv"))
    
    ## save dimReducedData
    dimReducedData <- analysis_results$dimReducedRes
    for(i in 1:length(dimReducedData)){
        methodi <- names(dimReducedData)[i]
        if(!is.null(dimReducedData[[i]])){
            write.csv(dimReducedData[[i]], paste(projectName, methodi,"dimension_reduced_data.csv", sep="_"))
        }
    }
    
    ## save clusterData
    clusterData <- analysis_results$clusterRes
    if(!is.null(clusterData) && length(clusterData) > 0){
        for(j in 1:length(clusterData)){
            methodj <- names(clusterData)[j]
            dataj <- clusterData[[j]]
            if(!is.null(dataj)){
                write.csv(dataj, paste(projectName, methodj, "clusters.csv", sep="_"))
                
                exprs_cluster_sample <- data.frame(exprs, cluster = dataj, check.names = FALSE)
                ## cluster mean 
                cluster_mean <- cytof_clusterStat(data= exprs_cluster_sample, cluster = "cluster", statMethod = "mean")
                write.csv(cluster_mean, paste(projectName, methodj, "cluster_mean_data.csv", sep = "_"))
                pdf(paste(projectName, methodj, "cluster_mean_heatmap.pdf", sep = "_"))
                cytof_heatmap(cluster_mean, paste(projectName, methodj, "\ncluster mean", sep = " "))
                dev.off()
                ## cluster median
                cluster_median <- cytof_clusterStat(data= exprs_cluster_sample, cluster = "cluster", statMethod = "median")
                write.csv(cluster_median, paste(projectName, methodj, "cluster_median_data.csv", sep = "_"))
                pdf(paste(projectName, methodj, "cluster_median_heatmap.pdf", sep = "_"))
                cytof_heatmap(cluster_median, paste(projectName, methodj, "\ncluster median", sep = " "))
                dev.off()
                
                ## cluster percentage
                if (ifMultiFCS) {
                    cluster_percentage <- cytof_clusterStat(data= exprs_cluster_sample, cluster = "cluster", statMethod = "percentage")
                    write.csv(cluster_percentage, paste(projectName, methodj, "cluster_cell_percentage.csv", sep = "_"))
                    pdf(paste(projectName, methodj, "cluster_percentage_heatmap.pdf", sep = "_"))
                    cytof_heatmap(cluster_percentage, paste(projectName, methodj, "cluster\ncell percentage", sep = " "))
                    dev.off()
                }
            }
        }
    }  
    
    ## visualizationData x clusterData plot
    visualizationData <- analysis_results$dimReducedRes[analysis_results$visualizationMethods]
    for(i in 1:length(visualizationData)){
        if(!is.null(visualizationData[[i]])){
            methodi <- names(visualizationData)[i]
            datai <- as.data.frame(visualizationData[[i]])
            if(!is.null(clusterData) && length(clusterData) > 0){
                for(j in 1:length(clusterData)){
                    if(!is.null(clusterData[[j]])){
                        methodj <- names(clusterData)[j]
                        dataj <- clusterData[[j]]
                        
                        # combine datai and dataj
                        xlab <- colnames(datai)[1]
                        ylab <- colnames(datai)[2]
                        dataij <- datai
                        dataij$sample <- sub("_[0-9]*$", "", row.names(dataij))
                        dataij$cluster <- factor(dataj)
                        cluster <- "cluster"
                        sample <- "sample"
                        
                        ## cluster plot
                        figName <- paste(projectName, methodi, methodj, sep=" ")
                        cluster_plot <- cytof_clusterPlot(dataij, xlab, ylab, cluster, sample, figName, 1)
                        ggsave(filename = paste(projectName, methodi, methodj, "cluster_scatter_plot.pdf", sep = "_"), 
                               cluster_plot, width = 12, height = 10)
                        
                        ## cluster grid plot if multiple files
                        if (ifMultiFCS) {
                            figName <- paste(projectName, methodi, methodj, sep=" ")
                            cluster_grid_plot <- cytof_clusterPlot(dataij, xlab, ylab, cluster, sample, figName, 2)
                            ggsave(filename = paste(projectName, methodi, methodj, "cluster_grid_scatter_plot.pdf", sep = "_"), cluster_grid_plot)
                        }
                    }
                }
            }  
        }
    }
    
    ## save data to FCS
    if(saveToFCS == TRUE){
        tcols <- do.call(cbind, dimReducedData)
        ctols <- do.call(cbind, clusterData)
        dataToAdd <- cbind(tcols, ctols)
        row.names(dataToAdd) <- row.names(exprs)
        trans_col_names <- colnames(tcols)
        cluster_col_names <- colnames(ctols)
        cytof_addToFCS(dataToAdd, rawFCSdir = rawFCSdir, 
                       analyzedFCSdir = paste(projectName, "analyzedFCS", sep = "_"), 
                       transformed_cols = trans_col_names, cluster_cols = cluster_col_names)
    }
    setwd(curwd)
    
    message(paste0("Writing results Done! Results are saved under path: ",
                   resultDir))
    
    return(NULL)
}


#' Scatter plot of the cluster results
#' 
#' Dot plot visualization of the cluster results, with color indicating different clusters, 
#' and shape of different samples.
#' 
#' @param data The data frame of cluster results, which should contains at least xlab, ylab and cluster.
#' @param xlab The column name of the x axis in input \code{data}.
#' @param ylab The column name of the y axis in input \code{data}.
#' @param cluster The column name of cluster in input \code{data}.
#' @param sample the column name of the sample in input \code{data}.
#' @param title the title of the plot.
#' @param type plot type, 1 indicates combined plot, 2 indicated grid facet plot seperated by samples.
#' @param point_size the size of the dot.
#' @param addLabel Boolean, if add cluster labels.
#' @param labelSize the size of cluster labels.
#' @param sampleLabel If use point shapes to represent different samples.
#' @param labelRepel If repel the cluste labels to avoid label overlapping.
#' @param fixCoord If fix the Cartesian coordinates.
#' @return the ggplot object of the scatter cluster plot.
#' @export
#' @importFrom ggplot2 element_text element_rect element_blank element_line element_text annotate
#' @examples
#' x <- c(rnorm(100, mean = 1), rnorm(100, mean = 3), rnorm(100, mean = 9))
#' y <- c(rnorm(100, mean = 2), rnorm(100, mean = 8), rnorm(100, mean = 5))
#' c <- c(rep(1,100), rep(2,100), rep(3,100))
#' rnames <- paste(paste('sample_', c('A','B','C'), sep = ''), rep(1:100,each = 3), sep='_') 
#' data <- data.frame(dim1 = x, dim2 = y, cluster = c)
#' rownames(data) <- rnames
#' data$sample <- "data"
#' cytof_clusterPlot(data, xlab="dim1", ylab="dim2", cluster="cluster", sample = "sample")
cytof_clusterPlot <- function(data, xlab, ylab, cluster, sample, title = "cluster", 
                              type = 1, point_size = NULL, addLabel=TRUE, 
                              labelSize=10, sampleLabel=TRUE, 
                              labelRepel = FALSE, fixCoord=TRUE) {
    
    if(!is.data.frame(data))
        data <- as.data.frame(data)
    
    if(missing(sample)){
        sample <- "sample"
        data$sample <- "data"
    }
    
    paraCheck <- c(xlab, ylab, cluster, sample) %in% colnames(data)
    if(any(!paraCheck)){
        stop("Undefined parameters found: ",
             c(xlab, ylab, cluster, sample)[!paraCheck])
    }
    
    data[[cluster]] <- as.factor(data[[cluster]])
    data[[sample]] <- as.factor(data[[sample]])
    cluster_num <- length(unique(data[[cluster]]))
    sample_num <- length(unique(data[[sample]]))
    col_legend_row <- ceiling(cluster_num/15)
    size_legend_row <- ceiling(sample_num/4)
    grid_col_num <- round(sqrt(sample_num))
    if (sample_num >= 8) {
        shape_value <- LETTERS[1:sample_num]
    } else {
        shape_value <- c(1:sample_num) + 15
    }
    if (is.null(point_size)) {
        point_size <- ifelse(nrow(data) > 10000, 1, 1.5)
    }
    
    if(type == 1){
        if(sampleLabel){
            cp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = cluster, shape = sample)) + 
                geom_point(size = point_size) + scale_shape_manual(values = shape_value) + 
                scale_colour_manual(values = rainbow(cluster_num)) + 
                xlab(xlab) + ylab(ylab) + ggtitle(paste(title, "Scatter Plot", sep = " ")) + 
                theme_bw() + theme(legend.position = "bottom") + 
                theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold")) +
                guides(colour = guide_legend(nrow = col_legend_row, override.aes = list(size = 4)), 
                       shape = guide_legend(nrow = size_legend_row, override.aes = list(size = 4)))
        }else{
            cp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = cluster)) + 
                geom_point(size = point_size) + scale_shape_manual(values = shape_value) + 
                scale_colour_manual(values = rainbow(cluster_num)) + 
                xlab(xlab) + ylab(ylab) + ggtitle(paste(title, "Scatter Plot", sep = " ")) + 
                theme_bw() + theme(legend.position = "bottom") + 
                theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold")) +
                guides(colour = guide_legend(nrow = col_legend_row, override.aes = list(size = 4)))
        }
        
        if(addLabel){
            edata <- data[ ,c(xlab, ylab, cluster)]
            colnames(edata) <- c('x', "y", "z")
            center <- aggregate(cbind(x,y) ~ z, data = edata, median)
            
            if(labelRepel && !sampleLabel){
                cp <- cp + geom_text_repel(data=center, aes_string(x = "x", y = "y", label = "z"),
                                           size = labelSize, fontface = 'bold', color = "black",
                                           box.padding = unit(0.5, 'lines'),
                                           point.padding = unit(1.6, 'lines'),
                                           segment.color = '#555555',
                                           segment.size = 0.5,
                                           arrow = arrow(length = unit(0.02, 'npc')))
            }else{
                cp <- cp + annotate("text", label = center[,1], x=center[,2], y = center[,3],
                                    size = labelSize, colour = "black")
            }
        }
        
    }else if (type == 2){
        cp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = cluster)) + 
            facet_wrap(~sample, ncol = grid_col_num, scales = "fixed") + 
            geom_point(size = point_size - 0.05 * sample_num) + 
            scale_colour_manual(values = rainbow(cluster_num)) + 
            xlab(xlab) + ylab(ylab) + ggtitle(paste(title, "Grid Plot", sep = " ")) + 
            theme_bw() + theme(legend.position = "bottom") + 
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold")) +
            guides(colour = guide_legend(nrow = col_legend_row, override.aes = list(size = 4)), 
                   shape = guide_legend(nrow = size_legend_row, override.aes = list(size = 4)))
    }else{ 
        stop("Undefined type, only 1 or 2.") 
    }
    
    if(fixCoord){
        cp <- cp + coord_fixed()
    }
    
    return(cp)
}



#' Heatmap plot of cluster mean value results
#' 
#' @param data a matrix with rownames and colnames
#' @param baseName The name as a prefix in the title of the heatmap.
#' @param scaleMethod Method indicating if the values should be centered and scaled in either the row direction or the column direction, or none. The default is 'none'.
#' @param cex_row_label Text size for row labels.
#' @param cex_col_label Text size for column labels.
#' @param key.par graphical parameters for the color key. 
#' @param keysize numeric value indicating the size of the key.
#' @param margins numeric vector of length 2 containing the margins (see par(mar= *)) for column and row names, respectively.
#' @return a heatmap object from \code{gplots}
#' 
#' @importFrom gplots heatmap.2 bluered   
#' 
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
#' cluster_mean <- aggregate(. ~ cluster, data = exprs_cluster, mean)
#' rownames(cluster_mean) <- paste("cluster_", cluster_mean$cluster, sep = "")
#' cytof_heatmap(cluster_mean[, -which(colnames(cluster_mean) == "cluster")])
cytof_heatmap <- function(data, baseName = "Cluster", scaleMethod = "none",
                          cex_row_label = NULL, cex_col_label = NULL, 
                          key.par = list(mgp=c(1.5, 0.5, 0), mar=c(2.5, 2.5, 3, 1.5)), 
                          keysize = 1.4,
                          margins = c(5, 5)) {
    
    data <- as.matrix(data)
    
    if(is.null(cex_row_label)){
        cex_row_label <- (11 - ceiling(nrow(data)/10))/10
    }
    if(is.null(cex_col_label)){
        cex_col_label <- (11 - ceiling(ncol(data)/10))/10
    }
    
    heatmap.2(x = data, 
              col = bluered, 
              trace = "none", 
              symbreaks = FALSE, 
              scale = scaleMethod, 
              cexRow = cex_row_label, 
              cexCol = cex_col_label, 
              srtCol = 30, symkey = FALSE, 
              key.par=key.par, 
              margins = margins,
              keysize = keysize,
              main = paste(baseName, "Heat Map"))
}


#' Plot the data with color-coded marker values
#' 
#' @param data A dataframe containing the xlab, ylab and zlab.
#' @param xlab The column name of data for x lab.
#' @param ylab The column name of data for y lab.
#' @param zlab The column name of data for z lab.
#' @param colorPalette Color Palette. 
#' @param pointSize The size of the point.
#' @param removeOutlier If remove the outliers.
#' @return A ggplot object.
#' 
#' @importFrom ggplot2 scale_colour_gradientn
#' @importFrom grDevices colorRampPalette topo.colors heat.colors terrain.colors cm.colors
#' 
#' @export
#' @examples
#' x <- c(rnorm(100, mean = 1), rnorm(100, mean = 3), rnorm(100, mean = 9))
#' y <- c(rnorm(100, mean = 2), rnorm(100, mean = 8), rnorm(100, mean = 5))
#' c <- rnorm(300, 10, 5)
#' data <- data.frame(dim1 = x, dim2 = y, marker = c)
#' cytof_colorPlot(data = data, xlab = "dim1", ylab = "dim2", zlab = "marker")
cytof_colorPlot <- function(data, xlab, ylab, zlab, 
                            colorPalette = c("bluered", "topo", "heat", "terrain", "cm"), 
                            pointSize=1, 
                            removeOutlier = TRUE){
    
    remove_outliers <- function(x, na.rm = TRUE, ...) {
        qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
        H <- 1.5 * IQR(x, na.rm = na.rm)
        y <- x
        y[x < (qnt[1] - H)] <- qnt[1] - H
        y[x > (qnt[2] + H)] <- qnt[2] + H
        y
    }
    
    data <- as.data.frame(data)
    title <- paste(zlab, "Expression Level Plot")
    data <- data[,c(xlab, ylab, zlab)]
    
    if(removeOutlier)
        data[,zlab] <- remove_outliers(data[,zlab])

    colorPalette <- match.arg(colorPalette)
    switch(colorPalette,
           bluered = {
               myPalette <- colorRampPalette(c("blue", "white", "red"))
           },
           topo = {
               myPalette <- colorRampPalette(topo.colors(50))
           },
           heat = {
               myPalette <- colorRampPalette(heat.colors(50))
           },
           terrain = {
               myPalette <- colorRampPalette(terrain.colors(50))
           },
           cm = {
               myPalette <- colorRampPalette(cm.colors(50))
           }
    )
    zlength <- nrow(data)
    exp <- "Expression"
    colnames(data) <- c(xlab, ylab, exp)
    gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = exp)) + 
        scale_colour_gradientn(name = zlab, colours = myPalette(zlength)) +
        geom_point(size = pointSize) + theme_bw() + coord_fixed() +
        theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
    
    return(gp)
}


#' Statistics of the cluster relusts
#' 
#' Calculate the mean or median expression level of each marker for each cluster, or percentage
#' of cell numbers of each cluster for each sample.
#' 
#' @param data Input data frame.
#' @param markers The names of markers used for calcualtion.
#' @param cluster The column name contatining cluster labels.
#' @param sample The samples used for calculation.
#' @param statMethod Statistics contatining mean, median or percentage.
#' 
#' @return A matrix of the satatistics results
#'
#' @importFrom stats aggregate
#' @importFrom reshape2 dcast
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
#' row.names(exprs_cluster) <- rnames
#' cytof_clusterStat(data = exprs_cluster, cluster = "cluster", statMethod = "mean")
cytof_clusterStat <- function(data, markers, cluster = "cluster", sample, 
                              statMethod = c("mean", "median", "percentage", "NULL")){
    
    data <- as.data.frame(data)
    statMethod <- match.arg(statMethod)
    
    if(missing(sample)){
        sample <- "sample"
        data$sample <- sub("_[0-9]*$", "", row.names(data))
    }
    
    if(missing(markers)){
        markers <- setdiff(colnames(data), c(cluster, sample))
    }
    
    exprs_cluster <- data[ ,markers, drop=FALSE]
    exprs_cluster$cluster <- data[[cluster]]
    switch(statMethod,
           mean = {
               statData <- aggregate(. ~ cluster, data = exprs_cluster, mean)
           },
           median = {
               statData <- aggregate(. ~ cluster, data = exprs_cluster, median)
           }, 
           percentage = {
               cluster_sample <- data.frame(cluster = data[[cluster]], sample = data[[sample]])
               cluster_sample$value <- 1
               clust_sample_count <- dcast(cluster_sample, cluster ~ sample, fun.aggregate = length)
               statData <- apply(clust_sample_count[,-1], 2, function(x){round(x/sum(x)*100, 2)})
               statData <- data.frame(statData, cluster = clust_sample_count$cluster, check.names = FALSE)
           })
    
    rownames(statData) <- paste0("cluster_", statData$cluster)
    statData$cluster <- NULL  ## remove cluster column
    
    return(as.matrix(statData))
}
    
    
#' Progression plot
#' 
#' Plot the expression trend along the estimated cell progressing order
#' 
#' @param data The data frame for progression plot.
#' @param markers The column names of the selected markers for visualization.
#' @param clusters Selecte clusters for plotting, defauls select all.
#' @param orderCol The column name of the estimated cell progression order.
#' @param clusterCol The column name of the cluster results.
#' @param reverseOrder If reverse the value of orderCol.
#' @param addClusterLabel If add the cluster label on the plot.
#' @param clusterLabelSize The size of the cluster label.
#' @param segmentSize The size of the cluster label arrow.
#' @param min_expr the threshold of the minimal expression value for markers.
#' @param trend_formula a symbolic description of the model to be fit.
#' 
#' @return a ggplot2 object
#' @importFrom VGAM vgam sm.ns
#' @importFrom ggplot2 arrow unit
#' @importFrom reshape2 melt
#' @importFrom plyr ddply .
#' @importFrom ggrepel geom_text_repel
#' 
#' @export
#' 
#' @examples
#' m1 <- c(rnorm(300, 10, 2), rnorm(400, 4, 2), rnorm(300, 7))
#' m2 <- c(rnorm(300, 4), rnorm(400, 16), rnorm(300, 10, 3))
#' m3 <- c(rnorm(300, 16), rnorm(400, 40, 3), rnorm(300, 10))
#' m4 <- c(rnorm(300, 7, 3), rnorm(400, 30, 2), rnorm(300, 10))
#' m5 <- c(rnorm(300, 27), rnorm(400, 40, 1),rnorm(300, 10))
#' c <- c(rep(1,300), rep(2,400), rep(3,300))
#' rnames <- paste(paste('sample_', c('A','B','C','D'), sep = ''), 
#' rep(1:250,each = 4), sep='_') 
#' exprs_cluster <- data.frame(cluster = c, m1 = m1, m2 = m2, m3 = m3, m4 = m4, isomap_1 = m5)
#' row.names(exprs_cluster) <- sample(rnames, 1000)
#' cytof_progressionPlot(exprs_cluster, markers = c("m1","m2","m3","m4"))
cytof_progressionPlot <- function(data, markers, clusters, 
                                  orderCol="isomap_1", 
                                  clusterCol = "cluster", 
                                  reverseOrder = FALSE,
                                  addClusterLabel = TRUE,
                                  clusterLabelSize = 5,
                                  segmentSize = 0.5,
                                  min_expr = NULL, 
                                  trend_formula="expression ~ sm.ns(Pseudotime, df=3)"){
    
    if(!is.data.frame(data)) data <- data.frame(data, check.names = FALSE)
    if(!all(markers %in% colnames(data))) stop("Unmatching markers found!")
    if(!(length(orderCol)==1 && orderCol %in% colnames(data)))
        stop("Can not find orderCol in data!")
    if(!(length(clusterCol)==1 && clusterCol %in% colnames(data)))
        stop("Can not find clusterCol in data!")
    if(!missing(clusters)){
        if(!all(clusters %in% data[[clusterCol]]))
            stop("Wrong clusters selected!")
        data <- data[data[[clusterCol]] %in% clusters, , drop=FALSE]
    }
    
    if(reverseOrder){
        newOrderCol <- paste0(orderCol, "(reverse)")
        data[[newOrderCol]] <- -data[[orderCol]]
        orderCol <- newOrderCol
    }
    orderValue <- data[[orderCol]]
    data <- data[order(orderValue), c(markers, clusterCol)]
    data$Pseudotime <- sort(orderValue)
    
    mdata <- melt(data, id.vars = c("Pseudotime", clusterCol), 
                  variable.name = "markers", value.name= "expression")
    colnames(mdata) <- c("Pseudotime", clusterCol, "markers", "expression")
    mdata$markers <- factor(mdata$markers)
    mdata[[clusterCol]] <- factor(mdata[[clusterCol]])
    min_expr <- min(mdata$expression)
    
    ## tobit regression
    vgamPredict <- ddply(mdata, .(markers), function(x) { 
        fit_res <- tryCatch({
            vg <- suppressWarnings(vgam(formula = as.formula(trend_formula), 
                                        family = VGAM::tobit(Lower = min_expr, lmu = "identitylink"), 
                                        data = x, maxit=30, checkwz=FALSE))
            res <- VGAM::predict(vg, type="response")
            res[res < min_expr] <- min_expr
            res
        }
        ,error = function(e) {
            print("Error!")
            print(e)
            res <- rep(NA, nrow(x))
            res
        }
        )
        expectation = fit_res
        data.frame(Pseudotime=x[["Pseudotime"]], expectation=expectation)
    })
    
    color_by <- clusterCol
    plot_cols <- round(sqrt(length(markers)))
    cell_size <- 1
    x_lab <- orderCol
    y_lab <- "Expression"
    legend_title <- "Cluster"
    
    ## copied from monocle package
    monocle_theme_opts <- function(){
        theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
            #theme(panel.border = element_blank(), axis.line = element_line()) +
            theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
            theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
            theme(panel.background = element_rect(fill='white')) +
            theme(legend.position = "right") +
            theme(axis.title = element_text(size = 15)) +
            theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))}
    
    q <- ggplot(aes_string(x="Pseudotime", y="expression"), data=mdata) 
    q <- q + geom_point(aes_string(color=color_by), size=I(cell_size))
    q <- q + geom_line(aes_string(x="Pseudotime", y="expectation"), data=vgamPredict, size = 0.9)
    q <- q + facet_wrap(~markers, ncol=plot_cols, scales="free_y")
    q <- q + ylab(y_lab) + xlab(x_lab) + theme_bw()
    q <- q + guides(colour = guide_legend(title = legend_title, override.aes = list(size = cell_size*3)))
    q <- q + monocle_theme_opts() 
    
    if(addClusterLabel){
        # edata <- data[ ,c("Pseudotime", clusterCol)]
        # colnames(edata) <- c('x', "z")
        # center <- aggregate(x ~ z, data = edata, median)
        # center$y <- -0.5 ## add to the botom
        # q <- q + geom_text_repel(data=center, aes(x=x, y=y, label=z), parse=TRUE)
        mdata$cluster <- mdata[[clusterCol]]
        center <- aggregate(cbind(Pseudotime, expression) ~ cluster + markers, data = mdata, median)
        q <- q + geom_text_repel(data=center, aes(x=Pseudotime, y=expression, label=cluster),
                                 size = clusterLabelSize, fontface = 'bold',
                                 box.padding = unit(0.5, 'lines'),
                                 point.padding = unit(1.6, 'lines'),
                                 segment.color = '#555555',
                                 segment.size = segmentSize,
                                 arrow = arrow(length = unit(0.02, 'npc')))
    }

    q
}

#' Add data to the original FCS files
#' 
#' Store the new dimension transformed data and cluster data into the exprs 
#' matrix in new fcs files under \code{analyzedFCSdir}
#' 
#' @param data The new data matrix to be added in.
#' @param rawFCSdir The directory containing the original fcs files.
#' @param analyzedFCSdir The directory to store the new fcs files.
#' @param transformed_cols the column name of the dimension transformend data in \code{data}.
#' @param cluster_cols the column name of the cluster data in \code{data}.
#' @param inLgclTrans Boolean value decides if apply the inverse lgcl transformation to the data before saving
#' 
#' @export
#' 
#' @return new fcs files stored under \code{analyzedFCSdir}
#' @importMethodsFrom flowCore keyword
#' @importFrom Biobase  exprs exprs<- description description<- pData pData<- 
cytof_addToFCS <- function(data, rawFCSdir, analyzedFCSdir, transformed_cols = c("tsne_1", 
    "tsne_2"), cluster_cols = c("cluster"), inLgclTrans = TRUE) {
    
    lgcl <- logicleTransform(w = 0.1, t = 4000, m = 4.5, a = 0)
    ilgcl <- inverseLogicleTransform(trans = lgcl)
    
    if (!file.exists(analyzedFCSdir)) {
        dir.create(analyzedFCSdir)
    }
    
    ## transform transformed_cols
    if (!is.null(transformed_cols)) {
        transformed <- data[, transformed_cols, drop=FALSE]
        # if(inLgclTrans){
        #     N_transformed <- apply(transformed, 2, function(x)
        #         ((x-min(x))/(max(x)-min(x)))*4.4) 
        #     R_N_transformed <- apply(N_transformed,2,ilgcl)
        # }else{
        #     R_N_transformed <- apply(transformed, 2, function(x) (x - min(x)) + 0.1)
        # }
        
        R_N_transformed <- apply(transformed, 2, function(x) (x - min(x)) + 0.1)  ## no inLgclTrans
        row.names(R_N_transformed) <- row.names(data)
    }
    
    ## transform cluster_cols
    if (!is.null(cluster_cols)) {
        if(inLgclTrans){
            for(i in 1:length(cluster_cols)){
                cCol <- data[, cluster_cols[i]]
                clust_cor_1 <- as.numeric(cCol)%%10
                clust_cor_2 <- floor(as.numeric(cCol)/10)
                clust_cor_1 <- clust_cor_1 + runif(length(clust_cor_1), 0, 0.2)
                clust_cor_2 <- clust_cor_2 + runif(length(clust_cor_2), 0, 0.2)
                cluster_cor12 <- cbind(clust_cor_1, clust_cor_2)
                N_clust_cor <- apply(cluster_cor12, 2, function(x) ((x - min(x))/(max(x) - min(x))) * 4.4)
                clust_cor <- apply(N_clust_cor, 2, ilgcl)
                colnames(clust_cor) <- paste(cluster_cols[i], c("cor_1", "cor_2"), sep="_") 
                if(i == 1){
                    R_N_clust_cor <- clust_cor
                }else{
                    R_N_clust_cor <- cbind(R_N_clust_cor, clust_cor)
                }
            }
            
        }else{
            R_N_clust_cor <- data[, cluster_cols, drop = FALSE]
        }
        
        row.names(R_N_clust_cor) <- row.names(data)
    }
    
    ## write data to FCS
    if ((!is.null(transformed_cols)) && (!is.null(cluster_cols))) {
        to_add <- cbind(R_N_transformed, R_N_clust_cor)
    } else if (!is.null(transformed_cols)) {
        to_add <- R_N_transformed
    } else if (!is.null(cluster_cols)) {
        to_add <- R_N_clust_cor
    }
    addColNames <- colnames(to_add)
    sample <- unique(sub("_[0-9]*$", "", row.names(to_add)))
    
    for (i in 1:length(sample)) {
    	fn <- paste0(rawFCSdir, .Platform$file.sep, sample[i], ".fcs")
    	if(!file.exists(fn)){
    	    ## stop the writing if cannot find the file
    	    message(paste("rawFCSdir is not correct, can not find raw FCS file:", fn))
    	    return(NULL)
    	}
        cat("Save to file:", fn, "\n")
        fcs <- read.FCS(fn, transformation = FALSE)
        pattern <- paste(sample[i], "_", sep = "")
        to_add_i <- as.data.frame(to_add[grep(pattern, row.names(to_add), fixed = TRUE), ])
        m <- regexpr("_[0-9]*$", row.names(to_add_i))
        cellNo_i <- as.integer(substring(regmatches(row.names(to_add_i), m), 2))
        sub_exprs <- fcs@exprs[cellNo_i, ]
        params <- parameters(fcs)
        pd <- pData(params)
        keyval <- keyword(fcs)
        
        for (j in 1:length(addColNames)) {
            ## update the parameters
            if (addColNames[j] %in% colnames(fcs@exprs)) {
                addColNames[j] <- paste(addColNames[j], "_new", sep = "")
            }
            
            addColName <- addColNames[j]
            channel_number <- nrow(pd) + 1
            channel_id <- paste("$P", channel_number, sep = "")
            channel_name <- addColName
            minRange <- ceiling(min(to_add_i[[j]]))
            maxRange <- ceiling(max(to_add_i[[j]]))
            channel_range <- maxRange - minRange
            
            plist <- matrix(c(channel_name, "<NA>", channel_range, 
                minRange, maxRange))
            rownames(plist) <- c("name", "desc", "range", "minRange", 
                "maxRange")
            colnames(plist) <- c(channel_id)
            pd <- rbind(pd, t(plist))
            
            ## update the expression value
            out_col_names <- colnames(sub_exprs)
            sub_exprs <- cbind(sub_exprs, to_add_i[[j]])
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
        out_frame <- flowFrame(exprs = sub_exprs, parameters = params, description = keyval)
        
        suppressWarnings(write.FCS(out_frame, paste0(analyzedFCSdir, "/cytofkit_", sample[i], ".fcs")))
    }
} 
