## loading pacakge
require(cytofkit)
require(ggplot2)


## Main function for scatter plot
scatterPlot <- function(obj, plotMethod, plotFunction, pointSize=1, 
                      addLabel=TRUE, labelSize=1, sampleLabel = TRUE,
                      FlowSOM_k = 40, selectSamples, facetPlot = FALSE, 
                      colorPalette = "bluered", labelRepel = FALSE, removeOutlier = TRUE){
    
    data <- cbind(obj$expressionData, 
                  obj$dimReducedRes[[plotMethod]], 
                  do.call(cbind, obj$clusterRes))
    data <- as.data.frame(data)
    xlab <- colnames(obj$dimReducedRes[[plotMethod]])[1]
    ylab <- colnames(obj$dimReducedRes[[plotMethod]])[2]
    row.names(data) <- row.names(obj$expressionData)
    
    clusterMethods <- names(obj$clusterRes)
    samples <- sub("_[0-9]*$", "", row.names(obj$expressionData))
    data <- data[samples %in% selectSamples, ]
    nsamples <- samples[samples %in% selectSamples]
    data$sample <- nsamples
    sample_num <- length(unique(nsamples))

    if(plotFunction == "DensityPlot"){
        colPalette <- colorRampPalette(c("blue", "turquoise", "green", 
                                         "yellow", "orange", "red"))
        densCol <- densCols(data[, c(xlab, ylab)], colramp = colPalette)
        data$densCol <- densCol
        gp <- ggplot(data, aes_string(x=xlab, y=ylab)) + 
            geom_point(colour=densCol, size = pointSize) + ggtitle("Density Plot") +
            theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
    }else if(plotFunction == "DotPlot"){
        gp <- ggplot(data, aes_string(x=xlab, y=ylab)) + 
            geom_point(size = pointSize) + ggtitle("Dot Plot") +
            xlab(xlab) + ylab(ylab) + theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
    }else if(plotFunction == "ColorBySample"){
        size_legend_row <- ceiling(sample_num/4)
        sample <- "sample"
        gp <- ggplot(data, aes_string(x=xlab, y=ylab, colour = sample)) +
            geom_point(size = pointSize) + ggtitle("Color By Sample") +
            xlab(xlab) + ylab(ylab) + theme_bw() + theme(legend.position = "bottom") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold")) +
            guides(colour = guide_legend(nrow = size_legend_row, override.aes = list(size = 4)))
    }else if(plotFunction %in% clusterMethods){
        gp <- cytof_clusterPlot(data = data, 
                                xlab = xlab, 
                                ylab = ylab, 
                                cluster = plotFunction, 
                                sample = "sample",
                                title = plotFunction, 
                                type = ifelse(facetPlot, 2, 1),
                                point_size = pointSize, 
                                addLabel = addLabel, 
                                labelSize = labelSize, 
                                sampleLabel = sampleLabel,
                                labelRepel = labelRepel,
                                fixCoord = FALSE)
    }else{
        gp <- cytof_colorPlot(data = data, 
                              xlab = xlab, 
                              ylab = ylab, 
                              zlab = plotFunction, 
                              colorPalette = colorPalette,
                              pointSize = pointSize, 
                              removeOutlier = TRUE)
    }
    
    return(gp)
}


heatMap <- function(data, clusterMethod = "DensVM", type = "mean", selectSamples,
                    cex_row_label = 1, cex_col_label = 1, scaleMethod = "none") {
    exprs <- data$expressionData
    samples <- sub("_[0-9]*$", "", row.names(exprs))
    mySamples <- samples %in% selectSamples
    exprs <- exprs[mySamples, , drop = FALSE]
    dataj <- data$clusterRes[[clusterMethod]][mySamples]
    exprs_cluster <- data.frame(exprs, cluster = dataj, check.names = FALSE )
    
    cluster_stat <- cytof_clusterStat(data = exprs_cluster,
                             cluster = "cluster", 
                             statMethod = type)
    
    cytof_heatmap(data = as.matrix(cluster_stat), 
                  baseName = paste(clusterMethod, type, "Heat Map"), 
                  scaleMethod = scaleMethod, 
                  cex_row_label = cex_row_label, 
                  cex_col_label = cex_col_label,
                  margins = c(8, 8), 
                  keysize = 1, 
                  key.par=list(mgp=c(2, 1, 0), mar=c(4, 3, 4, 0))) 
}

## function for opening the results directory
opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
        shell.exec(dir)
    } else {
        system(paste(Sys.getenv("R_BROWSER"), dir))
    }
}


