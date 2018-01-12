## loading package
require(cytofkit)
require(ggplot2)
require(reshape2)
require(plyr)
require(VGAM)
require(colourpicker)
require(gplots)


## Main function for scatter plot
scatterPlot <- function(obj, plotMethod, plotFunction, pointSize=1, alpha = 1,
                      addLabel=TRUE, labelSize=1, sampleLabel = TRUE,
                      FlowSOM_k = 40, selectCluster=NULL, selectSamples, 
                      facetPlot = FALSE, colorPalette = "bluered", labelRepel = FALSE, 
                      removeOutlier = TRUE, clusterColor, globalScale = TRUE, centerScale = FALSE){
    
    data <- data.frame(obj$expressionData, 
                       obj$dimReducedRes[[plotMethod]], 
                       do.call(cbind, obj$clusterRes), 
                       check.names = FALSE,
                       stringsAsFactors = FALSE)
    
    Markers <- obj$allMarkers
    
    xlab <- colnames(obj$dimReducedRes[[plotMethod]])[1]
    ylab <- colnames(obj$dimReducedRes[[plotMethod]])[2]
    row.names(data) <- row.names(obj$expressionData)
    
    clusterMethods <- names(obj$clusterRes)
    samples <- sub("_[0-9]*$", "", row.names(obj$expressionData))
    data <- data[samples %in% selectSamples, ,drop=FALSE]
    nsamples <- samples[samples %in% selectSamples]
    data$sample <- nsamples
    sample_num <- length(unique(nsamples))
    
    if(plotFunction == "Density"){
        colPalette <- colorRampPalette(c("blue", "turquoise", "green", 
                                         "yellow", "orange", "red"))
        densCol <- densCols(data[, c(xlab, ylab)], colramp = colPalette)
        data$densCol <- densCol
        gp <- ggplot(data, aes_string(x=xlab, y=ylab)) + 
            geom_point(colour=densCol, size = pointSize) + ggtitle("Density Plot") +
            theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
    }else if(plotFunction == "None"){
        gp <- ggplot(data, aes_string(x=xlab, y=ylab)) + 
            geom_point(size = pointSize) + ggtitle("Dot Plot") +
            xlab(xlab) + ylab(ylab) + theme_bw() + 
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold"))
    }else if(plotFunction == "Sample"){
        size_legend_row <- ceiling(sample_num/4)
        sample <- "sample"
        gp <- ggplot(data, aes_string(x=xlab, y=ylab, colour = sample)) +
            geom_point(size = pointSize) + ggtitle("Color By Sample") +
            xlab(xlab) + ylab(ylab) + theme_bw() + theme(legend.position = "bottom") +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
            theme(axis.text=element_text(size=14), axis.title=element_text(size=18,face="bold")) +
            guides(colour = guide_legend(nrow = size_legend_row, override.aes = list(size = 4)))
    }else if(plotFunction == "All Markers"){
        gp <- cytof_wrap_colorPlot(data = data, 
                              xlab = xlab, 
                              ylab = ylab, 
                              markers = colnames(obj$expressionData), 
                              colorPalette = colorPalette,
                              limits = NULL,
                              pointSize = pointSize, 
                              removeOutlier = TRUE)
        
    }else if(plotFunction == "All Markers(scaled)"){
        gp <- cytof_wrap_colorPlot(data = data, 
                                   xlab = xlab, 
                                   ylab = ylab, 
                                   markers = colnames(obj$expressionData), 
                                   scaleMarker = TRUE,
                                   colorPalette = colorPalette,
                                   limits = NULL,
                                   pointSize = pointSize, 
                                   removeOutlier = TRUE)
        
    }else if(plotFunction %in% clusterMethods){
        
        if(!is.null(selectCluster)){
            clusterIDs <- as.character(data[,plotFunction])
            selectCluster <- as.character(selectCluster)
            data <- data[clusterIDs %in% selectCluster, ,drop=FALSE]
        }
        clusterVec <- obj$clusterRes[[plotFunction]]
        ## make sure they are not factors before transforming to factors
        selectColors <- match(levels(as.factor(data[,plotFunction])), levels(as.factor(clusterVec)))
        clusterColor <- clusterColor[selectColors]
        
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
                                fixCoord = FALSE,
                                clusterColor = clusterColor)
    }else{
        limits <- NULL
        if(globalScale){
          exprData <- obj$expressionData
          markers <- colnames(exprData)
          glimits <- quantile(exprData, probs=c(.02, .98), na.rm = TRUE)
          local.bounds <- as.data.frame(lapply(markers, function(x) quantile(exprData[,x], probs=c(.02, .98), na.rm = TRUE)), col.names = markers)
          gmax <- ifelse(max(local.bounds[2,]) < glimits[2], glimits[2], max(local.bounds[2,]))
          gmin <- ifelse(min(local.bounds[1,]) > glimits[1],min(local.bounds[1,]), glimits[1])
          limits <- c(gmin, gmax)
        }
        if(length(plotFunction > 1)){
          gp <- cytof_wrap_colorPlot(data = data, 
                                     xlab = xlab, 
                                     ylab = ylab, 
                                     markers = plotFunction, 
                                     colorPalette = colorPalette,
                                     limits = limits,
                                     scaleMarker = centerScale,
                                     pointSize = pointSize,
                                     alpha = alpha,
                                     removeOutlier = TRUE)
        }else{
          gp <- cytof_colorPlot(data = data, 
                                xlab = xlab, 
                                ylab = ylab, 
                                zlab = plotFunction, 
                                colorPalette = colorPalette,
                                limits = limits,
                                pointSize = pointSize,
                                alpha = alpha,
                                removeOutlier = TRUE)
        }
    }
    
    return(gp)
}

## Facet wrap plot of marker expression
cytof_wrap_colorPlot <- function(data, xlab, ylab, markers, scaleMarker = FALSE,
                            colorPalette = c("bluered", "spectral1", "spectral2", "heat"),
                            limits = NA,
                            pointSize=1,
                            alpha = 1,
                            removeOutlier = TRUE){
    
    remove_outliers <- function(x, na.rm = TRUE, ...) {
        qnt <- quantile(x, probs=c(.02, .98), na.rm = na.rm, ...)
        x[x <= qnt[1]] <- qnt[1]
        x[x >= qnt[2]] <- qnt[2]
        x
    }
    
    data <- as.data.frame(data)
    title <- "Marker Expression Level Plot"
    data <- data[,c(xlab, ylab, markers)]
    
    if(removeOutlier){
        for(m in markers){
            data[[m]] <- remove_outliers(data[ ,m])
        }
    }
    
    if(scaleMarker){
        data[ ,markers] <- scale(data[ ,markers], center = TRUE, scale = TRUE)
        ev <- "ScaledExpression"
        data <- melt(data, id.vars = c(xlab, ylab), 
                     measure.vars = markers,
                     variable.name = "markers", 
                     value.name = ev)
    }else{
        ev <- "Expression"
        data <- melt(data, id.vars = c(xlab, ylab), 
                     measure.vars = markers,
                     variable.name = "markers", 
                     value.name = ev)
    }
    

    colorPalette <- match.arg(colorPalette)
    switch(colorPalette,
           bluered = {
               myPalette <- colorRampPalette(c("blue", "white", "red"))
           },
           spectral1 = {
               myPalette <- colorRampPalette(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4",
                                               "#E6F598", "#FFFFBF", "#FEE08B", "#FDAE61",
                                               "#F46D43", "#D53E4F", "#9E0142"))
           },
           spectral2 = {
               myPalette <- colorRampPalette(rev(c("#7F0000","red","#FF7F00","yellow","white", 
                                                   "cyan", "#007FFF", "blue","#00007F")))
           },
           heat = {
               myPalette <- colorRampPalette(heat.colors(50))
           }
    )
    zlength <- nrow(data)
    grid_row_num <- round(sqrt(length(markers)))
    gp <- ggplot(data, aes_string(x = xlab, y = ylab, colour = ev)) + 
        facet_wrap(~markers, nrow = grid_row_num, scales = "fixed") +
        scale_colour_gradientn(limits = limits, name = ev, colours = myPalette(zlength * 2)) +
        geom_point(size = pointSize, alpha = alpha) + theme_bw() + coord_fixed() +
        theme(legend.position = "right") + xlab(xlab) + ylab(ylab) + ggtitle(title) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(axis.text=element_text(size=8), axis.title=element_text(size=12,face="bold"))
    
    return(gp)
}

## Heat Map
heatMap <- function(data, clusterMethod = "DensVM", type = "mean", 
                    dendrogram = "both", colPalette = "bluered", selectSamples, selectMarkers = NULL,
                    cex_row_label = 1, cex_col_label = 1, scaleMethod = "none") {
    exprs <- data$expressionData
    samples <- sub("_[0-9]*$", "", row.names(exprs))
    if(!(is.null(selectMarkers))) {
      marker_id <- selectMarkers
    }else{
      marker_id <- colnames(exprs)
    }
    markers <- colnames(exprs)
    mySamples <- samples %in% selectSamples
    myMarkers <- markers %in% marker_id
    exprs <- exprs[mySamples, , drop = FALSE]
    exprs <- exprs[, myMarkers, drop = FALSE]
    dataj <- data$clusterRes[[clusterMethod]][mySamples]
    exprs_cluster <- data.frame(exprs, cluster = dataj, check.names = FALSE )
    
    cluster_stat <- cytof_clusterStat(data = exprs_cluster,
                             cluster = "cluster", 
                             statMethod = type)
    
    cytof_heatmap(data = as.matrix(cluster_stat), 
                  baseName = paste(clusterMethod, type), 
                  scaleMethod = scaleMethod, 
                  dendrogram = dendrogram,
                  colPalette = colPalette,
                  cex_row_label = cex_row_label, 
                  cex_col_label = cex_col_label,
                  margins = c(8, 8), 
                  keysize = 1, 
                  key.par=list(mgp=c(1.5, 0.5, 0), mar=c(3, 2.5, 3.5, 1))) 
}

## density plot

#' @param densData Data frame.
#' @param stackRotation Rotation degree of density plot to the right side, range (0-90).
#' @param stackSeperation Control factor for stack seperation interval, numeric value from 0-1, or auto.
#'
#' @importFrom plyr ldply
#' @importFrom reshape2 melt
#' @import ggplot2
stackDenistyPlot <- function(data, densityCols, stackFactor,
                             kernel = c("gaussian", "epanechnikov", "rectangular",
                                        "triangular", "biweight",
                                        "cosine", "optcosine"),
                             bw = "nrd0", adjust = 1,
                             reomoveOutliers = FALSE, 
                             stackRotation = 0, 
                             stackSeperation = "auto",
                             x_text_size = 2, 
                             strip_text_size = 7,
                             legend_text_size = 0.5, 
                             legendRow = 1,
                             legend_title = "stackName",
                             stackFactorColours = NULL){
    
    if(!is.numeric(stackRotation)){
        stop("stackRotation must be a numeric number")
    }else if(stackRotation < 0 || stackRotation > 90){
        stop("stackRotation must be a numeric number in range 0-90")
    }
    
    if(missing(densityCols)){
        densityCols <- colnames(data)
    }else if(any(!(densityCols %in% colnames(data)))){
        stop("Unmatch densityCols found:", paste(densityCols[!(densityCols %in% colnames(data))], collapse = " "))
    }
    
    if(missing(stackFactor)){
        warning("no stackFactor was provided!")
        stackFactor <- rep("stack", length = nrow(data))
    }else if(length(stackFactor) != nrow(data)){
        stop("Length of stackFactor unequal row number of input data")
    }
    kernel <- match.arg(kernel)
    
    stackCount <- length(unique(stackFactor))
    densityCount <- length(densityCols)
    
    if(missing(stackFactorColours) || is.null(stackFactorColours)){
        stackFactorColours <- rainbow(stackCount)
    }else if(length(stackFactorColours) == 0 || length(stackFactorColours) != stackCount){
        stackFactorColours <- rainbow(stackCount)
    }
    
    data <- data.frame(data[ ,densityCols, drop=FALSE], stackFactor = stackFactor, check.names = FALSE)
    
    densData <- .densityCal(data, kernel = kernel, bw = bw, adjust = adjust, reomoveOutliers = reomoveOutliers)
    ## dataframe densData contains {stackName, x , y , densityName}
    xStat <- aggregate(x ~ stackName + densityName, densData, max)
    yStat <- aggregate(y ~ stackName + densityName, densData, max)
    
    if(stackSeperation == "auto"){
        stackIntervals <- aggregate(y ~ densityName, yStat, function(x){0.8*median(x) * (1-(stackRotation/90)^0.2)^2})
    }else if(stackSeperation < 0 || stackSeperation > 1){
        stop("stackSeperation must be value in range 0-1")
    }else{
        stackIntervals <- aggregate(y ~ densityName, yStat, function(x){median(x)*stackSeperation})
    }
    
    stackShifts <- aggregate(x ~ densityName, xStat, function(x){max(x) * (stackRotation/90)})
    
    densData$stack_x <- densData$x + (as.numeric(densData$stackName)-1) * stackShifts$x[match(densData$densityName, stackShifts$densityName)]
    densData$stack_y <- densData$y + (as.numeric(densData$stackName)-1) * stackIntervals$y[match(densData$densityName, stackIntervals$densityName)]
    
    ## segment lines, x tick, x label
    alignSegments <- ldply(split(densData$x, densData$densityName),
                           function(x){seq(min(x), max(x), length.out=5)},
                           .id = "densityName")
    alignSegments <- melt(alignSegments, id.vars="densityName", variable.name="x_tick", value.name = "x")
    alignSegments$y <- min(densData$y)
    alignSegments$xend <- alignSegments$x + (length(unique(densData$stackName))-1) * stackShifts$x[match(alignSegments$densityName, stackShifts$densityName)]
    alignSegments$yend <- min(densData$y) + (length(unique(densData$stackName))-1) * stackIntervals$y[match(alignSegments$densityName, stackIntervals$densityName)]
    
    densityHeights <- aggregate(y ~ densityName, yStat, max)
    alignSegments$tickXend <- alignSegments$x
    alignSegments$tickYend <- alignSegments$y - densityHeights$y[match(alignSegments$densityName, densityHeights$densityName)] * 0.01
    alignSegments$tickText <- format(alignSegments$x,scientific=TRUE, digits=3)
    alignSegments$textY <- alignSegments$y - densityHeights$y[match(alignSegments$densityName, densityHeights$densityName)] * 0.03
    
    cat(" Plotting ...\n")
    stackDensityPlot_theme <- theme(legend.position = "top",
                                    legend.title = element_text(size = rel(1)),
                                    legend.text = element_text(size = rel(legend_text_size)),
                                    strip.text = element_text(size=strip_text_size, lineheight=1, hjust = 0.5, vjust = 0.5),
                                    axis.text.x = element_blank(),
                                    axis.ticks.x = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.ticks.y = element_blank(),
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    panel.border = element_blank(),
                                    strip.background=element_rect(fill = "grey90", colour = NA))
    
    gp <- ggplot(densData, aes(x=stack_x, y=stack_y)) +
        geom_segment(data = alignSegments,
                     aes(x = x, y = y, xend = xend, yend = yend),
                     color = "grey80", size=0.3) +
        geom_segment(data = alignSegments,
                     aes(x = x, y = y, xend = tickXend, yend = tickYend),
                     color = "grey20", size=0.3) +
        geom_text(data = alignSegments, aes(x = x, y = textY, label = tickText),
                  hjust = 0.3, vjust = 1.1, size = x_text_size) +
        geom_polygon(aes(fill=stackName, color=stackName), alpha = 0.15) + 
        scale_colour_manual(values = stackFactorColours) + 
        scale_fill_manual(values = stackFactorColours) +
        facet_wrap(~densityName, scale = "free") +
        xlab("") + ylab("") +
        guides(col = guide_legend(title = legend_title, nrow = legendRow, byrow = TRUE),
               fill = guide_legend(title = legend_title, nrow = legendRow, byrow = TRUE)) +
        theme_bw() + stackDensityPlot_theme
    
    gp
}


#' Internal density calculation function serves for \code{stackDenistyPlot}
#'
#' Output data frame with columns: stackName, x , y , densityName
.densityCal <- function(data, kernel, bw, adjust, reomoveOutliers = FALSE){
    cat("  Calculating Density for each stack column...\n")
    print(table(data$stackFactor))
    dataBystackFactor <- split(subset(data, select = -stackFactor), data$stackFactor)
    densityWrap <- function(d, ...){
        resOut <- NULL
        for(i in colnames(d)){
            x <- d[,i]
            if(reomoveOutliers){
                cat("  Remove outliers...\n")
                x_IQR <- IQR(x)
                x_lowLimit <- quantile(x, 0.25) - 1.5*x_IQR
                x_highLimit <- quantile(x, 0.75) + 1.5*x_IQR
                x <- x[x >= x_lowLimit && x <= x_highLimit]
            }
            dens <- density(x, ...)
            densOut <- data.frame(x=dens$x, y=dens$y, densityName = i)
            resOut <- rbind(resOut, densOut)
        }
        return(resOut)
    }
    
    r <- ldply(dataBystackFactor, densityWrap,
               kernel = kernel, bw = bw, adjust = adjust,
               .progress = "text",
               .id = "stackName")
    return(r)
}


## Combined marker expression trend
cytof_expressionTrends <- function(data, markers, clusters, 
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
    
    q <- ggplot(data=vgamPredict, aes_string(x="Pseudotime", y="expectation", col="markers")) + geom_line(size = 1.5)
    q <- q + ylab(y_lab) + xlab(x_lab) + theme_bw()
    q <- q + guides(colour = guide_legend(title = legend_title, override.aes = list(size = cell_size*3)))
    q <- q + monocle_theme_opts() 
    
    # if(addClusterLabel){
    #     # edata <- data[ ,c("Pseudotime", clusterCol)]
    #     # colnames(edata) <- c('x', "z")
    #     # center <- aggregate(x ~ z, data = edata, median)
    #     # center$y <- -0.5 ## add to the botom
    #     # q <- q + geom_text_repel(data=center, aes(x=x, y=y, label=z), parse=TRUE)
    #     mdata$cluster <- mdata[[clusterCol]]
    #     center <- aggregate(cbind(Pseudotime, expression) ~ cluster + markers, data = mdata, median)
    #     q <- q + geom_text_repel(data=center, aes(x=Pseudotime, y=expression, label=cluster),
    #                              size = clusterLabelSize, fontface = 'bold',
    #                              box.padding = unit(0.5, 'lines'),
    #                              point.padding = unit(1.6, 'lines'),
    #                              segment.color = '#555555',
    #                              segment.size = segmentSize,
    #                              arrow = arrow(length = unit(0.02, 'npc')))
    # }
    
    q
}

## function for opening the results directory
opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
        shell.exec(dir)
    } else {
        system(paste(Sys.getenv("R_BROWSER"), dir))
    }
}


