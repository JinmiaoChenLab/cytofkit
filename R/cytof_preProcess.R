#' Merge the expression matrix from multiple FCS files with preprocessing
#' 
#' Apply preprocessing on each FCS file including compensation (for FCM data only) and transformation
#' with selected markers, then expression matrix are extracted and merged using one of the methods, 
#' \code{all}, \code{min}, \code{fixed} or \code{ceil}
#' 
#' @param fcsFiles A vector of FCS file names.
#' @param comp Either boolean value tells if do compensation (compensation matrix contained in FCS), or a compensation matrix to be applied.
#' @param markers Selected markers for analysis, either marker names/descriptions or marker IDs.
#' @param transformMethod Data Transformation method, including \code{cytofAsinh} (suggest for CyTOF data), \code{autoLgcl} (suggest for FCM data), \code{logicle} and \code{arcsinh}.
#' @param scaleTo Scale the expression to a specified range c(a, b), default is NULL.
#' @param mergeMethod Merge method for mutiple FCS expression data. cells can be combined using one of the four different methods including \code{ceil}, \code{all}, \code{min}, \code{fixed}. The default option is 
#' \code{ceil}, up to a fixed number (specified by \code{fixedNum}) of cells are sampled without replacement from each fcs file and combined for analysis.
#' \code{all}: all cells from each fcs file are combined for analysis. 
#' \code{min}: The minimum number of cells among all the selected fcs files are sampled from each fcs file and combined for analysis. 
#' \code{fixed}: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each fcs file and combined for analysis.
#' @param fixedNum The fixed number of cells to be extracted from each FCS file.
#' @param sampleSeed A sampling seed for reproducible expression matrix merging.
#' @param ... Other arguments passed to \code{cytof_exprsExtract}
#' 
#' @return A matrix containing the merged expression data, with selected markers, row names added as \code{filename_cellID}, column mamed added as \code{name<desc>}.
#' @seealso \code{\link{cytof_exprsExtract}}
#' @export
#' @examples
#' d<-system.file('extdata',package='cytofkit')
#' fcsFiles <- list.files(d,pattern='.fcs$',full=TRUE)
#' parameters <- list.files(d, pattern='.txt$', full=TRUE)
#' markers <- as.character(read.table(parameters, sep = "\t", header = TRUE)[, 1])
#' merged <- cytof_exprsMerge(fcsFiles, markers = markers)
cytof_exprsMerge <- function(fcsFiles, 
                             comp = FALSE, 
                             markers = NULL, 
                             transformMethod = "cytofAsinh", 
                             scaleTo = NULL, 
                             mergeMethod = c("ceil", "all", "fixed", "min"), 
                             fixedNum = 10000, 
                             sampleSeed = 123, ...) {
    
    exprsL <- mapply(cytof_exprsExtract, fcsFiles, 
                     MoreArgs = list(comp = comp, markers = markers, 
                                     transformMethod = transformMethod, 
                                     scaleTo = scaleTo, ...), 
                     SIMPLIFY = FALSE)

    mergeMethod <- match.arg(mergeMethod)
    if(is.numeric(sampleSeed))
        set.seed(sampleSeed)
    switch(mergeMethod,
           ceil = {
               mergeFunc <- function(x) {
                   if (nrow(x) < fixedNum) {
                       x
                   } else {
                       x[sample(nrow(x), size = fixedNum, replace = FALSE), , drop = FALSE]
                   }
               }
               merged <- do.call(rbind, lapply(exprsL, mergeFunc))
           },
           all = {
               merged <- do.call(rbind, exprsL)
           },
           fixed = {
               mergeFunc <- function(x) {
                   x[sample(nrow(x), size = fixedNum, replace = ifelse(nrow(x) < fixedNum, TRUE, FALSE)), , drop=FALSE]
               }
               merged <- do.call(rbind, lapply(exprsL, mergeFunc))
           },
           min = {
               minSize <- min(sapply(exprsL, nrow))
               mergeFunc <- function(x) {
                   x[sample(nrow(x), size = minSize, replace = FALSE), , drop=FALSE]
               }
               merged <- do.call(rbind, lapply(exprsL, mergeFunc))
           })
    
    return(merged)
}


#' Extract the expression data from a FCS file with preprocessing
#' 
#' Extract the FCS expresssion data with preprocessing of compensation (for FCM data only)
#' and transformation. Transformtion methods includes \code{cytofAsinh} (suggest for CyTOF data), 
#' \code{autoLgcl} (suggest for FCM data), \code{logicle} (customizable) and \code{arcsinh} (customizable).
#' 
#' @param fcsFile The name of the FCS file.
#' @param verbose Boolean value detecides if print the massage details of FCS loading.
#' @param comp Either boolean value tells if do compensation (compensation matrix contained in FCS), or a compensation matrix to be applied.
#' @param markers Selected markers for analysis, either marker names/descriptions or marker IDs.
#' @param transformMethod Data Transformation method, including \code{cytofAsinh} (suggest for CyTOF data), \code{autoLgcl} (suggest for FCM data), \code{logicle} and \code{arcsinh}.
#' @param scaleTo Scale the expression to a specified range c(a, b), default is NULL.
#' @param q quantile of negative values removed for auto w estimation, default is 0.05, parameter for autoLgcl transformation.
#' @param l_w Linearization width in asymptotic decades, parameter for logicle transformation.
#' @param l_t Top of the scale data value, parameter for logicle transformation.
#' @param l_m Full width of the transformed display in asymptotic decades, parameter for logicle transformation.
#' @param l_a Additional negative range to be included in the display in asymptotic decades, parameter for logicle transformation.
#' @param a_a Positive double that corresponds to the base of the arcsinh transformation, \code{arcsinh} = asinh(a + b * x) + c).
#' @param a_b Positive double that corresponds to a scale factor of the arcsinh transformation, \code{arcsinh} = asinh(a + b * x) + c).
#' @param a_c Positive double that corresponds to another scale factor of the arcsinh transformation, \code{arcsinh} = asinh(a + b * x) + c).
#' 
#' @return A transformend expression data matrix with selected markers, row names added as \code{filename_cellID}, column mamed added as \code{name<desc>}.
#' @importFrom flowCore read.FCS compensate estimateLogicle logicleTransform parameters transformList arcsinhTransform biexponentialTransform
#' @importMethodsFrom flowCore transform
#' @importClassesFrom flowCore transformList
#' @export
#' @examples
#' d <- system.file('extdata',package='cytofkit')
#' fcsFile <- list.files(d,pattern='.fcs$',full=TRUE)
#' parameters <- list.files(d, pattern='.txt$', full=TRUE)
#' markers <- as.character(read.table(parameters, sep = "\t", header = TRUE)[, 1])
#' transformed <- cytof_exprsExtract(fcsFile, markers = markers)
cytof_exprsExtract <- function(fcsFile, 
                               verbose = FALSE, 
                               comp = FALSE, 
                               markers = NULL, 
                               transformMethod = c("cytofAsinh", "autoLgcl", "logicle", "arcsinh"), 
                               scaleTo = NULL, 
                               q = 0.05,
                               l_w = 0.1, l_t = 4000, l_m = 4.5, l_a = 0,
                               a_a = 1, a_b = 1, a_c =0) {
    
    ## load FCS files
    name <- sub(".fcs", "", basename(fcsFile))
    if (verbose) {
        fcs <- read.FCS(fcsFile, transformation = FALSE)
    } else {
        fcs <- suppressWarnings(read.FCS(fcsFile, transformation = FALSE))
    }
    
    ## compensation
    if(is.matrix(comp)){
        fcs <- applyComp(fcs, comp)
        cat("    Compensation is applied on", fcsFile, "\n")
    }else if(isTRUE(comp)) {
        if(!is.null(fcs@description$SPILL)) {
            fcs <- applyComp(fcs, fcs@description[["SPILL"]])
            cat("    Compensation is applied on ", fcsFile, "\n")
        }else if(!is.null(fcs@description$SPILLOVER)) {
            fcs <- applyComp(fcs, fcs@description[["SPILLOVER"]])
            cat("    Compensation is applied on ", fcsFile, "\n")
        }else if(!is.null(fcs@description$COMP)) {
            fcs <- applyComp(fcs, fcs@description[["COMP"]])
            cat("    Compensation is applied on ", fcsFile, "\n")
        }else{
            warning("Cannot find compensation matrix in the FCS files!
                    Please CHECK the keyword of 'SPILL', 'SPILLOVER', or 'COMP'
                    in the FCS file and make sure it stores the compensation matrix.")
        }
    }

    ## match marker names to get marker ID, use all if NULL 
    pd <- fcs@parameters@data
    if (!(is.null(markers))) {
        if(is.numeric(markers)){
            if(all(markers %in% 1:ncol(fcs@exprs))){
                marker_id <- markers
            }else{
                stop("Marker ID out of range!")
            }
        }else if(is.character(markers)){
            right_marker <- markers %in% pd$desc || markers %in% pd$name
            if (!(right_marker)) {
                stop("\n Selected marker(s) is not in the input fcs files \n please check your selected markers! \n")
            } else {
                desc_id <- match(markers, pd$desc)
                name_id <- match(markers, pd$name)
                mids <- c(desc_id, name_id)
                marker_id <- unique(mids[!is.na(mids)])
            }
        }else{
            stop("Sorry, input markers cannot be recognized!")
        }
    }else{
        marker_id <- 1:ncol(fcs@exprs)
    }
   
    ## exprs transformation
    transformMethod <- match.arg(transformMethod)
    switch(transformMethod,
           cytofAsinh = {
               exprs <- apply(as.matrix(fcs@exprs[, marker_id]), 2, cytofAsinh)
           },
           autoLgcl = {
               trans <- autoLgcl(fcs, channels = colnames(fcs@exprs)[marker_id], q = q)
               transformed <- flowCore::transform(fcs, trans)
               exprs <- transformed@exprs[, marker_id]
           },
           logicle = {
               trans <- flowCore::logicleTransform(w = l_w, t = l_t, m = l_m, a = l_a)
               exprs <- apply(fcs@exprs[, marker_id], 2, trans)
           },
           arcsinh = {
               trans <- flowCore::arcsinhTransform(a = a_a, b = a_b, c = a_c)
               exprs <- apply(fcs@exprs[, marker_id], 2, trans)
           })
    
    ## rescale data
    if (!is.null(scaleTo)) {
        exprs <- apply(exprs, 2, function(x) scaleData(x, scaleTo))
    }
     
    ## add rownames and colnames   
    col_names <- paste0(pd$name, "<", pd$desc,">")
    colnames(exprs) <- col_names[marker_id]
    row.names(exprs) <- paste(name, 1:nrow(exprs), sep = "_")
    
    return(exprs)
}


#' apply compensation on the FCS expression data
#' 
#' @param fcs FCS file.
#' @param compMatrix Compensation matrix.
#' @noRd
#' @return compensated expression value
applyComp <- function(fcs, compMatrix) {
    comp_fcs <- compensate(fcs, compMatrix)
}

#' rescale the data
#' 
#' @param x data.
#' @param range The range of the data.
#' @noRd
#' @return scaled data
scaleData <- function(x, range = c(0, 4.5)) {
    (x - min(x))/(max(x) - min(x)) * (range[2] - range[1]) + range[1]
}


#' Noise reduced arsinh transformation 
#' 
#' Inverse hyperbolic sine transformation (arsinh) with a cofactor of 5, reduce noise from negative values
#' Adopted from Plos Comp reviewer
#' 
#' @param value A vector of numeric values.
#' @param cofactor Cofactor for asinh transformation, default 5 for mass cytometry data.
#' @noRd
#' @return transformed value
cytofAsinh <- function(value, cofactor = 5) {
    value <- value-1
    loID <- which(value < 0)
    if(length(loID) > 0)
        value[loID] <- rnorm(length(loID), mean = 0, sd = 0.01)
    value <- value / cofactor
    value <- asinh(value) # value <- log(value + sqrt(value^2 + 1))
    return(value)
}


#' a modified version of "estimateLogicle" from flowCore
#' 
#' 
#' @param x Data.
#' @param channels Channel names.
#' @param m Para m.
#' @param q Para q.
#' @importFrom methods is
#' @noRd
#' @return a list of transformations
autoLgcl <- function(x, channels, m = 4.5, q = 0.05) {
    if (!is(x, "flowFrame")) 
        stop("x has to be an object of class \"flowFrame\"")
    if (missing(channels)) 
        stop("Please specify the channels to be logicle transformed")
    indx <- channels %in% colnames(x@exprs)
    if (!all(indx)) 
        stop(paste("Channels", channels[!indx], "were not found in the FCS file.\n ", 
            sep = " "))
    rng <- range(x)
    trans <- lapply(channels, function(p) {
        lgclParaEstimate(x@exprs[, p], p = p, m = m, q = q)
    })
    transformList(channels, trans)
}

lgclParaEstimate <- function(data, p, m = 4.5, q = 0.05, type = "instrument") {
    t <- max(data)
    ndata <- data[data < 0]
    w <- 0
    a <- 0
    transId <- paste(p, "autolgclTransform", sep = "_")
    
    if (missing(m)) {
        if (type == "instrument") 
            m <- 4.5
        else m <- log10(t) + 1
    }
    
    if (length(ndata)) {
        r <- .Machine$double.eps + quantile(ndata, q)
        if (10^m * abs(r) <= t) {
            w <- 0  ## this special check to avoid failure of negative w
        } else {
            w <- (m - log10(t/abs(r)))/2
            if(w>2) {
                w <- 0.5
                #cat("w is coerced to 0.5\n")
            }
        }
    }
    # cat(paste('sign_autoLgcl parameters:', ' w=', w, ' t=',t,'
    # m=',m,' a=',a,'\n', sep = ''))
    logicleTransform(transformationId = transId, 
                     w = w, t = t, m = m, a = a)
} 
