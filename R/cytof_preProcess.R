#' merge the transformed expression data of FCS file(s) of selected markers
#' 
#' Apply logicle transformation of selected markers of each FCS file, auto logicle
#' transformation and fixed logicle transformation are provided, then mutilple 
#' FCS files are merged using method \code{all}, \code{min}, \code{fixed} or \code{ceil}
#' 
#' @param fcsFiles the input fcsFiles (usually more than 1 file)
#' @param comp Boolean tells if do compensation
#' @param verbose Boolean
#' @param markers Selected markers for analysis, either from names or from description
#' @param transformationMethod transformationMethod transformation method, three logicle transformation methods includes: \code{auto}, \code{sign_auto} or \code{fixed} for FCM data, and \code{arcsin} for CyTOF data
#' @param scaleTo scale the expression to same scale, default is NULL, should be a vector of two numbers if scale
#' @param w Linearization width in asymptotic decades
#' @param t Top of the scale data value
#' @param m Full width of the transformed display in asymptotic decades
#' @param a Additional negative range to be included in the display in asymptotic decades
#' @param q quantile of negative values removed for auto w estimation, default is 0.05
#' @param mergeMethod merge method for mutiple FCS expression data, default is all
#' @param fixedNum the fixed number of cells for merging multiple FCSs
#' @return Merged FCS expression data matrix of selected markers with logicle transformation
#' @export
#' @examples
#' d<-system.file('extdata',package='cytofkit')
#' fcsFile <- list.files(d,pattern='.fcs$',full=TRUE)
#' merged <- fcs_trans_merge(fcsFile)
fcs_trans_merge <- function(fcsFiles, comp = FALSE, verbose = FALSE, 
    markers = NULL, transformationMethod = "arcsin", scaleTo = NULL, w = 0.1, t = 4000, 
    m = 4.5, a = 0, q = 0.05, mergeMethod = "ceil", fixedNum = 10000) {
    
    exprsL <- mapply(fcs_trans, fcsFiles, MoreArgs = list(comp = comp, 
        verbose = verbose, markers = markers, transformationMethod = transformationMethod, 
        scaleTo = scaleTo, w = w, t = t, m = m, a = a, q = q), SIMPLIFY = FALSE)
    
    if (mergeMethod == "all") {
        merged <- do.call(rbind, exprsL)
    } else if (mergeMethod == "min") {
        minSize <- min(sapply(exprsL, nrow))
        mergeFunc <- function(x) {
            x[sample(nrow(x), size = minSize, replace = FALSE), ]
        }
        merged <- do.call(rbind, lapply(exprsL, mergeFunc))
    } else if (mergeMethod == "fixed") {
        mergeFunc <- function(x) {
            x[sample(nrow(x), size = fixedNum, replace = ifelse(nrow(x) < 
                fixedNum, TRUE, FALSE)), ]
        }
        merged <- do.call(rbind, lapply(exprsL, mergeFunc))
    } else if (mergeMethod == "ceil") {
        mergeFunc <- function(x) {
            if (nrow(x) < fixedNum) {
                x
            } else {
                x[sample(nrow(x), size = fixedNum, replace = FALSE), ]
            }
        }
        merged <- do.call(rbind, lapply(exprsL, mergeFunc))
    }
    return(merged)
}


#' Logicle transformation of the FCS data
#' 
#' Read the FCS expresssion data and apply Logicle transformation 
#' 
#' @param fcsFile The name of the FCS file
#' @param comp Boolean tells if do compensation
#' @param verbose Boolean
#' @param markers Selected markers for analysis, either from names or from description
#' @param transformationMethod transformationMethod transformation method, three logicle transformation methods includes: \code{auto}, \code{sign_auto} or \code{fixed} for FCM data, and \code{arcsin} for CyTOF data
#' @param scaleTo scale the expression to same scale, default is NULL, should be a vector of two numbers if scale
#' @param w Linearization width in asymptotic decades
#' @param t Top of the scale data value
#' @param m Full width of the transformed display in asymptotic decades
#' @param a Additional negative range to be included in the display in asymptotic decades
#' @param q quantile of negative values removed for auto w estimation, default is 0.05
#' @return The logicle transformend expression data matrix of selected markers
#' @importFrom flowCore read.FCS compensate estimateLogicle logicleTransform parameters transformList
#' @importMethodsFrom flowCore transform
#' @importClassesFrom flowCore transformList
#' @export
#' @examples
#' d<-system.file('extdata',package='cytofkit')
#' fcsFile <- list.files(d,pattern='.fcs$',full=TRUE)
#' transformed <- fcs_trans(fcsFile)
fcs_trans <- function(fcsFile, comp = FALSE, verbose = FALSE, 
                     markers = NULL, transformationMethod = "arcsin", scaleTo = NULL,
                     w = 0.1, t = 4000, m = 4.5, a = 0, q = 0.05) {
        
        ## load FCS data
        name <- sub(".fcs", "", basename(fcsFile))
        if (verbose) {
                fcs <- read.FCS(fcsFile)
        } else {
                fcs <- suppressWarnings(read.FCS(fcsFile))
        }
        
        ## compensation if provided in the data
        apply.comp <- function(fcs, keyword) {
                comp_fcs <- compensate(fcs, fcs@description[[keyword]])
        }
        if (comp == TRUE) {
                if (comp && !is.null(fcs@description$SPILL)) {
                        fcs <- apply.comp(fcs, "SPILL")
                } else if (comp && !is.null(fcs@description$SPILLOVER)) {
                        fcs <- apply.comp(fcs, "SPILLOVER")
                } else if (comp && !is.null(fcs@description$COMP)) {
                        fcs <- apply.comp(fcs, "COMP")
                }
        }

        pd <- fcs@parameters@data
        
        ## marker check, allow mix marker names
        if (!(is.null(markers))) {
                right_marker <- markers %in% pd$desc || markers %in% pd$name
                if (!(right_marker)) {
                        stop("\n Selected marker(s) is not in the input fcs files \n please check your selected markers! \n")
                } else {
                        desc_id <- match(markers, pd$desc)
                        name_id <- match(markers, pd$name)
                        mids <- c(desc_id, name_id)
                        marker_id <- unique(mids[!is.na(mids)])
                }
        } else {
                marker_id <- 1:ncol(fcs@exprs)
        }
        
        ## logicle transformation
        if (transformationMethod == "auto") {
                lgcl <- estimateLogicle(fcs, channels = colnames(exprs(fcs))[marker_id])
                lgcl_transformed <- transform(fcs, lgcl)
                exprs <- lgcl_transformed@exprs
        } else if (transformationMethod == "sign_auto") {
                lgcl <- sign_auto(fcs, channels = colnames(fcs@exprs)[marker_id])
                lgcl_transformed <- transform(fcs, lgcl)
                exprs <- lgcl_transformed@exprs
        } else if (transformationMethod == "fixed") {
                lgcl <- logicleTransform(w = w, t = t, m = m, a = a)
                exprs_raw <- fcs@exprs
                exprs <- apply(exprs_raw, 2, lgcl)
        } else if (transformationMethod == "arcsin") {
                a <- 1
                b <- 1
                c <- 0
                exprs <- fcs@exprs
                for(m in colnames(exprs)[marker_id]){
                        exprs[,m] <- asinh(a + b * exprs[,m]) + c
                }
        }
        
        ## scale the range of data
        if (!is.null(scaleTo)){
                exprs <- apply(exprs[, marker_id], 2, function(x) {scaleRange(x, scaleTo)})
        }else{
                exprs <- exprs[, marker_id]
        }
        
        col_names <- pd$desc
        col_names[is.na(pd$desc)] <- pd$name[is.na(pd$desc)]
        colnames(exprs) <- col_names[marker_id]
        row.names(exprs) <- paste(name, 1:nrow(exprs), sep = "_")
        return(exprs)
} 

scaleRange <- function(x, range = c(0,4.5)){
        (x-min(x))/(max(x)-min(x)) * (range[2] - range[1]) + range[1]
}

sign_auto <- function (x, channels, m = 4.5, q = 0.05) {
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
                sign_autoLgcl(x@exprs[ ,p], m = m, q = q)
        })
        transformList(channels, trans)
}

sign_autoLgcl <- function(data, m=m, q=q) {
        t <- max(data)
        ndata <- data[data<0]
        w <- 0
        a <- 0
        if(length(ndata)) {
                r <- .Machine$double.eps + quantile(ndata, q)
                if (10^m * abs(r) <= t){
                        w <-  0
                }else{
                        w <- (m-log10(t/abs(r))) / 2
                }      
        } 
        #cat(paste("sign_autoLgcl parameters:", " w=", w, " t=",t," m=",m," a=",a,"\n", sep = ""))
        logicleTransform( w=w, t=t, m=m, a=a)
}
