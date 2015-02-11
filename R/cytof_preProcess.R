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
#' @param lgclMethod Logicle transformation method, \code{auto}, \code{sign_auto} or \code{fixed}
#' @param w Linearization width in asymptotic decades
#' @param t Top of the scale data value
#' @param m Full width of the transformed display in asymptotic decades
#' @param a Additional negative range to be included in the display in asymptotic decades
#' @param mergeMethod merge method for mutiple FCS expression data, default is all
#' @param fixedNum the fixed number of cells for merging multiple FCSs
#' @return Merged FCS expression data matrix of selected markers with logicle transformation
#' @export
#' @examples
#' d<-system.file('extdata',package='cytofkit')
#' fcsFile <- list.files(d,pattern='.fcs$',full=TRUE)
#' merged <- fcs_lgcl_merge(fcsFile)
fcs_lgcl_merge <- function(fcsFiles, comp = FALSE, verbose = FALSE, 
    markers = NULL, lgclMethod = "fixed", w = 0.1, t = 4000, 
    m = 4.5, a = 0, mergeMethod = "ceil", fixedNum = 10000) {
    
    exprsL <- mapply(fcs_lgcl, fcsFiles, MoreArgs = list(comp = comp, 
        verbose = verbose, markers = markers, lgclMethod = lgclMethod, 
        w = w, t = t, m = m, a = a), SIMPLIFY = FALSE)
    
    if (mergeMethod == "all") {
        merged <- do.call(rbind, exprsL)
    } else if (mergeMethod == "min") {
        minSize <- min(sapply(exprsL, nrow))
        mergeFunc <- function(x) {
            x[sample(nrow(x), size = minSize, replace = FALSE), 
                ]
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
                x[sample(nrow(x), size = fixedNum, replace = FALSE), 
                  ]
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
#' @param lgclMethod Logicle transformation method, \code{auto}, \code{sign_auto} or \code{fixed}
#' @param w Linearization width in asymptotic decades
#' @param t Top of the scale data value
#' @param m Full width of the transformed display in asymptotic decades
#' @param a Additional negative range to be included in the display in asymptotic decades
#' @return The logicle transformend expression data matrix of selected markers
#' @importFrom flowCore read.FCS compensate estimateLogicle logicleTransform parameters transformList
#' @importMethodsFrom flowCore transform
#' @importClassesFrom flowCore transformList
#' @export
#' @examples
#' d<-system.file('extdata',package='cytofkit')
#' fcsFile <- list.files(d,pattern='.fcs$',full=TRUE)
#' transformed <- fcs_lgcl(fcsFile)
fcs_lgcl <- function(fcsFile, comp = FALSE, verbose = FALSE, 
    markers = NULL, lgclMethod = "fixed", w = 0.1, t = 4000, 
    m = 4.5, a = 0) {
    
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
            marker_id <- unique(desc_id, name_id)
        }
    } else {
        marker_id <- 1:ncol(fcs@exprs)
    }
    
    ## logicle transformation
    if (lgclMethod == "auto") {
        lgcl <- estimateLogicle(fcs, channels = colnames(exprs(fcs))[marker_id])
        lgcl_transformed <- transform(fcs, lgcl)
        exprs <- lgcl_transformed@exprs
    } else if (lgclMethod == "sign_auto") {
            lgcl <- sign_auto(fcs, channels = colnames(exprs(fcs))[marker_id])
            lgcl_transformed <- transform(fcs, lgcl)
            exprs <- lgcl_transformed@exprs
    } else if (lgclMethod == "fixed") {
        lgcl <- logicleTransform(w = w, t = t, m = m, a = a)
        exprs_raw <- fcs@exprs
        exprs <- apply(exprs_raw, 2, lgcl)
    }
    
    ## return the transformed expression data of selected markers
    col_names <- pd$desc
    col_names[is.na(pd$desc)] <- pd$name[is.na(pd$desc)]
    exprs <- exprs[, marker_id]
    colnames(exprs) <- col_names[marker_id]
    row.names(exprs) <- paste(name, 1:nrow(exprs), sep = "_")
    return(exprs)
} 


sign_auto <- function (x, channels, ...) {
        if (!is(x, "flowFrame")) 
                stop("x has to be an object of class \"flowFrame\"")
        if (missing(channels)) 
                stop("Please specify the channels to be logicle transformed")
        indx <- channels %in% colnames(x)
        if (!all(indx)) 
                stop(paste("Channels", channels[!indx], "were not found in x ", 
                           sep = " "))
        sign_autoLgcl <- function(data, m=4.5, q=0.05) {
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
                cat(paste("sign_autoLgcl parameters:", " w=", w, " t=",t," m=",m," a=",a,"\n", sep = ""))
                logicleTransform( w=w, t=t, m=m, a=a)
        }
        
        rng <- range(x)
        trans <- lapply(channels, function(p) {
                sign_autoLgcl(x)
        })
        transformList(channels, trans)
}
