#' cytofkit: an integrated analysis pipeline for mass cytometry data
#' 
#' This package is designed to facilitate the analysis workflow of mass cytometry data with 
#' automatic subset identification and population boundary detection. Both command line and 
#' a GUI are provided for runing the workflow easily.
#' 
#' This package integrates merge methods of multiple FCS files, dimension reduction (PCA, t-SNE and ISOMAP) 
#' with density-based clustering (DensVM) for rapid subset detection. Subset-clustering scatter plot 
#' and heat map will be generated for objective comparative analysis and statistical testing. This workflow can be 
#' easily done using the main function \code{\link{cytof_tsne_densvm}} or a GUI for the main function 
#' \code{\link{cytof_tsne_densvm_GUI}}.
#' 
#' Pre-processing
#' 
#' Using function \code{\link{fcs_trans_merge}}, one or multiple FCS files were imported via the *read.FCS* 
#' function in the *flowCore* package. Then logicle transformation was applied to the expression value 
#' of selected markers of each FCS file. Auto logicle transformation and fixed logicle transformation 
#' are provided, then mutilple FCS files are merged using method \code{all}, \code{min}, \code{fixed} 
#' or \code{ceil}.
#' 
#' Dimensionality reduction
#' 
#' Using function \code{\link{cytof_dimReduction}}, t-Distributed Stochastic Neighbor Embedding (\code{tsne}) 
#' is suggested for dimensionality reduction although we also provide methods like \code{isomap} and \code{pca}.
#' 
#' Cluster analysis using DensVM
#' 
#' Density-based clustering aided by support Vector Machine (\code{\link{densVM_cluster}}) are used to automate 
#' subset detection from the dimension-reducted map. By using DensVM, we are able to objectively assign every 
#' cell to an appropriate cluster.
#' 
#' Post-processing
#' 
#' Cluster results are annotated by using scatter plot and heatmap. Scatter plot visualize the cell points 
#' with colour indicating their assigned clusters and point shape representing their belonging samples
#' (\code{\link{cluster_plot}} and \code{\link{cluster_gridPlot}}). Cell events are also grouped by clusters 
#' and samples, and mean intensity values per cluster for every marker is calculated
#' (\code{\link{clust_mean_heatmap}} and \code{\link{clust_percentage_heatmap}}). 
#' Heat map visualizing the mean expression of every marker in every cluster is generated with no scaling on 
#' the row or column direction. Hierarchical clustering was generated using Euclidean distance and complete 
#' agglomeration method. We used the heat maps to interrogate marker expression to identify each cluster's 
#' defining markers. All intermediate files and the plots can be saved using the function \code{\link{cytof_write_results}}.
#' 
#' @examples
#' 
#' ## Run on GUI
#' #cytof_tsne_densvm_GUI()  # remove the hash symbol to launch the GUI
#' 
#' ## Run on command
#' dir <- system.file('extdata',package='cytofkit')
#' file <- list.files(dir, pattern='.fcs$', full=TRUE)
#' parameters <- list.files(dir, pattern='.txt$', full=TRUE)
#' ## remove the hash symbol to run the following command
#' #cytof_tsne_densvm(fcsFile = file, paraFile = parameters, rawFCSdir = dir, baseName = 'test')
#' 
#' ## Checking the vignettes for more details 
#' if(interactive()) browseVignettes(package = 'cytofkit')
#' 
#' @seealso \code{\link{cytof_tsne_densvm}}, \code{\link{cytof_tsne_densvm_GUI}}
#' @references \url{http://signbioinfo.github.io/cytofkit/}
#' @docType package
#' @name cytofkit
#' 
NULL


#' CyTOF data analysis for subpopulation detection 
#' 
#' \code{cytof_tsne_densvm} provides a workflow for one or multiple CyToF data analysis, 
#' including data preprocess with merging methods of multiple fcs file, logicle transformation, 
#' dimension reduction with PCA, isomap or tsne(default), and a kernal-based local maxima 
#' clustering combined with SVM for subpopulation detection. The intermediate results can be saved 
#' into seperate files and the cluster results can be visualized in heatmaps and scatter plots.
#' 
#' @param rawFCSdir the directory that contains fcs files to be analysed.
#' @param fcsFile a vector containing names of fcs files to be analyzed. One or multiple fcs files are allowed.
#' @param resDir the directory where result files will be generated.
#' @param baseName a prefix that will be added to the names of result files.
#' @param comp Boolean tells if do compensation. This will be applied to flow cytometry data.
#' @param verbose Boolean.
#' @param mergeMethod when multiple fcs files are selected, cells can be combined using 
#' one of the four different methods including \code{ceil}, \code{all}, \code{min}, \code{fixed}. 
#' The default option is \code{ceil}, up to a fixed number (specified by \code{fixedNum}) of cells are sampled 
#' without replacement from each fcs file and combined for analysis.
#' \code{all}: all cells from each fcs file are combined for analysis. 
#' \code{min}: The minimum number of cells among all the selected fcs files are sampled from each fcs file and combined for analysis. 
#' \code{fixed}: a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than 
#' fixedNum) from each fcs file and combined for analysis.
#' @param fixedNum up to fixedNum of cells from each fcs file are used for analysis.
#' @param transformationMethod transformation method, three logicle transformation methods includes: \code{auto}, \code{sign_auto} or \code{fixed} for FCM data, and \code{arcsin} for CyTOF data.
#' @param scaleTo scale the expression to same scale, default is NULL, should be a vector of two numbers if scale
#' @param q quantile of negative values removed for auto w estimation, default is 0.05
#' @param para the vector of selected makers. This can be provided in the \code{paraFile}.
#' @param paraFile a text file that specifies the list of makers to be used for analysis.
#' @param ifTransform a boolean to decide if dimensionality reduction will be performed. Default is TRUE.
#' @param dimReductionMethod the method used for dimensionality reduction, including \code{tsne}, \code{pca} and \code{isomap}.
#' @param ifCluster a boolean to determine if cluster will be conducted.
#' @param visualizationMethods the method(s) used for visualize the cluster data, multiple selection are accepted, including \code{tsne}, \code{pca} and \code{isomap}
#' @param writeResults if save the results, and the post-processing results including scatter plot, heatmap, and statistical results.
#' @param ... more arguments contral the logicle transformation
#' @return a list containing \code{transMergedExprs}, \code{transData} and \code{clustersRes}. If choose 'writeResults = TRUE', results will be saved into files under \code{resDir}
#' @author Chen Jinmiao 
#' @references \url{http://signbioinfo.github.io/cytofkit/}
#' @seealso \code{\link{cytofkit}}, \code{\link{cytof_tsne_densvm_GUI}}
#' @export
#' @examples
#' dir <- system.file('extdata',package='cytofkit')
#' file <- list.files(dir, pattern='.fcs$', full=TRUE)
#' parameters <- list.files(dir, pattern='.txt$', full=TRUE)
#' ## remove the hash symbol to run the following command
#' #cytof_tsne_densvm(fcsFile = file, paraFile = parameters, rawFCSdir = dir, baseName = 'test')   
cytof_tsne_densvm <- function(rawFCSdir = getwd(), fcsFile = NULL, 
    resDir = getwd(), baseName = "cytofkit_analysis", para = NULL, 
    paraFile = "./parameter.txt", comp = FALSE, verbose = FALSE, 
    transformationMethod = "arcsin", scaleTo = NULL, q = 0.05, mergeMethod = "ceil", fixedNum = 10000, 
    ifTransform = TRUE, dimReductionMethod = "tsne", ifCluster = TRUE,
    visualizationMethods = "tsne", writeResults = TRUE, ...) {
    
    ## para checking
    if (is.null(fcsFile))
        fcsFile <- list.files(path = rawFCSdir, pattern = ".fcs$", 
            full.names = TRUE)
    if (is.null(fcsFile) || length(fcsFile) < 1)
            stop("No FCS file selected!")
    if (is.null(para))
        para <- as.character(read.table(paraFile, sep = "\t", 
            header = TRUE)[, 1])
    if (is.null(para) || length(para) < 1)
            stop("no parameter selected!")
    if(!(mergeMethod %in% c("ceil", "all", "min", "fixed")))
            stop("wrong mergeMethod selected!")  
    if (!(transformationMethod %in% c("auto", "sign_auto", "fixed")))
            stop("wrong transformationMethod selected!")
    if(!(dimReductionMethod %in% c("tsne", "pca", "isomap")))
            stop("wrong dimReductionMethod selected!")
    if(!(all(visualizationMethods %in% c("tsne", "pca", "isomap"))))
            stop("wrong visualizationMethods selected")
    
    ## get transformed, combined, marker-filtered exprs data
    para <- sort(para)
    exprs <- fcs_trans_merge(fcsFile, comp = FALSE, verbose = FALSE, 
        markers = para, transformationMethod = transformationMethod, scaleTo = scaleTo, 
        q = q, mergeMethod = mergeMethod, fixedNum = fixedNum)
    
    ## dimension reduction
    transformed <- NULL
    if (ifTransform){
            transformed <- cytof_dimReduction(exprs, method = dimReductionMethod)      
    }else{
            transformed <- NULL
    }
        
    
    ## cluster
    cluster_output <- NULL
    if (ifCluster){
            cluster_output <- densVM_cluster(transformed, exprs)
    }else{
            cluster_output <- NULL
    }
        
    
    ## write results
    analysis_results <- list(transMergedExprs = exprs, transData = transformed, 
        clustersRes = cluster_output)
    if (writeResults == TRUE) {
        cytof_write_results(analysis_results, visualizationMethods, baseName, rawFCSdir, resDir)
    } else {
        return(analysis_results)
    }
} 
