#' The user friendly GUI for function \code{cytof_tsne_densvm}
#' 
#' This GUI provides an easy way for CyToF data analysis using \code{cytofkit} package. All parameters 
#' for running 'cytof_tsne_densvm' were integrated in this GUI, each parameter has help button 
#' in the GUI to help user get details of the information of each parameter, and launch the 
#' \code{cytof_tsne_densvm} analysis after submitting.
#' 
#' @author Chen Hao
#' @return the GUI for the main function \code{cytof_tsne_densvm}
#' @import tcltk
#' @export
#' @seealso \code{\link{cytof_tsne_densvm}}, \code{\link{cytofkit}}
#' @references \url{http://signbioinfo.github.io/cytofkit/}
#' @examples
#' #cytof_tsne_densvm_GUI()  # remove the comment hash to run
#' 
cytof_tsne_densvm_GUI <- function() {
    
    ## input parameters
    cur_dir <- getwd()
    mergeMethod_array <- c("all", "min", "ceil", "fixed")
    fcsFiles <- ""
    
    rawFCSdir <- tclVar(cur_dir)
    fcsFile <- tclVar("")
    resDir <- tclVar(cur_dir)
    baseName <- tclVar("cytof")
    mergeMethod <- tclVar(mergeMethod_array[3])
    fixedNum <- tclVar("10000")
    lgclMethod <- tclVar("fixed")
    paraFile <- tclVar("")
    ifTransform <- tclVar("TRUE")
    transformMethod <- tclVar("tsne")
    ifCluster <- tclVar("TRUE")
    ret_var <- tclVar("")
    ## visualizationMethods <- tclVar("tsne")
    vizMethods <- c("pca","isomap","tsne")
    vizSelect <- c()
    i <- 1
    j <- 3
    while (i<=j){
            aux<-paste("vizSelect",i,sep="")
            vizSelect<-c(vizSelect,as.character(aux))
            tclvalue(vizSelect[i]) <- "0"
            i<-i+1
    }
    tclvalue(vizSelect[3]) <- "1"  # default tsne 
    
    ## button functions
    reset_rawFCS_dir <- function() {
        rawFCS_dir <- ""
        rawFCS_dir <- tclvalue(tkchooseDirectory(title = "chose your rawFCS dircetory ..."))
        if (rawFCS_dir != "") {
            tclvalue(rawFCSdir) <- rawFCS_dir
            tclvalue(resDir) <- rawFCS_dir
        }
    }
    
    reset_res_dir <- function() {
        res_dir <- ""
        res_dir <- tclvalue(tkchooseDirectory(title = "chose your result dircetory ..."))
        if (res_dir != "") {
            tclvalue(resDir) <- res_dir
        }
    }
    
    reset_fcs_data <- function() {
            fnames <- ""
            fnames <- tk_choose.files(default = paste(tclvalue(rawFCSdir), "fcs", sep = .Platform$file.sep), 
                                      caption = "Select FCS files", multi = TRUE, 
                                      filters = matrix(c("{fcs files}", "{.fcs}"), 1, 2), index = 1)
            if (length(fnames) >= 1) {
                    fnames <- fnames[!(grepl(paste0(.Platform$file.sep, ".fcs$"), fnames))]  # remove empty .fcs files
                    tclvalue(fcsFile) <- paste(fnames, collapse = "}{")
            }
    }
    
    reset_para_data <- function() {
        
        if (tclvalue(fcsFile) == "" && tclvalue(rawFCSdir) == 
            "") {
            tkmessageBox(title = "cytofkit: an integrated analysis pipeline for mass cytometry data", 
                message = "Please input your \"rawFCSdirectory\" or \"fcsFile\".", 
                icon = "info", type = "ok")
        }
        if (tclvalue(fcsFile) != "") {
            # fcsFiles <- sub('[{}]', '', strsplit(tclvalue(fcsFile), '}
            # {', fixed = TRUE)[[1]])
            fcsFiles <- strsplit(tclvalue(fcsFile), "}{", fixed = TRUE)[[1]]
        }
        fcsDir <- tclvalue(rawFCSdir)
        selectMarkers <- getParameters_GUI(fcsFiles, fcsDir)
        
        if (length(selectMarkers) > 0) {
            tclvalue(paraFile) <- paste(selectMarkers, collapse = "}{")
        }
    }
    
    reset_num2null <- function() {
        tclvalue(fixedNum) <- "NULL"
    }
    
    reset_num2any <- function() {
        tclvalue(fixedNum) <- "1000"
    }
    
    rawFCSdir_help <- function() {
        tkmessageBox(title = "rawFCSdir", message = "The directory that contains fcs files.", 
            icon = "info", type = "ok")
    }
    
    fcsFile_help <- function() {
        tkmessageBox(title = "fcsFile", message = "The fcs files to be analyzed. One or multiple fcs files are allowed. When multiple fcs files are selected, cells from each fcs file are combined for analysis.", 
            icon = "info", type = "ok")
    }
    
    para_help <- function() {
        tkmessageBox(title = "paraFile", message = "Select the list of makers to be used for analysis.", 
            icon = "info", type = "ok")
    }
    
    baseName_help <- function() {
        tkmessageBox(title = "baseName", message = "A prefix that will be added to the names of result files.", 
            icon = "info", type = "ok")
    }
    
    mergeMethod_help <- function() {
        tkmessageBox(title = "mergeMethod", message = "When multiple fcs files are selected, cells can be combined using one of the four different methods including \"ceil\",\"all\", \"min\",\"fixed\". \n\n\"ceil\" (the default option): up to a fixed number (specified by fixedNum) of cells are sampled without replacement from each fcs file and combined for analysis. \n\n\"all\": all cells from each fcs file are combined for analysis. \n\n\"min\": The minimum number of cells among all the selected fcs files are sampled from each fcs file and combined for analysis. \n\n\"fixed\": a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each fcs file and combined for analysis.", 
            icon = "info", type = "ok")
    }
    
    fixedNum_help <- function() {
        tkmessageBox(title = "fixedNum", message = "Up to fixedNum of cells from each fcs file are used for analysis.", 
            icon = "info", type = "ok")
    }
    
    reset_ifTransform <- function() {
        tclvalue(ifTransform) <- "FALSE"
    }
    
    transformMethod_help <- function() {
        tkmessageBox(title = "transformMethod", message = "The method used for dimensionality reduction, including \"pca\", \"isomap\", \"tsne\". \n\nIf \"NULL\" was selected, no dimension reduction will be performed.", 
            icon = "info", type = "ok")
    }
    
    ifCluster_help <- function() {
        tkmessageBox(title = "ifCluster", message = "If clustering analysis to be execute. The clustering analysis will takes couple of hours.", 
            icon = "warning", type = "ok")
    }
    
    visualizationMethods_help <- function(){
            tkmessageBox(title = "visualizationMethods", message = "The method used for visualizing the clustering results(s), multiple selections are allowed. Including \"pca\", \"isomap\", \"tsne\". \n\nWARNING: \"tsne\" is the default selection, \"isomap\" may take long time.", 
                         icon = "info", type = "ok")
    }
    
    resDir_help <- function() {
        tkmessageBox(title = "resDir", message = "The directory where result files will be generated", 
            icon = "info", type = "ok")
    }
    
    reset <- function() {
        tclvalue(rawFCSdir) = cur_dir
        tclvalue(fcsFile) = ""
        tclvalue(resDir) = cur_dir
        tclvalue(baseName) = "cytof"
        tclvalue(mergeMethod) = mergeMethod_array[3]
        tclvalue(fixedNum) = "10000"
        tclvalue(lgclMethod) = "fixed"
        tclvalue(paraFile) = ""
        tclvalue(ifTransform) = "TRUE"
        tclvalue(transformMethod) = "tsne"
        tclvalue(vizSelect[3]) <- "1"
        tclvalue(ifCluster) = "TRUE"
    }
    
    submit <- function() {
        has_error = FALSE
        if (tclvalue(paraFile) == "") {
            tkmessageBox(title = "cytofkit: an integrated analysis pipeline for mass cytometry data", 
                message = "Please select the markers for your analysis.", 
                icon = "info", type = "ok")
            has_error = TRUE
        }
        
        if (has_error == FALSE) {
            tclvalue(ret_var) <- "OK"
            tkdestroy(tt)
        }
    }
    
    quit <- function() {
        # q(save = 'no')
        tkdestroy(tt)
    }
    
    ## build the GUI widgets
    
    ## head line
    tt <- tktoplevel(borderwidth = 20)
    tkwm.title(tt, "cytofkit: an integrated analysis pipeline for mass cytometry data")
    
    box_length <- 55
    cell_width <- 3
    bt_width <- 8
    hb_width <- 2
    
    ## rawFCSdir
    rawFCSdir_label <- tklabel(tt, text = "Raw_FCS_directory :")
    rawFCSdir_entry <- tkentry(tt, textvariable = rawFCSdir, 
        width = box_length)
    rawFCSdir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, 
        command = reset_rawFCS_dir)
    rawFCSdir_hBut <- tkbutton(tt, text = "?", width = hb_width, 
        command = rawFCSdir_help)
    
    ## resDir
    resDir_label <- tklabel(tt, text = "Result_directory :")
    resDir_entry <- tkentry(tt, textvariable = resDir, width = box_length)
    resDir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, 
        command = reset_res_dir)
    resDir_hBut <- tkbutton(tt, text = "?", width = hb_width, 
        command = resDir_help)
    
    ## fcsFile
    fcsFile_label <- tklabel(tt, text = "FCS_File(s) :")
    fcsFile_entry <- tkentry(tt, textvariable = fcsFile, width = box_length)
    fcsFile_button <- tkbutton(tt, text = " Select... ", width = bt_width, 
        command = reset_fcs_data)
    fcsFile_hBut <- tkbutton(tt, text = "?", width = hb_width, 
        command = fcsFile_help)
    
    ## paraFile
    paraFile_label <- tklabel(tt, text = "Parameters :")
    paraFile_entry <- tkentry(tt, textvariable = paraFile, width = box_length)
    paraFile_button <- tkbutton(tt, text = " Select... ", width = bt_width, 
        command = reset_para_data)
    paraFile_hBut <- tkbutton(tt, text = "?", width = hb_width, 
        command = para_help)
    
    ## baseName
    baseName_label <- tklabel(tt, text = "BaseName :")
    baseName_entry <- tkentry(tt, textvariable = baseName, width = box_length)
    baseName_hBut <- tkbutton(tt, text = "?", width = hb_width, 
        command = baseName_help)
    
    ## mergeMethod && fixedNum
    mergeMethod_label <- tklabel(tt, text = "Merge_method :")
    mergeMethod_hBut <- tkbutton(tt, text = "?", width = hb_width, 
        command = mergeMethod_help)
    merge_method_rbuts <- tkframe(tt)
    tkpack(tklabel(merge_method_rbuts, text = ""), side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethod_array[1], 
        variable = mergeMethod, value = mergeMethod_array[1], 
        command = reset_num2null), side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethod_array[2], 
        variable = mergeMethod, value = mergeMethod_array[2], 
        command = reset_num2null), side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethod_array[3], 
        variable = mergeMethod, value = mergeMethod_array[3], 
        command = reset_num2any), side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethod_array[4], 
        variable = mergeMethod, value = mergeMethod_array[4], 
        command = reset_num2any), side = "left")
    tkpack(tklabel(merge_method_rbuts, text = "                  Fixed_num :"), 
        side = "left")
    tkpack(tkentry(merge_method_rbuts, textvariable = fixedNum, 
        width = 8), side = "right")
    fixedNum_hBut <- tkbutton(tt, text = "?", width = hb_width, 
        command = fixedNum_help)
    
    
    ## transformMethod
    transformMethod_label <- tklabel(tt, text = "Dimention_reduction :")
    transformMethod_hBut <- tkbutton(tt, text = "?", width = hb_width, 
        command = transformMethod_help)
    transformMethod_rbuts <- tkframe(tt)
    tkpack(tklabel(transformMethod_rbuts, text = ""), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = "pca", 
        variable = transformMethod, value = "pca"), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = "isomap", 
        variable = transformMethod, value = "isomap"), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = "tsne", 
        variable = transformMethod, value = "tsne"), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = "NULL", 
        variable = transformMethod, command = reset_ifTransform, 
        value = "NULL"), side = "left")
    
    ## visualizationMethods
    visualizationMethods_label <- tklabel(tt, text = "Visualization_method(s) :")
    visualizationMethods_hBut <- tkbutton(tt, text = "?", width = hb_width, 
                                          command = visualizationMethods_help)
    visualizationMethods_cbuts <- tkframe(tt)
    tkpack(tklabel(visualizationMethods_cbuts, text = ""), side = "left")
    i<-1
    while (i<=j){
            t<-tkcheckbutton(visualizationMethods_cbuts,text=vizMethods[i],variable=eval(vizSelect[i]))
            tkpack(t, side = "left")
            i<-i+1
    }
    
    ## ifCluster
    ifCluster_label <- tklabel(tt, text = "If_cluster :")
    ifCluster_hBut <- tkbutton(tt, text = "?", width = hb_width, 
        command = ifCluster_help)
    ifCluster_rbuts <- tkframe(tt)
    tkpack(tklabel(ifCluster_rbuts, text = ""), side = "left")
    tkpack(tkradiobutton(ifCluster_rbuts, text = "YES", variable = ifCluster, 
        value = "TRUE"), side = "left")
    tkpack(tkradiobutton(ifCluster_rbuts, text = "NO", variable = ifCluster, 
        value = "FALSE"), side = "left")
    
    ## submit / reset / quit
    submit_button <- tkbutton(tt, text = "Submit", command = submit)
    reset_button <- tkbutton(tt, text = "Reset", command = reset)
    quit_button <- tkbutton(tt, text = "Quit", command = quit)
    
    ## display GUI
    tkgrid(rawFCSdir_label, rawFCSdir_hBut, rawFCSdir_entry, 
        rawFCSdir_button, padx = cell_width)
    tkgrid.configure(rawFCSdir_label, rawFCSdir_entry, rawFCSdir_button, 
        sticky = "e")
    tkgrid.configure(rawFCSdir_hBut, sticky = "w")
    
    tkgrid(fcsFile_label, fcsFile_hBut, fcsFile_entry, fcsFile_button, 
        padx = cell_width)
    tkgrid.configure(fcsFile_label, fcsFile_entry, fcsFile_button, 
        sticky = "e")
    tkgrid.configure(fcsFile_hBut, sticky = "w")
    
    tkgrid(paraFile_label, paraFile_hBut, paraFile_entry, paraFile_button, 
        padx = cell_width)
    tkgrid.configure(paraFile_label, paraFile_entry, paraFile_button, 
        sticky = "e")
    tkgrid.configure(paraFile_hBut, sticky = "w")
    
    tkgrid(resDir_label, resDir_hBut, resDir_entry, resDir_button, 
        padx = cell_width)
    tkgrid.configure(resDir_label, resDir_entry, resDir_button, 
        sticky = "e")
    tkgrid.configure(resDir_hBut, sticky = "w")
    
    tkgrid(baseName_label, baseName_hBut, baseName_entry, padx = cell_width)
    tkgrid.configure(baseName_label, baseName_entry, sticky = "e")
    tkgrid.configure(baseName_hBut, sticky = "w")
    
    tkgrid(mergeMethod_label, mergeMethod_hBut, merge_method_rbuts, 
        fixedNum_hBut, padx = cell_width)
    tkgrid.configure(mergeMethod_label, sticky = "e")
    tkgrid.configure(merge_method_rbuts, sticky = "w")
    tkgrid.configure(mergeMethod_hBut, sticky = "w")
    tkgrid.configure(fixedNum_hBut, sticky = "w")
    
    tkgrid(transformMethod_label, transformMethod_hBut, transformMethod_rbuts, 
        padx = cell_width)
    tkgrid.configure(transformMethod_label, sticky = "e")
    tkgrid.configure(transformMethod_rbuts, sticky = "w")
    tkgrid.configure(transformMethod_hBut, sticky = "w")
    
    tkgrid(ifCluster_label, ifCluster_hBut, ifCluster_rbuts, 
        padx = cell_width)
    tkgrid.configure(ifCluster_label, sticky = "e")
    tkgrid.configure(ifCluster_rbuts, sticky = "w")
    tkgrid.configure(ifCluster_hBut, sticky = "w")
    
    tkgrid(visualizationMethods_label, visualizationMethods_hBut, visualizationMethods_cbuts, 
           padx = cell_width)
    tkgrid.configure(visualizationMethods_label, sticky = "e")
    tkgrid.configure(visualizationMethods_cbuts, sticky = "w")
    tkgrid.configure(visualizationMethods_hBut, sticky = "w")
    
    tkgrid(tklabel(tt, text = "\n"), padx = cell_width)  # leave blank line
    
    tkgrid(tklabel(tt, text = ""), reset_button, submit_button, 
        quit_button, padx = cell_width)
    tkgrid.configure(reset_button, sticky = "e")
    tkgrid.configure(quit_button, sticky = "w")
    
    tkwait.window(tt)
    
    if (tclvalue(ret_var) != "OK") {
        stop("Analysis is cancelled.")
    }
    
    ## return the inputs fcsFiles <- sub('[{}]', '',
    ## strsplit(tclvalue(fcsFile), '} {', fixed = TRUE)[[1]])
    fcsFiles <- strsplit(tclvalue(fcsFile), "}{", fixed = TRUE)[[1]]
    parameters <- strsplit(tclvalue(paraFile), "}{", fixed = TRUE)[[1]]
    
    ## get visualization options
    vizCheck <- c()
    i<-1
    while (i<=j){
            v <- as.numeric(tclvalue(vizSelect[i])) > 0
            vizCheck <- c(vizCheck, v)
            i <- i + 1
    }
    
    inputs <- list()
    inputs[["rawFCSdir"]] <- tclvalue(rawFCSdir)
    inputs[["fcsFile"]] <- fcsFiles
    inputs[["resDir"]] <- tclvalue(resDir)
    inputs[["baseName"]] <- tclvalue(baseName)
    inputs[["mergeMethod"]] <- tclvalue(mergeMethod)
    inputs[["fixedNum"]] <- suppressWarnings(as.numeric(tclvalue(fixedNum)))
    inputs[["lgclMethod"]] <- tclvalue(lgclMethod)
    inputs[["paraFile"]] <- parameters
    inputs[["ifTransform"]] <- as.logical(tclvalue(ifTransform))
    inputs[["transformMethod"]] <- tclvalue(transformMethod)
    inputs[["ifCluster"]] <- as.logical(tclvalue(ifCluster))
    inputs[["visualizationMethods"]] <- vizMethods[vizCheck]
            
                    
           
    ## return(inputs)
    
    ## pass the parameters and run the cytof_tsne_densvm function
    cytof_tsne_densvm(rawFCSdir = inputs[["rawFCSdir"]], fcsFile = inputs[["fcsFile"]], 
        resDir = inputs[["resDir"]], baseName = inputs[["baseName"]], 
        mergeMethod = inputs[["mergeMethod"]], fixedNum = inputs[["fixedNum"]], 
        transformationMethod = inputs[["lgclMethod"]], para = inputs[["paraFile"]], 
        ifTransform = inputs[["ifTransform"]], dimReductionMethod = inputs[["transformMethod"]], 
        ifCluster = inputs[["ifCluster"]], visualizationMethods = inputs[["visualizationMethods"]])
}


getParameters_GUI <- function(fcsFile, rawFCSdir) {
    
    if (is.null(fcsFile)) {
        fcsFile <- list.files(path = rawFCSdir, pattern = ".fcs$")
    }
    
    fcs <- suppressWarnings(read.FCS(fcsFile[1]))
    pd <- fcs@parameters@data
    markers <- paste("<", pd$name, ">:", pd$desc, sep = "")
    channels <- pd$name
 
    if (length(markers) == 0) {
        stop("No markers found in the FCS file!")
    }
    
    # GUI
    markerChoice <- tclVar("")
    mm <- tktoplevel()
    tkwm.title(mm, "cytofkit: marker selection")
    scr <- tkscrollbar(mm, repeatinterval = 5, command = function(...) tkyview(tl, 
        ...))
    tl <- tklistbox(mm, height = 20, selectmode = "multiple", 
        yscrollcommand = function(...) tkset(scr, ...), background = "white")
    OnOK <- function() {
        tclvalue(markerChoice) <- paste(markers[as.numeric(tkcurselection(tl)) + 
            1], collapse = "}{")
        tkdestroy(mm)
    }
    OK.but <- tkbutton(mm, text = " OK ", command = OnOK)
    tkgrid(tklabel(mm, text = "Please select your markers:"))
    tkgrid(tl, scr)
    tkgrid.configure(scr, rowspan = 4, sticky = "nsw")
    for (i in (1:length(markers))) {
        tkinsert(tl, "end", markers[i])
    }
    tkgrid(OK.but)
    tkwait.window(mm)
    
    # return parameters
    paras <- strsplit(tclvalue(markerChoice), "}{", fixed = TRUE)[[1]]
    paras <- channels[match(paras, markers)]
    return(paras)
} 
