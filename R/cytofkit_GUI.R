#' The user friendly GUI client for \code{cytofkit-package}
#' 
#' This GUI provides an easy way for CyToF data analysis using \code{cytofkit} package. 
#' Main parameters for running 'cytofkit' were integrated in this GUI, 
#' and each parameter has a help button to show the instruction. 
#' \code{cytofkit} analysis will be launched after submitting.
#' 
#' @author Hao Chen
#' @return the GUI for \code{cytofkit-package}
#' @import tcltk
#' @export
#' @seealso \code{\link{cytofkit-package}}, \code{\link{cytofkit}}
#' @references \url{http://signbioinfo.github.io/cytofkit/}
#' @examples
#' #cytofkit_GUI()  # remove the hash symbol to run
cytofkit_GUI <- function() {
    
    ##--------------------------##
    ## parameter initialization ##
    ##--------------------------##
    
    fcsFiles <- ""
    cur_dir <- getwd()
    mergeMethods <- c("all", "min", "ceil", "fixed")
    transformMethods <- c("autoLgcl", "cytofAsinh", "none")
    vizMethods <- c("pca", "isomap", "tsne", "NULL")
    clusterMethods <- c("Rphenograph", "ClusterX", "DensVM", "NULL")
    progressionMethods <- c("diffusionmap", "isomap", "NULL")
    
    rawFCSdir <- tclVar(cur_dir)
    fcsFile <- tclVar("")
    resDir <- tclVar(cur_dir)
    projectName <- tclVar("cytofkit")
    mergeMethod <- tclVar("ceil")
    fixedNum <- tclVar("5000")
    markers <- tclVar("")
    transformMethod <- tclVar("autoLgcl")
    progressionMethod <- tclVar("NULL")
    
    clusterSelect <- c()
    i <- 1
    while (i <= length(clusterMethods)) {
        aux <- paste("clusterSelect", i, sep = "")
        clusterSelect <- c(clusterSelect, as.character(aux))
        tclvalue(clusterSelect[i]) <- "0"
        i <- i + 1
    }
    tclvalue(clusterSelect[1]) <- "1"  # default Rphenograph
    
    vizSelect <- c()
    i <- 1
    while (i <= length(vizMethods)) {
        aux <- paste("vizSelect", i, sep = "")
        vizSelect <- c(vizSelect, as.character(aux))
        tclvalue(vizSelect[i]) <- "0"
        i <- i + 1
    }
    tclvalue(vizSelect[3]) <- "1"  # default tsne 
    
    ret_var <- tclVar("")
    
    ##-------------------##
    ##  button functions ##
    ##-------------------##
    
    reset_rawFCS_dir <- function() {
        rawFCS_dir <- ""
        rawFCS_dir <- tclvalue(tkchooseDirectory(title = "Choose your rawFCS dircetory ..."))
        if (rawFCS_dir != "") {
            tclvalue(rawFCSdir) <- rawFCS_dir
            tclvalue(resDir) <- rawFCS_dir
        }
    }
    
    reset_res_dir <- function() {
        res_dir <- ""
        res_dir <- tclvalue(tkchooseDirectory(title = "Choose your result dircetory ..."))
        if (res_dir != "") {
            tclvalue(resDir) <- res_dir
        }
    }
    
    reset_fcs_data <- function() {
        fnames <- ""
        fnames <- tk_choose.files(default = paste(tclvalue(rawFCSdir), 
            "fcs", sep = .Platform$file.sep), caption = "Select FCS files", 
            multi = TRUE, filters = matrix(c("{fcs files}", "{.fcs}"), 
                1, 2), index = 1)
        if (length(fnames) >= 1) {
            fnames <- fnames[!(grepl(paste0(.Platform$file.sep, 
                "fcs$"), fnames))]  # remove empty .fcs files
            tclvalue(fcsFile) <- paste(fnames, collapse = "}{")
        }
    }
    
    reset_para_data <- function() {
        
        if (tclvalue(fcsFile) == "" && tclvalue(rawFCSdir) == "") {
            tkmessageBox(title = "cytofkit: an integrated mass cytometry data analysis pipeline", 
                message = "Please input your \"rawFCSdirectory\" or \"fcsFile\".", 
                icon = "info", type = "ok")
        }
        if (tclvalue(fcsFile) != "") {
            fcsFiles <- strsplit(tclvalue(fcsFile), "}{", fixed = TRUE)[[1]]
        }
        fcsDir <- tclvalue(rawFCSdir)
        selectMarkers <- getParameters_GUI(fcsFiles, fcsDir)
        
        if (length(selectMarkers) > 0) {
            tclvalue(markers) <- paste(selectMarkers, collapse = "}{")
        }
    }
    
    reset_num2null <- function() {
        tclvalue(fixedNum) <- "NULL"
    }
    
    reset_num2any <- function() {
        tclvalue(fixedNum) <- "5000"
    }
    
    rawFCSdir_help <- function() {
        tkmessageBox(title = "rawFCSdir", message = "The directory that contains fcs files.", 
            icon = "info", type = "ok")
    }
    
    fcsFile_help <- function() {
        tkmessageBox(title = "fcsFiles", message = "The fcs files to be analyzed. One or multiple fcs files are allowed. When multiple fcs files are selected, cells from each fcs file are combined for analysis.", 
            icon = "info", type = "ok")
    }
    
    para_help <- function() {
        tkmessageBox(title = "markers", message = "Select the list of makers to be used for analysis.", 
            icon = "info", type = "ok")
    }
    
    projectName_help <- function() {
        tkmessageBox(title = "projectName", message = "A prefix that will be added to the names of result files.", 
            icon = "info", type = "ok")
    }
    
    mergeMethod_help <- function() {
        tkmessageBox(title = "mergeMethod", message = "When multiple fcs files are selected, cell expression data can be merged using one of the four different methods including \"ceil\",\"all\", \"min\",\"fixed\". \n\n\"ceil\" (the default option): up to a fixed number (specified by fixedNum) of cells are sampled without replacement from each fcs file and combined for analysis. \n\n\"all\": all cells from each fcs file are combined for analysis. \n\n\"min\": The minimum number of cells among all the selected fcs files are sampled from each fcs file and combined for analysis. \n\n\"fixed\": a fixed num (specified by fixedNum) of cells are sampled (with replacement when the total number of cell is less than fixedNum) from each fcs file and combined for analysis.", 
            icon = "info", type = "ok")
    }
    
    fixedNum_help <- function() {
        tkmessageBox(title = "fixedNum", message = "Up to fixedNum of cells from each fcs file are used for analysis.", 
            icon = "info", type = "ok")
    }
    
    
    transformMethod_help <- function() {
        tkmessageBox(title = "transformationMethod", message = "Data Transformation method, including \"cytofAsinh\" and \"autoLgcl\".", 
            icon = "info", type = "ok")
    }
    
    cluster_help <- function() {
        tkmessageBox(title = "clusterMethods", message = "The method(s) for clustering, including \"DensVM\", \"ClusterX\", \"Rphenograph\". \n\nIf \"NULL\" was selected, no clustering will be performed.", 
            icon = "info", type = "ok")
    }
    
    visualizationMethods_help <- function() {
        tkmessageBox(title = "visualizationMethods", message = "The method(s) used for visualizing the clustering results, multiple selections are allowed. Including \"pca\", \"isomap\", \"tsne\". \n\nWARNING: \"tsne\" is the default selection, \"isomap\" may take long time.", 
            icon = "info", type = "ok")
    }
    
    progressionMethod_help <- function() {
        tkmessageBox(title = "progressionMethod", message = "The method used for cellular progression analysis including \"diffusion map\" and \"isomap\"\n\nIf \"NULL\" was selected, no progression estimation will be performed.", 
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
        tclvalue(projectName) = "cytofkit"
        tclvalue(mergeMethod) = mergeMethods[3]
        tclvalue(fixedNum) = "5000"
        tclvalue(markers) = ""
        tclvalue(transformMethod) = "autoLgcl"
        tclvalue(clusterSelect[1]) = "0"
        tclvalue(clusterSelect[2]) = "1"
        tclvalue(clusterSelect[3]) = "0"
        tclvalue(clusterSelect[4]) = "0"
        tclvalue(vizSelect[1]) <- "0"
        tclvalue(vizSelect[2]) <- "0"
        tclvalue(vizSelect[3]) <- "1"
        tclvalue(progressionMethod) <- "NULL"
    }
    
    submit <- function() {
        has_error = FALSE
        if (tclvalue(markers) == "") {
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
        tkdestroy(tt)
    }
    
    
    ##----------------##
    ##  build the GUI ##
    ##--------------- ##
    
    ## head line
    tt <- tktoplevel(borderwidth = 20)
    tkwm.title(tt, "cytofkit: an integrated analysis pipeline for mass cytometry data")
    
    if(.Platform$OS.type == "windows"){
        box_length <- 63
    }else{
        box_length <- 55 
    }
    cell_width <- 3
    bt_width <- 8
    #hb_width <- 8
    
    imgfile <- system.file("extdata", "help.gif", package = "cytofkit")
    image1 <- tclVar()
    tkimage.create("photo", image1, file = imgfile)
    image2 <- tclVar()
    tkimage.create("photo", image2)
    tcl(image2, "copy", image1, subsample = 6)
    
    ## rawFCSdir
    rawFCSdir_label <- tklabel(tt, text = "Raw FCS Directory :")
    rawFCSdir_entry <- tkentry(tt, textvariable = rawFCSdir, width = box_length)
    rawFCSdir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, command = reset_rawFCS_dir)
    rawFCSdir_hBut <- tkbutton(tt, image = image2, command = rawFCSdir_help)
    
    ## resDir
    resDir_label <- tklabel(tt, text = "Result Directory :")
    resDir_entry <- tkentry(tt, textvariable = resDir, width = box_length)
    resDir_button <- tkbutton(tt, text = " Choose... ", width = bt_width, 
        command = reset_res_dir)
    resDir_hBut <- tkbutton(tt, image = image2, command = resDir_help)
    
    ## fcsFiles
    fcsFile_label <- tklabel(tt, text = "FCS File(s) :")
    fcsFile_entry <- tkentry(tt, textvariable = fcsFile, width = box_length)
    fcsFile_button <- tkbutton(tt, text = " Select... ", width = bt_width, 
        command = reset_fcs_data)
    fcsFile_hBut <- tkbutton(tt, image = image2, command = fcsFile_help)
    
    ## markers
    markers_label <- tklabel(tt, text = "Markers :")
    markers_entry <- tkentry(tt, textvariable = markers, width = box_length)
    markers_button <- tkbutton(tt, text = " Select... ", width = bt_width, 
        command = reset_para_data)
    markers_hBut <- tkbutton(tt, image = image2, command = para_help)
    
    ## projectName
    projectName_label <- tklabel(tt, text = "Project Name :")
    projectName_entry <- tkentry(tt, textvariable = projectName, width = box_length)
    projectName_hBut <- tkbutton(tt, image = image2, command = projectName_help)
    
    ## mergeMethod && fixedNum
    mergeMethod_label <- tklabel(tt, text = "Merge Method :")
    mergeMethod_hBut <- tkbutton(tt, image = image2, command = mergeMethod_help)
    merge_method_rbuts <- tkframe(tt)
    tkpack(tklabel(merge_method_rbuts, text = ""), side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethods[1], 
        variable = mergeMethod, value = mergeMethods[1], command = reset_num2null), 
        side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethods[2], 
        variable = mergeMethod, value = mergeMethods[2], command = reset_num2null), 
        side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethods[3], 
        variable = mergeMethod, value = mergeMethods[3], command = reset_num2any), 
        side = "left")
    tkpack(tkradiobutton(merge_method_rbuts, text = mergeMethods[4], 
        variable = mergeMethod, value = mergeMethods[4], command = reset_num2any), 
        side = "left")
    tkpack(tkentry(merge_method_rbuts, textvariable = fixedNum, 
                   width = 9), side = "right")
    tkpack(tklabel(merge_method_rbuts, text = "Fixed Number :"), 
        side = "right")
    tkpack(tklabel(merge_method_rbuts, text = "                 "), 
        side = "left")
    fixedNum_hBut <- tkbutton(tt, image = image2, command = fixedNum_help)
    
    ## transformMethod
    transformMethod_label <- tklabel(tt, text = "Transformation Method :")
    transformMethod_hBut <- tkbutton(tt, image = image2,
        command = transformMethod_help)
    transformMethod_rbuts <- tkframe(tt)
    tkpack(tklabel(transformMethod_rbuts, text = ""), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = transformMethods[1], 
        variable = transformMethod, value = transformMethods[1]), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = transformMethods[2],
        variable = transformMethod, value = transformMethods[2]), side = "left")
    tkpack(tkradiobutton(transformMethod_rbuts, text = transformMethods[3],
        variable = transformMethod, value = transformMethods[3]), side = "left")
    
    ## cluster method
    cluster_label <- tklabel(tt, text = "Cluster Method(s) :")
    cluster_hBut <- tkbutton(tt, image = image2, command = cluster_help)
    
    clusterMethods_cbuts <- tkframe(tt)
    tkpack(tklabel(clusterMethods_cbuts, text = ""), side = "left")
    i <- 1
    while (i <= length(clusterMethods)) {
        t <- tkcheckbutton(clusterMethods_cbuts, text = clusterMethods[i], 
                           variable = eval(clusterSelect[i]))
        tkpack(t, side = "left")
        i <- i + 1
    }
    
    ## visualizationMethods
    visualizationMethods_label <- tklabel(tt, text = "Visualization Method(s) :")
    visualizationMethods_hBut <- tkbutton(tt, image = image2,
        command = visualizationMethods_help)
    visualizationMethods_cbuts <- tkframe(tt)
    tkpack(tklabel(visualizationMethods_cbuts, text = ""), side = "left")
    i <- 1
    while (i <= length(vizMethods)) {
        t <- tkcheckbutton(visualizationMethods_cbuts, text = vizMethods[i], 
            variable = eval(vizSelect[i]))
        tkpack(t, side = "left")
        i <- i + 1
    }
    
    ## progressionMethod
    progressionMethod_label <- tklabel(tt, text = "Cellular Progression :")
    progressionMethod_hBut <- tkbutton(tt, image = image2, 
                                     command = progressionMethod_help)
    progressionMethod_rbuts <- tkframe(tt)
    tkpack(tklabel(progressionMethod_rbuts, text = ""), side = "left")
    tkpack(tkradiobutton(progressionMethod_rbuts, text = progressionMethods[1], 
                         variable = progressionMethod, value = progressionMethods[1]), side = "left")
    tkpack(tkradiobutton(progressionMethod_rbuts, text = progressionMethods[2], 
                         variable = progressionMethod, value = progressionMethods[2]), side = "left")
    tkpack(tkradiobutton(progressionMethod_rbuts, text = progressionMethods[3], 
                         variable = progressionMethod, value = progressionMethods[3]), side = "left")
    
    ## submit / reset / quit
    submit_button <- tkbutton(tt, text = "Submit", command = submit)
    reset_button <- tkbutton(tt, text = "Reset", command = reset)
    quit_button <- tkbutton(tt, text = "Quit", command = quit)
    
    ## display GUI
    tkgrid(rawFCSdir_label, rawFCSdir_hBut, rawFCSdir_entry, rawFCSdir_button, 
        padx = cell_width)
    tkgrid.configure(rawFCSdir_label, rawFCSdir_entry, rawFCSdir_button, 
        sticky = "e")
    tkgrid.configure(rawFCSdir_hBut, sticky = "e")
    
    tkgrid(fcsFile_label, fcsFile_hBut, fcsFile_entry, fcsFile_button, 
        padx = cell_width)
    tkgrid.configure(fcsFile_label, fcsFile_entry, fcsFile_button, 
        sticky = "e")
    tkgrid.configure(fcsFile_hBut, sticky = "e")
    
    tkgrid(markers_label, markers_hBut, markers_entry, markers_button, 
        padx = cell_width)
    tkgrid.configure(markers_label, markers_entry, markers_button, 
        sticky = "e")
    tkgrid.configure(markers_hBut, sticky = "e")
    
    tkgrid(resDir_label, resDir_hBut, resDir_entry, resDir_button, 
        padx = cell_width)
    tkgrid.configure(resDir_label, resDir_entry, resDir_button, 
        sticky = "e")
    tkgrid.configure(resDir_hBut, sticky = "e")
    
    tkgrid(projectName_label, projectName_hBut, projectName_entry, padx = cell_width)
    tkgrid.configure(projectName_label, projectName_entry, sticky = "e")
    tkgrid.configure(projectName_hBut, sticky = "e")
    
    tkgrid(mergeMethod_label, mergeMethod_hBut, merge_method_rbuts, 
        fixedNum_hBut, padx = cell_width)
    tkgrid.configure(mergeMethod_label, sticky = "e")
    tkgrid.configure(mergeMethod_hBut, sticky = "e")
    tkgrid.configure(merge_method_rbuts, sticky = "w")
    tkgrid.configure(fixedNum_hBut, sticky = "w")
    
    tkgrid(transformMethod_label, transformMethod_hBut, transformMethod_rbuts,
        padx = cell_width)
    tkgrid.configure(transformMethod_label, sticky = "e")
    tkgrid.configure(transformMethod_rbuts, sticky = "w")
    tkgrid.configure(transformMethod_hBut, sticky = "e")
    
    tkgrid(cluster_label, cluster_hBut, clusterMethods_cbuts, padx = cell_width)
    tkgrid.configure(cluster_label, sticky = "e")
    tkgrid.configure(clusterMethods_cbuts, sticky = "w")
    tkgrid.configure(cluster_hBut, sticky = "e")
    
    tkgrid(visualizationMethods_label, visualizationMethods_hBut, 
        visualizationMethods_cbuts, padx = cell_width)
    tkgrid.configure(visualizationMethods_label, sticky = "e")
    tkgrid.configure(visualizationMethods_cbuts, sticky = "w")
    tkgrid.configure(visualizationMethods_hBut, sticky = "e")
    
    tkgrid(progressionMethod_label, progressionMethod_hBut, progressionMethod_rbuts, 
           padx = cell_width)
    tkgrid.configure(progressionMethod_label, sticky = "e")
    tkgrid.configure(progressionMethod_rbuts, sticky = "w")
    tkgrid.configure(progressionMethod_hBut, sticky = "e")
    
    tkgrid(tklabel(tt, text = "\n"), padx = cell_width)  # leave blank line
    
    tkgrid(reset_button, tklabel(tt, text = ""), submit_button, 
        quit_button, padx = cell_width)
    tkgrid.configure(reset_button, sticky = "e")
    tkgrid.configure(quit_button, sticky = "w")
    
    tkwait.window(tt)
    
    ##-------------------##
    ## Return parameters ##
    ##-------------------##
    
    if (tclvalue(ret_var) != "OK") {
        okMessage <- "Analysis is cancelled."
    }else{
        fcsFiles <- strsplit(tclvalue(fcsFile), "}{", fixed = TRUE)[[1]]
        parameters <- strsplit(tclvalue(markers), "}{", fixed = TRUE)[[1]]
        
        clusterCheck <- c()
        i <- 1
        while (i <= length(clusterMethods)) {
            v <- as.numeric(tclvalue(clusterSelect[i])) > 0
            clusterCheck <- c(clusterCheck, v)
            i <- i + 1
        }
        
        vizCheck <- c()
        i <- 1
        while (i <= length(vizMethods)) {
            v <- as.numeric(tclvalue(vizSelect[i])) > 0
            vizCheck <- c(vizCheck, v)
            i <- i + 1
        }
        
        inputs <- list()
        inputs[["fcsFiles"]] <- fcsFiles
        inputs[["markers"]] <- parameters
        inputs[["mergeMethod"]] <- tclvalue(mergeMethod)
        inputs[["fixedNum"]] <- suppressWarnings(as.numeric(tclvalue(fixedNum)))
        inputs[["transformMethod"]] <- tclvalue(transformMethod)
        inputs[["dimReductionMethod"]] <- "tsne"
        inputs[["clusterMethods"]] <- clusterMethods[clusterCheck]
        inputs[["visualizationMethods"]] <- vizMethods[vizCheck]
        inputs[["progressionMethod"]] <- tclvalue(progressionMethod)
        inputs[["projectName"]] <- tclvalue(projectName)
        inputs[["resultDir"]] <- tclvalue(resDir)
        
        
        ## pass the parameters and run the cytofkit function
        cytofkit(fcsFiles = inputs[["fcsFiles"]], 
                 markers = inputs[["markers"]],
                 projectName = inputs[["projectName"]], 
                 mergeMethod = inputs[["mergeMethod"]], 
                 fixedNum = inputs[["fixedNum"]], 
                 transformMethod = inputs[["transformMethod"]], 
                 dimReductionMethod = inputs[["dimReductionMethod"]], 
                 clusterMethods = inputs[["clusterMethods"]], 
                 visualizationMethods = inputs[["visualizationMethods"]], 
                 progressionMethod = inputs[["progressionMethod"]], 
                 clusterSampleSize = 500,
                 resultDir = inputs[["resultDir"]], 
                 saveResults = TRUE, 
                 saveObject = TRUE)
        
        okMessage <- paste0("Analysis Done, results are saved under ",
                            inputs[["resultDir"]])
    }
    
    launchShinyAPP_GUI(message = okMessage, dir = inputs[["resultDir"]])
}


#' GUI for launching shiny APP
#' 
#' A shiny APP for interactive exploration of the analysis results
#' 
#' @param message A message to determine if open the shiny APP
#' @param dir Result direcroty.
#' 
#' @export
#' @examples
#' # launchShinyAPP_GUI()
launchShinyAPP_GUI <- function(message="cytofkit", dir = getwd()){
    ifAPP <- tclVar("n")
    ss <- tktoplevel(borderwidth = 10)
    tkwm.title(ss, "cytofkit: Analysis Done")
    
    onYes <- function() {
        tclvalue(ifAPP) <- "y"
        tkdestroy(ss)
    }
    
    onNo <- function() {
        tclvalue(ifAPP) <- "n"
        tkdestroy(ss)
    }
    yesBut <- tkbutton(ss, text = " YES ", command = onYes)
    noBut <- tkbutton(ss, text = " NO ", command = onNo)
    openDirBut <- tkbutton(ss, text = "Open", command = function(){opendir(dir)})
    okBut <- tkbutton(ss, text = "OK", command = function(){tkdestroy(ss)})
    tkgrid(tklabel(ss, text = message))
    
    if(message != "Analysis is cancelled."){
        tkgrid(openDirBut)
        tkgrid(tklabel(ss, text = "\n"))
        tkgrid(tklabel(ss, text = "Launch Shiny APP to check your reuslts:"))
        tkgrid(noBut, tklabel(ss, text = "    "), yesBut)
        tkgrid.configure(noBut, sticky = "e")
        tkgrid.configure(yesBut, sticky = "e")
    }else{
        tkgrid(tklabel(ss, text = "\n"))
        tkgrid(okBut)
    }
    
    tkwait.window(ss)
    
    if(tclvalue(ifAPP) == "y"){
        cytofkitShinyAPP()
    }
}

## function for opening the results directory
opendir <- function(dir = getwd()){
    if (.Platform['OS.type'] == "windows"){
        shell.exec(dir)
    } else {
        system(paste(Sys.getenv("R_BROWSER"), dir))
    }
}


#' GUI for marker selection 
#' 
#' Extract the markers from the fcsfiles
#' 
#' @param fcsFile The name of the FCS file
#' @param rawFCSdir The path of the FCS file
#' @examples 
#' #getParameters_GUI()
getParameters_GUI <- function(fcsFile, rawFCSdir) {
    
    if (missing(fcsFile)) {
        fcsFile <- list.files(path = rawFCSdir, pattern = ".fcs$", full.names = TRUE)
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
    tl <- tklistbox(mm, height = 30, width = 40, selectmode = "multiple", yscrollcommand = function(...) tkset(scr, 
        ...), background = "white")
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
