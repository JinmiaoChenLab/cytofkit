## max data size
options(shiny.maxRequestSize=1024^10) 

shinyServer(function(input, output, session) {
    
    ##------------------Reactive Values and Reactive Objects-------------------
    
    v <- reactiveValues(data = NULL, sampleInfo = NULL)
    
    ## Scatter plot methods
    visualizationMethods <- reactive({
        if(is.null(v$data)){
            return(NULL)
        }else{
            return(v$data$visualizationMethods)
        }
    })
    
    ## Scatter plot functions
    visualizationFunctions <- reactive({
        if(is.null(v$data)){
            return(NULL)
        }else{
            return(c(names(v$data$clusterRes), 
                     "ColorBySample",
                     "DensityPlot",
                     "DotPlot"))
        }
    })
    
    ## cluster methods
    clusterMethods <- reactive({
        if(is.null(v$data))
            return(NULL)
        cMethods <- names(v$data$clusterRes)
        return(cMethods)
    })
    
    ## progression labs
    progressionLabs <- reactive({
        if(is.null(v$data))
            return(NULL)
        if(is.null(v$data$progressionRes))
            return(NULL)
        progressionLabs <- colnames(v$data$progressionRes[[3]])
        return(progressionLabs)
    })
    
    
    ##--------------------------------Side Panel-------------------------------
    
    ## Load cytofkit RData object
    observeEvent(input$goButton, {
        cytofkitObj <- input$cytofkitObj
        if (is.null(cytofkitObj)){
            v$data <- NULL
        }else{
            cat(cytofkitObj$datapath)
            load(cytofkitObj$datapath)
            v$data <- analysis_results
            
            if(is.null(v$data$projectName)){
                v$data$projectName <- "cytofkit_shinyAPP_output"
            }
            
            # Need modification later
            # currently doesn't update sampleInfo with v$data$sampleInfo
            v$sampleInfo <- data.frame(cellID = row.names(v$data$expressionData),
                                       cellSample = factor(sub("_[0-9]*$", "", row.names(v$data$expressionData))),
                                       stringsAsFactors = FALSE)
            v$data$sampleInfo <- v$sampleInfo
        }
    })

    output$sampleSelect <- renderUI({
        if(is.null(v$data) || is.null(v$sampleInfo)){
            return(NULL)
        }else{
            sampleNames <- unique(as.character(v$sampleInfo$cellSample))
            checkboxGroupInput('samples', NULL, 
                               sampleNames, selected = sampleNames)
        }   
    })
    
    output$summaryText1 <- renderText({
        if(is.null(v$data))
            return(NULL)
        paste0("-- ", nrow(v$data[[1]]), "cells x ", ncol(v$data[[1]]), "markers")
    })
    
    output$summaryText2 <- renderText({
        if(is.null(v$data))
            return(NULL)
        paste0("-- ", paste(names(v$data$clusterRes), collapse = " | "))
    })
    
    output$summaryText3 <- renderText({
        if(is.null(v$data))
            return(NULL)
        paste0("-- ", paste(v$data$visualizationMethods, collapse =  " | "))
    })
    
    output$summaryText4 <- renderText({
        if(is.null(v$data))
            return(NULL)
        paste0("-- ", ifelse(is.null(v$data$progressionRes), "NULL", 
                             sub("_[0-9]*$", "", colnames(v$data$progressionRes$progressionData)[1])))
    })
    
    ## Save and parse cytofkit RData object
    observeEvent(input$saveButton, {
        if (!is.null(v$data)){
            withProgress(message='Saving Results ', value=0, {
                ## check results saving path
                if(is.null(v$data$resultDir) || !dir.exists(v$data$resultDir)){
                    v$data$resultDir <- path.expand("~")  ## default save to home if not specified
                }
                saveToFCS <- TRUE
                if(is.null(v$data$rawFCSdir)){
                    saveToFCS <- FALSE
                    warning("Path for original FCS files is not provided, 
                            data cannnot be saved to new copies of FCS files.")
                }else if(!dir.exists(v$data$rawFCSdir)){
                    saveToFCS <- FALSE
                    warning(paste0("Path for original FCS files doesn't exits, 
                                   data cannnot be saved to new copies of FCS files.", 
                                   "Please check path: ", v$data$rawFCSdir))
                }
                incProgress(1/2, message = paste0("To ", v$data$resultDir))
                analysis_results <<- v$data
                cytof_writeResults(analysis_results,
                                   saveToFCS = saveToFCS)
                incProgress(1/2)
                ## open the results directory
                opendir(v$data$resultDir)
            })
        }
    })
    
    ##------------------------------Cluster Panel------------------------------
    
    ## cluster plot
    
    output$C_PlotMethod <- renderUI({
        if(is.null(v$data) || is.null(visualizationMethods())){
            return(NULL)
        }else{
            selectInput('c_PlotMethod', 'Plot Data:', choices = visualizationMethods(), 
                        selected = visualizationMethods()[1], width = "100%")
        }   
    })
    
    output$C_PlotFunction <- renderUI({
        if(is.null(v$data) || is.null(visualizationFunctions())){
            return(NULL)
        }else{
            selectInput('c_PlotFunction', 'Plot Option:', choices = visualizationFunctions(), 
                        selected = visualizationFunctions()[1], width = "100%")
        }   
    })
    
    output$C_markerSelect <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            markerNames <- colnames(v$data$expressionData)
            checkboxGroupInput('c_markerSelect', strong('Select Markers:'),
                               markerNames, selected = markerNames, inline = TRUE)
        }   
    })
    
    output$C_ScatterPlot <- renderPlot({
        if(is.null(v$data) || is.null(input$c_PlotMethod) || is.null(input$c_PlotFunction)){
            return(NULL)
        }else{
            withProgress(message="Generating Cluster Scatter Plot", value=0, {
                gp <- scatterPlot(obj = v$data,
                                  plotMethod = input$c_PlotMethod,
                                  plotFunction = input$c_PlotFunction,
                                  pointSize = input$S_PointSize,
                                  addLabel = input$addLabel,
                                  labelSize = input$S_LabelSize,
                                  sampleLabel = FALSE,
                                  FlowSOM_k = input$S_FlowSOM_k, 
                                  selectSamples = input$samples, 
                                  facetPlot = input$facetPlot,
                                  colorPalette = input$colorPalette,
                                  labelRepel = input$labelRepel,
                                  removeOutlier = TRUE)
                incProgress(1/2)
                plot(gp)
                incProgress(1/2)
            })
        }
    }, height = 700, width = 750)
    
    ## annotate clusters
    
    output$C_labelCluster <- renderUI({
        if(is.null(v$data) || is.null(v$data$clusterRes)){
            return(NULL)
        }else{
            clusterMethods <- c(names(v$data$clusterRes)) 
            #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
            selectInput('c_labelCluster', 'Choose Cluster Results to Annotate:', 
                        choices = clusterMethods, 
                        selected = clusterMethods[1], width = "30%")
        }   
    })
    
    ## currently use 100 as a limit for cluster numbers 
    ## --- TODO: use reactiveValues to automatically retrive cluster numbers --- ## 
    lapply(1:100, function(i) {
        output[[paste0('Cluster', i)]] <- renderUI({
            if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_labelCluster)){
                return(NULL)
            }
            
            clusters <- sort(unique(v$data$clusterRes[[input$c_labelCluster]]))
            if (i <= length(clusters)){
                x <- clusters[i]
                textInput(paste0('cluster', i), paste0('Cluster ', x," :"), 
                          value = "", width = "30%", placeholder = "Type in the cell type")
            }
        })
    })
    
    ## update cluster labels
    observeEvent(input$updatelabel, {
        if(!is.null(v$data) && !is.null(input$c_labelCluster)){
            obj <- v$data
            clusterMethod <- input$c_labelCluster
            clusterVec<- obj$clusterRes[[clusterMethod]]
            clusterLabels <- clusterVec
            clusters <- sort(unique(clusterVec))
            
            for (i in 1:length(clusters)){
                clusteri <- clusters[i]
                ilabel <- input[[paste0('cluster', i)]]
                if(ilabel == ""){
                    clusterLabels[clusterLabels==clusteri] <- "Unknown"
                }else{
                    clusterLabels[clusterLabels==clusteri] <- ilabel
                }
            }
            ## either add new cluste or update
            labelName <- clusterMethod
            labelName <- ifelse(grepl("Annotated_", labelName),
                                clusterMethod,
                                paste0("Annotated_", clusterMethod))
            obj$clusterRes[[labelName]] <- clusterLabels
            
            ## update the project name
            obj$projectName <- paste0(obj$projectName, "_annotated_", clusterMethod)
            
            v$data <- obj
            
            ## jump to C_panel1
            updateTabsetPanel(session, "C_clusterTabs", selected = "C_panel1")
        }
    })
    
    ## result object which will be updated by C_runFlowSOM
    observeEvent(input$C_runFlowSOM, {
        if(!is.null(v$data) && !is.null(input$c_markerSelect)){
            obj <- v$data
            withProgress(message=paste0('Runing FlowSOM using k=', input$S_FlowSOM_k), value=0, {
                FlowSOM_cluster <- cytof_cluster(xdata = obj$expressionData[ ,input$c_markerSelect],
                                                 method = "FlowSOM",
                                                 FlowSOM_k = input$S_FlowSOM_k)
                incProgress(1/2)
                ## update FlowSOM cluster results
                obj$clusterRes[["FlowSOM"]] <- FlowSOM_cluster
                ## update the project name
                obj$projectName <- paste0(obj$projectName, "_add_FlowSOM")
                v$data <- obj
                incProgress(1/2)
            })
            
            ## jump to C_panel1
            updateTabsetPanel(session, "C_clusterTabs", selected = "C_panel1")
        }
    })
    
    ##------------------------------Marker Panel-------------------------------
    
    ## level plot
    
    output$M_PlotMethod <- renderUI({
        if(is.null(v$data) || is.null(visualizationMethods())){
            return(NULL)
        }else{
            selectInput('m_PlotMethod', 'Plot Method:', choices = visualizationMethods(), 
                        selected = visualizationMethods()[1], width = "100%")
        }   
    })
    
    output$M_PlotMarker <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            markers <- c(colnames(v$data$expressionData), "All Markers", "All Markers(scaled)")
            selectInput('m_PlotMarker', 'Plot Marker:', choices = markers, 
                        selected = markers[1], width = "100%")
        }   
    })
    
    output$M_markerExpressionPlot <- renderPlot({
        if(is.null(v$data) || is.null(input$m_PlotMethod) || is.null(input$m_PlotMarker)){
            return(NULL)
        }else{
            withProgress(message="Generating Marker Expression Plot", value=0, {
                gp <- scatterPlot(obj = v$data,
                                  plotMethod = input$m_PlotMethod,
                                  plotFunction = input$m_PlotMarker,
                                  pointSize = input$M_PointSize,
                                  addLabel = input$addLabel,
                                  labelSize = input$S_LabelSize,
                                  sampleLabel = FALSE,
                                  FlowSOM_k = input$S_FlowSOM_k, 
                                  selectSamples = input$samples, 
                                  facetPlot = input$facetPlot,
                                  colorPalette = input$colorPalette,
                                  labelRepel = input$labelRepel,
                                  removeOutlier = TRUE)
                incProgress(1/2)
                plot(gp)
                incProgress(1/2)
            })
        }
        
    }, height = 800, width = 850)
    
    ## histogram plot
    
    output$M_stackFactor <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            stackFactorChoice <- c(names(v$data$clusterRes), "sample") 
            selectInput('m_stackFactor', 'Stack Factor:', choices = stackFactorChoice, 
                        selected = stackFactorChoice[1], width = "100%")
        }   
    })
    
    output$M_markerSelect <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            markerNames <- colnames(v$data$expressionData)
            initNum <- ifelse(length(markerNames) >=4, 4, 1)
            selectizeInput('m_markerSelect', 'Select Markers:', 
                           choices = markerNames, selected = markerNames[1:initNum], 
                           multiple = TRUE, width = "100%")
            
            # checkboxGroupInput('m_markerSelect', strong('Select Markers:'),
            #                    markerNames, selected = markerNames[initNum], inline = TRUE)
        }   
    })
    
    observeEvent(input$M_updateDensityPlot, {
        m_markerSelect <- isolate(input$m_markerSelect)
        output$M_stackDensityPlot <- renderPlot({
            if(is.null(v$data) || is.null(input$m_stackFactor) || is.null(m_markerSelect)){
                return(NULL)
            }else{
                withProgress(message="Generating Stack Density Plot", value=0, {
                    data <- data.frame(v$data$expressionData, check.names = FALSE)
                    samples <- as.character(v$sampleInfo$cellSample)
                    mySamples <- samples %in% input$samples
                    sfactors <- data.frame(do.call(cbind, v$data$clusterRes), 
                                           sample = samples, 
                                           stringsAsFactors = FALSE, 
                                           check.names = FALSE)
                    data <- data[mySamples, ,drop=FALSE]
                    stackFactor <- sfactors[mySamples, input$m_stackFactor]
                    incProgress(1/3)
                    gp <- stackDenistyPlot(data = data, 
                                           densityCols=m_markerSelect, 
                                           stackFactor = stackFactor,
                                           kernel = "gaussian",
                                           bw = "nrd0", 
                                           adjust = 1,
                                           stackRotation = 0, 
                                           stackSeperation = "auto",
                                           x_text_size = input$M_xlab_size, 
                                           strip_text_size = input$M_markerTextSize,
                                           legend_text_size = input$M_legendTextSize, 
                                           legendRow = input$M_legendRow,
                                           legend_title = input$m_stackFactor)
                    incProgress(1/3)
                    plot(gp)
                    incProgress(1/3)
                })
            }
        }, height = 800, width = 850)
    })
    
    ## heat map plot
    
    output$H_plotCluster <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            selectInput('h_plotCluster', 'Cluster Method:', choices = clusterMethods(), 
                        selected = clusterMethods()[1], width = "100%")
        }   
    })
    
    output$H_heatmapPlot <- renderPlot({
        if(is.null(v$data) || is.null(input$h_plotCluster))
            return(NULL)
        heatMap(data = v$data, 
                clusterMethod = input$h_plotCluster, 
                type = input$H_plotMethod, 
                selectSamples = input$samples,
                cex_row_label= input$H_rowLabelSize, 
                cex_col_label= input$H_colLabelSize, 
                scaleMethod = input$H_scaleMethod)
    }, height = 800, width = 850)
    
    
    ##------------------------------Sample Panel-------------------------------
    
    ## cell counts
    
    output$S_clusterMethod <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            selectInput('s_clusterMethod', 'Cluster Method:', choices = clusterMethods(), 
                        selected = clusterMethods()[1], width = "100%")
        }   
    })
    
    output$S_clusterFilter <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod)){
            return(NULL)
        }else{
            clusterIDs <- sort(unique(v$data$clusterRes[[input$s_clusterMethod]]))
            selectizeInput('s_clusterFilter', 'Filter Clusters:', 
                           choices = clusterIDs, selected = clusterIDs, 
                           multiple = TRUE, width = "100%")
        }   
    })
    
    output$S_heatmapPlot <- renderPlot({
        if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod))
            return(NULL)
        heatMap(data = v$data, 
                clusterMethod = input$s_clusterMethod, 
                type = "percentage", 
                selectSamples = input$samples,
                cex_row_label= input$H_rowLabelSize, 
                cex_col_label= input$H_colLabelSize, 
                scaleMethod = input$H_scaleMethod)
    }, height = 800, width = 850)
    
    output$S_clusterTable <- renderTable({
        if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod)){
            return(NULL)
        }else{
            data <- data.frame(sample = v$sampleInfo$cellSample,
                               cluster = as.factor(v$data$clusterRes[[input$s_clusterMethod]]),
                               counts = 1)
            
            statData1 <- aggregate(counts ~ ., data = data, sum)
            statData2 <- aggregate(counts ~ sample, data = data, sum)
            statData <- merge(statData1, statData2, by="sample", suffixes = c("InAll","InSample"))
            if(is.numeric(statData$cluster)) statData$cluster <- as.integer(statData$cluster)
            statData$counts <- as.integer(statData$countsInAll)
            statData$percentageInAll <- round(statData$countsInAll/nrow(data), 4)
            statData$percentageInSample <- round(statData$countsInAll/statData$countsInSample, 2)
            statData[, c("sample", "cluster", "counts", "percentageInSample", "percentageInAll")]
        }   
    }) 
    
    output$S_rateChangePlot <- renderPlot({
        if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod) || is.null(input$s_clusterFilter))
            return(NULL)
        withProgress(message="Generating Rate Change Plot", value=0, {
            ## percentage stat
            data <- data.frame(sample = v$sampleInfo$cellSample,
                               cluster = as.factor(v$data$clusterRes[[input$s_clusterMethod]]),
                               counts = 1)
            statData1 <- aggregate(counts ~ ., data = data, sum)
            statData2 <- aggregate(counts ~ sample, data = data, sum)
            statData <- merge(statData1, statData2, by="sample", suffixes = c("InAll","InSample"))
            statData$percentageInSample <- statData$countsInAll/statData$countsInSample
            incProgress(1/3)
            ## filter clusters
            usedClusters <- input$s_clusterFilter
            clusterCheck <- as.character(statData$cluster) %in% usedClusters
            statData <- statData[clusterCheck, ,drop=FALSE]
            incProgress(1/3)
            gp <- ggplot(data = statData, aes_string(x="sample", 
                                               y="percentageInSample", 
                                               color = "cluster",
                                               group = "cluster")) + 
                geom_point(size = 2) + geom_line(size = 1.5) + 
                xlab("Cell Group") + ylab("Percentage of Cells in Group") + theme_bw() + 
                theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))
            incProgress(1/3)
            plot(gp)
        })
    }, height = 800, width = 850)
    
    
    ## group samples
    
    output$S_groupSamples <- renderUI({
        if(is.null(v$data) || is.null(v$data$clusterRes)){
            return(NULL)
        }else{
            clusterMethods <- c(names(v$data$clusterRes)) 
            #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
            selectInput('c_labelCluster', 'Choose Cluster Results to Annotate:', 
                        choices = clusterMethods, 
                        selected = clusterMethods[1], width = "30%")
        }   
    })
    
    ## currently use 100 as a limit for sample numbers 
    ## --- TODO: use reactiveValues to automatically retrive sample numbers --- ## 
    lapply(1:100, function(i) {
        output[[paste0('S_sample', i)]] <- renderUI({
            if(is.null(v$data) || is.null(v$sampleInfo)){
                return(NULL)
            }
            
            uniqueSampleNames <- sort(unique(v$sampleInfo$cellSample))
            if (i <= length(uniqueSampleNames)){
                x <- uniqueSampleNames[i]
                textInput(paste0('Sample', i), paste0(x," :"), 
                          value = "", width = "40%", 
                          placeholder = "Type in the group name for this sample")
            }
        })
    })
    
    ## update sample groups
    observeEvent(input$updateSampleGroups, {
        if(!is.null(v$data) && !is.null(v$sampleInfo)){
            v$sampleInfo$originalCellSample <- v$sampleInfo$cellSample
            uniqueSampleNames <- sort(unique(v$sampleInfo$originalCellSample))
            
            sampleGroupNames <- NULL
            for(i in 1:length(uniqueSampleNames)){
                sampleGroupNames <- c(sampleGroupNames, input[[paste0("Sample", i)]])
            }
            
            groupNameLevels <- strsplit(input$sampleGroupLevels, ";", fixed = TRUE)[[1]]
            
            if(groupNameLevels != "" && all(sampleGroupNames != "") 
               && length(groupNameLevels) == length(unique(sampleGroupNames))
               && all(as.character(groupNameLevels) %in% sampleGroupNames)){
                sampleMatchID <- match(v$sampleInfo$originalCellSample, uniqueSampleNames)
                v$sampleInfo$cellSample <- factor(sampleGroupNames[sampleMatchID],
                                                  levels = groupNameLevels)
            }else{
                sampleGroupNames[sampleGroupNames == ""] <- uniqueSampleNames[sampleGroupNames == ""]
                sampleMatchID <- match(v$sampleInfo$originalCellSample, uniqueSampleNames)
                v$sampleInfo$cellSample <- factor(sampleGroupNames[sampleMatchID])
            }
            
            cellID_number <- do.call(c, regmatches(v$sampleInfo$cellID, 
                                                   gregexpr("_[0-9]*$", v$sampleInfo$cellID, perl=T)))
            
            ## newCellID = "sampleGroup" + "_cellID" + "globalID" to avoid dumplicates
            v$sampleInfo$newCellID <- paste0(as.character(v$sampleInfo$cellSample), 
                                             cellID_number,
                                             1:length(cellID_number))
            
            
            ## update reactive object v$data
            expressionData <- v$data$expressionData
            row.names(expressionData) <- v$sampleInfo$newCellID
            v$data$expressionData <- expressionData
            
            ## update the project name
            v$data$projectName <- paste0(v$data$projectName, "_grouped_samples")

            ## update reactive object v$sampleInfo
            if(!is.null(v$data$progressionRes)){
                sampleExpressData <- v$data$progressionRes$sampleData
                row.names(sampleExpressData) <- v$sampleInfo$newCellID[match(row.names(sampleExpressData),
                                                                             v$sampleInfo$cellID)]
                v$data$progressionRes$sampleData <- sampleExpressData
            }
            
            ## jump to S_panel1
            updateTabsetPanel(session, "S_sampleTabs", selected = "S_panel1")
        }
    })
    
    ## revert old sample names
    observeEvent(input$revertSampleNames, {
        if(!is.null(v$data) && !is.null(v$sampleInfo)){
            if(!is.null(v$sampleInfo$originalCellSample)){
                v$sampleInfo$cellSample <- v$sampleInfo$originalCellSample
                v$sampleInfo$originalCellSample <- NULL
                
                ## update reactive object v$data
                expressionData <- v$data$expressionData
                row.names(expressionData) <- v$sampleInfo$cellID
                v$data$expressionData <- expressionData
                
                ## update the project name
                v$data$projectName <- sub("_grouped_samples", "", v$data$projectName)
                
                ## update reactive object v$sampleInfo
                if(!is.null(v$data$progressionRes)){
                    sampleExpressData <- v$data$progressionRes$sampleData
                    row.names(sampleExpressData) <- v$sampleInfo$cellID[match(row.names(sampleExpressData),
                                                                              v$sampleInfo$newCellID)]
                    v$data$progressionRes$sampleData <- sampleExpressData
                }
            }
            ## jump to S_panel1
            updateTabsetPanel(session, "S_sampleTabs", selected = "S_panel1")
        }
    })
    
    
    
    ##---------------------------Progression Panel------------------------------
    
    ## subset relationship plot
    
    output$P_xlab <- renderUI({
        if(is.null(v$data) || is.null(progressionLabs())){
            return(NULL)
        }else{
            selectInput('p_xlab', 'Plot X:', choices = progressionLabs(), 
                        selected = progressionLabs()[1], width = "100%")
        }   
    })
    
    output$P_ylab <- renderUI({
        if(is.null(v$data) || is.null(progressionLabs())){
            return(NULL)
        }else{
            selectInput('p_ylab', 'Plot Y:', choices = progressionLabs(), 
                        selected = progressionLabs()[2], width = "100%")
        }   
    })
    
    output$P_scatterPlot <- renderPlot({
        if(is.null(v$data) || is.null(v$data$progressionRes) || is.null(input$p_xlab) || is.null(input$p_ylab)){
            return(NULL)
        }else{
            withProgress(message="Generating Progression Scatter Plot", value=0, {
                obj <- v$data$progressionRes
                data <- data.frame(obj$progressionData, 
                                   cluster = obj$sampleCluster,
                                   sample = sub("_[0-9]*$", "", row.names(obj$sampleData)))
                incProgress(1/3)
                data <- data[data$sample %in% input$samples, ,drop=FALSE]
                gp <- cytof_clusterPlot(data = data, 
                                        xlab = input$p_xlab, 
                                        ylab = input$p_ylab, 
                                        cluster = "cluster", 
                                        sample = "sample",
                                        title = "Subset Relationship", 
                                        type = ifelse(input$facetPlot, 2, 1),
                                        point_size = input$P_PointSize, 
                                        addLabel = input$addLabel, 
                                        labelSize = input$P_LabelSize, 
                                        sampleLabel = FALSE, 
                                        labelRepel = input$labelRepel,
                                        fixCoord = FALSE)
                incProgress(1/3)
                plot(gp)
                incProgress(1/3)
            })
        }
    }, height = 700, width = 750)
    
    ## marker expression profile
    
    output$P_orderBy <- renderUI({
        if(is.null(v$data) || is.null(progressionLabs())){
            return(NULL)
        }else{
            selectInput('p_orderBy', 'Cell Order By:', choices = progressionLabs(), 
                        selected = progressionLabs()[1], width = "100%")
        }   
    })

    output$P_markerSelect <- renderUI({
        if(is.null(v$data) || is.null(v$data$progressionRes)){
            return(NULL)
        }else{
            markerNames <- colnames(v$data$progressionRes$sampleData)
            initNum <- ifelse(length(markerNames) >=4, 4, 1)
            selectizeInput('p_markerSelect', 'Select Markers:', 
                        choices = markerNames, selected = markerNames[1:initNum], 
                        multiple = TRUE, width = "100%")
            # checkboxGroupInput('p_markerSelect', strong('Select Markers:'), 
            #                    markerNames, selected = markerNames, inline = TRUE)
        }   
    })
    
    output$P_clusterSelect <- renderUI({
        if(is.null(v$data) || is.null(v$data$progressionRes)){
            return(NULL)
        }else{
            clusterIDs <- sort(unique(v$data$progressionRes$sampleCluster))
            selectizeInput('p_clusterSelect', 'Select Clusters:', 
                        choices = clusterIDs, selected = clusterIDs, 
                        multiple = TRUE, width = "100%")
            # checkboxGroupInput('p_clusterSelect', strong('Select Clusters:'), 
            #                    clusterIDs, selected = clusterIDs, inline = TRUE)
        }   
    })
    
    observeEvent(input$P_updateRegressionPlot, {
        p_markerSelect <- isolate(input$p_markerSelect)
        p_clusterSelect <- isolate(input$p_clusterSelect)
        output$P_markerPlot <- renderPlot({
            if(is.null(v$data) || is.null(v$data$progressionRes) || is.null(p_markerSelect) || is.null(p_clusterSelect) || is.null(input$p_orderBy))
                return(NULL)
            
            withProgress(message="Generating Marker Expression Profile", value=0, {
                data <- data.frame(v$data$progressionRes$sampleData,
                                   cluster = v$data$progressionRes$sampleCluster, 
                                   v$data$progressionRes$progressionData,
                                   check.names = FALSE)
                
                sampleNames <- sub("_[0-9]*$", "", row.names(v$data$progressionRes$sampleData))
                data <- data[sampleNames %in% input$samples, ,drop=FALSE]
                incProgress(1/3)
                if(input$P_combineTrends){
                    pp <- cytof_expressionTrends(data, 
                                                 markers = p_markerSelect, 
                                                 clusters = p_clusterSelect, 
                                                 orderCol = input$p_orderBy, 
                                                 clusterCol = "cluster", 
                                                 reverseOrder = input$P_reverseOrder,
                                                 addClusterLabel = input$addLabel,
                                                 clusterLabelSize = input$P_LabelSize2,
                                                 segmentSize = 0.5,
                                                 min_expr = NULL) 
                }else{
                    pp <- cytof_progressionPlot(data, 
                                                markers = p_markerSelect, 
                                                clusters = p_clusterSelect, 
                                                orderCol = input$p_orderBy, 
                                                clusterCol = "cluster", 
                                                reverseOrder = input$P_reverseOrder,
                                                addClusterLabel = input$addLabel,
                                                clusterLabelSize = input$P_LabelSize2,
                                                segmentSize = 0.5,
                                                min_expr = NULL) 
                }
                incProgress(1/3)
                plot(pp)
                incProgress(1/3)
            })
            
        }, height = 800, width = 850)  
    })
    
    ## Run Diffusionmap
    
    output$P_clusterTable <- renderTable({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            clusterTable <- t(as.matrix(table(v$data$clusterRes[[input$p_clusterMethod]])))
            out <- as.data.frame(clusterTable, row.names = "Cell Counts")
            colnames(out) <- paste("Cluster", colnames(out))
            out
        }   
    })
    
    output$P_clusterFilter <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            obj <- v$data
            clusterIDs <- sort(unique(obj$clusterRes[[input$p_clusterMethod]]))
            selectizeInput('p_clusterFilter', 'Filter Clusters:', 
                           choices = clusterIDs, selected = clusterIDs, 
                           multiple = TRUE, width = "100%")
        }   
    })
    
    output$P_clusterMethod <- renderUI({
        if(is.null(v$data) || is.null(clusterMethods())){
            return(NULL)
        }else{
            selectInput('p_clusterMethod', 'Cluster Method:', choices = clusterMethods(), 
                        selected = clusterMethods()[1], width = "100%")
        }   
    })
    
    ## result object which will be updated by P_runDiffusionmap
    observeEvent(input$P_runDiffusionmap, {
        
        if(!is.null(v$data)){
            obj <- v$data
            usedClusters <- input$p_clusterFilter
            clusterCheck <- obj$clusterRes[[input$p_clusterMethod]] %in% usedClusters
            mdata <- obj$expressionData[clusterCheck, ]
            mcluster <- obj$clusterRes[[input$p_clusterMethod]][clusterCheck]
            withProgress(message="Runing Diffusionmap", value=0, {
                diffmapRes <- cytof_progression(data = mdata, 
                                                cluster = mcluster, 
                                                method = "diffusionmap", 
                                                distMethod = input$P_distMethod,
                                                out_dim = input$P_outDim,
                                                clusterSampleMethod = input$P_sampleMethod,
                                                clusterSampleSize = input$P_clusterSampleSize)
                incProgress(1/2)
                ## update progressionRes results
                obj$progressionRes <- diffmapRes
                
                ## update the project name
                obj$projectName <- paste0(obj$projectName, "_added_diffusionmap")
                
                v$data <- obj
                incProgress(1/2)
            })
            
            ## jump to P_panel1
            updateTabsetPanel(session, "P_progressionTabs", selected = "P_panel1")
        }
    })
})



