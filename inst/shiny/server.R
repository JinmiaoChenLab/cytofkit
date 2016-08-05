## max data size
options(shiny.maxRequestSize=1024^10) 

shinyServer(function(input, output, session) {
   
    ##-----------------------------------Load Data--------------------------------------
    v <- reactiveValues(data = NULL,
                        FlowSOMstatus = "No",
                        DiffusionMapStatus = "No")
    
    
    ## Load cytofkit RData object
    observeEvent(input$goButton, {
        cytofkitObj <- input$cytofkitObj
        if (is.null(cytofkitObj)){
            v$data <- NULL
        }else{
            cat(cytofkitObj$datapath)
            load(cytofkitObj$datapath)
            v$data <- analysis_results
            if(is.null(v$data$resultDir)){
                v$data$resultDir <- path.expand("~")  ## default save to home if not specified
            }
            if(is.null(v$data$rawFCSdir)){
                v$data$rawFCSdir <- path.expand("~")  ## default to be home if not specified
            }
            if(is.null(v$data$projectName)){
                v$data$projectName <- "cytofkit_shinyAPP_output"
            }
        }
    })
    
    ## Save and parse cytofkit RData object
    observeEvent(input$saveButton, {
        if (!is.null(v$data)){
            withProgress(message=paste0('Saving Results to', v$data$resultDir), value=0, {
                analysis_results <<- v$data
                save(analysis_results, file = paste0(v$data$resultDir, .Platform$file.sep, v$data$projectName, ".RData"))
                cytof_writeResults(analysis_results)
            })
            ## open the results directory
            opendir(v$data$resultDir)
        }
    })
    
    ## result object which will be updated by S_runFlowSOM
    observeEvent(input$S_runFlowSOM, {
        if(!is.null(v$data) && !is.null(input$s_markerSelect)){
            obj <- v$data
            withProgress(message=paste0('Runing FlowSOM using k=', input$S_FlowSOM_k), value=0, {
                FlowSOM_cluster <- cytof_cluster(xdata = obj$expressionData[ ,input$s_markerSelect],
                                                 method = "FlowSOM",
                                                 FlowSOM_k = input$S_FlowSOM_k)
            })
            
            ## update FlowSOM cluster results
            obj$clusterRes[["FlowSOM"]] <- FlowSOM_cluster
            ## update the project name
            obj$projectName <- paste0(obj$projectName, "_cytofkit_ShinyApp_Output")
            v$data <- obj
            v$FlowSOMstatus <- "No"
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
            })
            
            ## update progressionRes results
            obj$progressionRes <- diffmapRes
            v$data <- obj
            v$DiffusionMapStatus <- "Yes"
        }
    })
    
    ## update cluster labels
    observeEvent(input$updatelabel, {
        if(!is.null(v$data) && !is.null(input$m_labelCluster)){
            obj <- v$data
            clusterMethod <- input$m_labelCluster
            clusterVec<- obj$clusterRes[[clusterMethod]]
            clusterLabels <- clusterVec
            clusters <- sort(unique(clusterVec))

            for (i in 1:length(clusters)){
                clusteri <- clusters[i]
                ilabel <- input[[paste0('cluster', i)]]
                if(ilabel == ""){
                    clusterLabels[clusterLabels==clusteri] <- "Unlabelled"
                }else{
                    clusterLabels[clusterLabels==clusteri] <- ilabel
                }
            }
            ## either add new cluste or update
            labelName <- clusterMethod
            labelName <- ifelse(grepl("_Subset", labelName),
                                clusterMethod,
                                paste0(clusterMethod, "_Subset"))
            obj$clusterRes[[labelName]] <- clusterLabels
            
            v$data <- obj
        }
    })
    
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
    
    
    ##--------------------------------Summary Panel--------------------------------
    
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
    
    output$sampleSelect <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            sampleNames <- unique(sub("_[0-9]*$", "", row.names(v$data$expressionData)))
            checkboxGroupInput('samples', NULL, 
                               sampleNames, selected = sampleNames)
        }   
    })
    
    
    ##--------------------------------Cluster Plot--------------------------------
    
    output$S_PlotMethod <- renderUI({
        if(is.null(v$data) || is.null(visualizationMethods())){
            return(NULL)
        }else{
            selectInput('s_PlotMethod', 'Plot Method:', choices = visualizationMethods(), 
                        selected = visualizationMethods()[1], width = "100%")
        }   
    })
    
    output$S_PlotFunction <- renderUI({
        if(is.null(v$data) || is.null(visualizationFunctions())){
            return(NULL)
        }else{
            selectInput('s_PlotFunction', 'Plot Function:', choices = visualizationFunctions(), 
                        selected = visualizationFunctions()[1], width = "100%")
        }   
    })
    
    output$S_markerSelect <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            markerNames <- colnames(v$data$expressionData)
            checkboxGroupInput('s_markerSelect', strong('Select Markers:'),
                               markerNames, selected = markerNames, inline = TRUE)
        }   
    })
    
    output$S_ifFlowSOM <- renderUI({
        if(is.null(v$data)){
            return(NULL)
        }else{
            selectInput('s_ifFlowSOM', 'Run FlowSOM:', choices = c("Yes", "No"), 
                        selected = v$FlowSOMstatus, width = "100%")
        }   
    })
    
    output$S_ScatterPlot <- renderPlot({
        if(is.null(v$data) || is.null(input$s_PlotMethod) || is.null(input$s_PlotFunction)){
            return(NULL)
        }else{
            withProgress(message="Generating Cluster Scatter Plot", value=0, {
                gp <- scatterPlot(obj = v$data,
                                  plotMethod = input$s_PlotMethod,
                                  plotFunction = input$s_PlotFunction,
                                  pointSize = input$S_PointSize,
                                  addLabel = input$addLabel,
                                  labelSize = input$S_LabelSize,
                                  sampleLabel = input$sampleLabel,
                                  FlowSOM_k = input$S_FlowSOM_k, 
                                  selectSamples = input$samples, 
                                  facetPlot = input$facetPlot,
                                  colorPalette = input$colorPalette,
                                  labelRepel = input$labelRepel,
                                  removeOutlier = TRUE)
            })
        }
        plot(gp)
    }, height = 700, width = 750)
    
    
    ##-------------------------------Marker Plot--------------------------------
    
    output$M_plotType <- renderUI({
        s <- "Heat Map"
        radioButtons("m_plotType", NULL,
                     c("Heat Map", "Expression Map", "Marker Density", "Label Clusters"), 
                     selected = s,
                     inline = TRUE)
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
    
    ## marker expression plot
    
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
                                  sampleLabel = input$sampleLabel,
                                  FlowSOM_k = input$S_FlowSOM_k, 
                                  selectSamples = input$samples, 
                                  facetPlot = input$facetPlot,
                                  colorPalette = input$colorPalette,
                                  labelRepel = input$labelRepel,
                                  removeOutlier = TRUE)
                
            })
        }
        plot(gp)
    }, height = 800, width = 850)
    
    ## density plot
    
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
            checkboxGroupInput('m_markerSelect', strong('Select Markers:'),
                               markerNames, selected = markerNames, inline = TRUE)
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
                    samples <- sub("_[0-9]*$", "", row.names(data))
                    mySamples <- samples %in% input$samples
                    sfactors <- data.frame(do.call(cbind, v$data$clusterRes), 
                                           sample = samples, stringsAsFactors = FALSE, check.names = FALSE)
                    data <- data[mySamples, ,drop=FALSE]
                    stackFactor <- sfactors[mySamples, input$m_stackFactor]
                    
                    gp <- stackDenistyPlot(data = data, 
                                           densityCols=m_markerSelect, 
                                           stackFactor = stackFactor,
                                           kernel = "gaussian",
                                           bw = "nrd0", 
                                           adjust = 1,
                                           stackRotation = input$M_rotationDegree, 
                                           stackSeperation = "auto",
                                           x_text_size = input$M_xlab_size, 
                                           strip_text_size = input$M_markerTextSize,
                                           legend_text_size = input$M_legendTextSize, 
                                           legendRow = input$M_legendRow,
                                           legend_title = input$m_stackFactor)
                    
                })
                plot(gp)
            }
        }, height = 800, width = 850)
    })
    
    ## label cluster
    
    output$M_labelCluster <- renderUI({
        if(is.null(v$data) || is.null(v$data$clusterRes)){
            return(NULL)
        }else{
            clusterMethods <- c(names(v$data$clusterRes)) 
            #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
            selectInput('m_labelCluster', 'Choose Cluster Results to Label:', 
                        choices = clusterMethods, 
                        selected = clusterMethods[1], width = "100%")
        }   
    })
    
    
    ## currently use 100 as a limit for cluster numbers 
    ## --- TODO: use reactiveValues to automatically retrive cluster numbers --- ## 
    lapply(1:100, function(i) {
        output[[paste0('Cluster', i)]] <- renderUI({
            if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$m_labelCluster)){
                return(NULL)
            }
            
            clusters <- sort(unique(v$data$clusterRes[[input$m_labelCluster]]))
            if (i <= length(clusters)){
                x <- clusters[i]
                textInput(paste0('cluster', i), paste0('Cluster ', x," :"), 
                          value = "", width = "30%", placeholder = "Type in the cell type")
            }
        })
    })
    
    
    
    ##-------------------------------Progression Plot--------------------------------
    output$P_plotType <- renderUI({
        s <- "Run Diffusionmap"
        if(!is.null(v$data)){
            if(!is.null(v$data$progressionRes) || v$DiffusionMapStatus == "Yes")
                s <- "Subset Relationship"
        }
        
        radioButtons("p_plotType", NULL,
                     c("Subset Relationship", "Marker Expression Profile", "Run Diffusionmap"), 
                     selected = s,
                     inline = TRUE)
    })
    
    
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
                                        sampleLabel = input$sampleLabel, 
                                        labelRepel = input$labelRepel,
                                        fixCoord = FALSE)
            })
        }
        plot(gp)
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
            
            withProgress(message="Generating Marker Regression Plot", value=0, {
                data <- data.frame(v$data$progressionRes$sampleData,
                                   cluster = v$data$progressionRes$sampleCluster, 
                                   v$data$progressionRes$progressionData,
                                   check.names = FALSE)
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
            })
            plot(pp)
            
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
})



