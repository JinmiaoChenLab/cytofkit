## max data size
options(shiny.maxRequestSize=1024^10) 

shinyServer(function(input, output, session) {
   
    ##-----------------------------------Load Data--------------------------------------
    
    ## cytofkit results RData
    resObj <- reactive({
        if (input$goButton == 0)
            return()
        isolate({
            cytofkitObj <- input$cytofkitObj
            if (is.null(cytofkitObj))
                resObj <- NULL
            cat(cytofkitObj$datapath)
            load(cytofkitObj$datapath)
            resObj <- analysis_results
            })
        return(resObj)
    })
    
    ## dimension names
    dimensionNames <- reactive({
        if(is.null(resObj()))
            return(NULL)
        dnames <- do.call(c, lapply(resObj()$dimReducedRes, colnames))
        return(dnames)
    })
    
    ## colour label
    colourLabels <- reactive({
        if(is.null(resObj()))
            return(NULL)
        clabels <- c(names(resObj()$clusterRes), colnames(resObj()$allExpressionData))
        return(clabels)
    })
    
    ## cluster methods
    clusterMethods <- reactive({
        if(is.null(resObj()))
            return(NULL)
        cMethods <- names(resObj()$clusterRes)
        return(cMethods)
    })
    
    ## progression methods
    progressionMethods <- reactive({
        if(is.null(resObj()))
            return(NULL)
        if(is.null(resObj()$progressionRes))
            return(NULL)
        progressionMethods <- colnames(resObj()$progressionRes[[3]])
        return(progressionMethods)
    })
    
    
    ##--------------------------------Creat Visualization--------------------------------
    
    output$summaryText1 <- renderText({
        if(is.null(resObj()))
            return(NULL)
        paste0("-- ", nrow(resObj()[[1]]), "cells x ", ncol(resObj()[[1]]), "markers")
    })
    
    output$summaryText2 <- renderText({
        if(is.null(resObj()))
            return(NULL)
        paste0("-- ", paste(names(resObj()$clusterRes), collapse = " | "))
    })
    
    output$summaryText3 <- renderText({
        if(is.null(resObj()))
            return(NULL)
        paste0("-- ", paste(resObj()$visualizationMethods, collapse =  " | "))
    })
    
    output$summaryText4 <- renderText({
        if(is.null(resObj()))
            return(NULL)
        paste0("-- ", ifelse(is.null(resObj()$progressionRes), "NULL", "isomap"))
    })
    
    output$sampleSelect <- renderUI({
        if(is.null(resObj())){
            return(NULL)
        }else{
            sampleNames <- unique(sub("_[0-9]*$", "", row.names(resObj()$expressionData)))
            checkboxGroupInput('samples', strong('Select Samples:'), 
                               sampleNames, selected = sampleNames)
        }   
    })
    
    output$x_choose <- renderUI({
        if(is.null(resObj()) || is.null(dimensionNames())){
            return(NULL)
        }else{
            selectInput('x_lab', 'Choose X:', choices = dimensionNames(), 
                        selected = dimensionNames()[1], width = "100%")
        }   
    })
    
    output$y_choose <- renderUI({
        if(is.null(resObj()) || is.null(dimensionNames())){
            return(NULL)
        }else{
            selectInput('y_lab', 'Choose Y:', choices = dimensionNames(), 
                        selected = dimensionNames()[2], width = "100%")
        }   
    })
    
    output$z_choose <- renderUI({
        if(is.null(resObj()) || is.null(colourLabels())){
            return(NULL)
        }else{
            selectInput('z_lab', 'Choose Z:', choices = colourLabels(), 
                        selected = colourLabels()[1], width = "100%")
        }   
    })
    
    output$cMethod_choose <- renderUI({
        if(is.null(resObj()) || is.null(clusterMethods())){
            return(NULL)
        }else{
            selectInput('clusterMethod', 'Cluster Method:', choices = clusterMethods(), 
                        selected = clusterMethods()[1], width = "100%")
        }   
    })
    
    output$pMethod_choose <- renderUI({
        if(is.null(resObj()) || is.null(progressionMethods())){
            return(NULL)
        }else{
            selectInput('progressionMethod', 'Progression Method:', choices = progressionMethods(), 
                        selected = progressionMethods()[1], width = "100%")
        }   
    })
    
    ## xyz plot
    output$xyzPlot <- renderPlot({
        if(is.null(resObj()))
            return(NULL)
        gp <- visuaPlot(obj = resObj(),
                        xlab = input$x_lab, 
                        ylab = input$y_lab, 
                        zlab = input$z_lab,
                        pointSize = input$pointSize,
                        addLabel = input$addLabel,
                        labelSize = input$labelSize,
                        sampleLabel = input$sampleLabel,
                        FlowSOM_k = input$FlowSOM_k, 
                        selectSamples = input$samples, 
                        removeOutlier = TRUE)
        plot(gp)
    }, height = 700, width = 750)
    
    ## heatmap plot
    output$heatmapPlot <- renderPlot({
        if(is.null(resObj()))
            return(NULL)
        heatMap(data=resObj(), clusterMethod=input$clusterMethod, 
                type=input$heatmapMethod, selectSamples=input$samples,
                cex_row_label=input$rLabelSize, cex_col_label=input$cLabelSize, 
                scaleMethod = input$scaleMethod)
    }, height = 800, width = 850)
    
    ## progression plot
    output$progressionPlot <- renderPlot({
        if(is.null(resObj()))
            return(NULL)
        progressionPlot(resObj(), orderCol = input$progressionMethod)
    }, height = 800, width = 850)
    
    
})



