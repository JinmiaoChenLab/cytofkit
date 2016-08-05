library(shiny)

shinyUI(fluidPage(
    titlePanel("Interactive Visualization of cytofkit Results"),
    
    fluidRow(
        column(3,
               
               fileInput('cytofkitObj', h4('cytofkit RData:'), multiple = FALSE, 
                         accept = c('text/RData', '.RData')),
               actionButton("goButton", "Submit", icon = icon("hand-o-right")),
               
               hr(),
               h4("Plot Annotation:"),
               wellPanel(
                   checkboxInput("addLabel", label = "Add Cluster Labels", value = TRUE),
                   checkboxInput("labelRepel", label = "Repel Cluster Labels", value = FALSE)
               ),
               
               hr(),
               h4("Sample Filter:"),
               wellPanel(uiOutput("sampleSelect")),
              
               hr(),
               h4("Data Summary:"),
               wellPanel(
                   h5("Expression Data:"),
                   textOutput("summaryText1"),
                   h5("Cluster Method(s):"),
                   textOutput("summaryText2"),
                   h5("Visualization Method(s):"),
                   textOutput("summaryText3"),
                   h5("Progression Method(s):"),
                   textOutput("summaryText4")
               ),
               
               hr(),
               actionButton("saveButton", "Save Results", icon = icon("download")),
               
               hr(),
               div(style = "margin-top: 30px; width: 200px; ", HTML("Developed by")),
               div(style = "margin-top: 10px; ", 
                   HTML("<img style='width: 150px;' src='http://archild.sign.a-star.edu.sg/images/logo.png'>"))
        ),
        column(9,
               tabsetPanel(type = "pills",
                           tabPanel("Cluster Plot", fluidPage(
                               hr(),
                               fluidRow(
                                   column(2,
                                          uiOutput("S_PlotMethod")
                                   ),
                                   column(4, 
                                          uiOutput("S_PlotFunction")
                                   ),
                                   column(2,
                                          numericInput("S_PointSize", "Point Size:", value = 1)
                                   ),
                                   column(2, 
                                          numericInput("S_LabelSize", "Label Size:", value = 12)
                                   ),
                                   column(2,
                                          uiOutput("S_ifFlowSOM")
                                   )
                               ),
                               conditionalPanel("input.s_ifFlowSOM == 'No'",
                                                fluidRow(
                                                    column(3,
                                                           checkboxInput("sampleLabel", label = "Label Samples by Shapes", value = FALSE)
                                                    ),
                                                    column(3,
                                                           checkboxInput("facetPlot", label = "Seperate Plot by Samples", value = FALSE)
                                                    ),
                                                    column(6)
                                                ),
                                                hr(),
                                                plotOutput("S_ScatterPlot", width = "80%")
                                                ),
                               conditionalPanel("input.s_ifFlowSOM == 'Yes'",
                                                hr(),
                                                h4("FlowSOM Clustering Setup:"),
                                                hr(),
                                                wellPanel(
                                                    numericInput("S_FlowSOM_k", "Cluster k", value = 10, width = "30%"),
                                                    uiOutput("S_markerSelect")
                                                ),
                                                hr(),
                                                actionButton("S_runFlowSOM", "Run FlowSOM", icon = icon("hand-pointer-o"))
                               )
                           )),
                           
                           tabPanel("Marker Plot", fluidPage(
                               hr(),
                               wellPanel(
                                   uiOutput("M_plotType")
                               ),
                               hr(),
                               conditionalPanel("input.m_plotType == 'Heat Map'",
                                                fluidRow(
                                                    column(4, 
                                                           uiOutput("H_plotCluster")
                                                    ),
                                                    column(2,
                                                           selectInput('H_plotMethod', strong('Heatmap Type:'), 
                                                                       choices = c("mean", "median", "percentage"), 
                                                                       selected = "mean", width = "100%")
                                                    ),
                                                    column(2,
                                                           selectInput('H_scaleMethod', strong('Scale Data:'), 
                                                                       choices = c("none", "row", "column"), 
                                                                       selected = "none", width = "100%")
                                                    ),
                                                    column(2,
                                                           numericInput("H_rowLabelSize", "Row Label Size:", value = 1, step = 0.5)
                                                    ),
                                                    column(2, 
                                                           numericInput("H_colLabelSize", "Col Label Size:", value = 1, step = 0.5)
                                                    )
                                                ),
                                                hr(),
                                                plotOutput("H_heatmapPlot", width = "100%")
                                                ),
                               conditionalPanel("input.m_plotType == 'Expression Map'",
                                                fluidRow(
                                                    column(3,
                                                           uiOutput("M_PlotMethod")
                                                    ),
                                                    column(3, 
                                                           uiOutput("M_PlotMarker")
                                                    ),
                                                    column(3,
                                                           numericInput("M_PointSize", "Point Size:", value = 1)
                                                    ),
                                                    column(3,
                                                           selectInput('colorPalette', label = "Color Palette:", 
                                                                       choices = c("bluered", "topo", "heat", "terrain", "cm"), 
                                                                       selected = "bluered", width = "100%")
                                                    )
                                                ),
                                                hr(),
                                                plotOutput("M_markerExpressionPlot", width = "100%")
                                                ),
                               conditionalPanel("input.m_plotType == 'Marker Density'",
                                                fluidRow(
                                                    column(2,
                                                           uiOutput("M_stackFactor")
                                                    ),
                                                    column(2, 
                                                           numericInput("M_rotationDegree", "Rotation Degree:", 
                                                                        value = 0, step = 1, min=0, max=90)
                                                    ),
                                                    column(2,
                                                           numericInput("M_markerTextSize", "Marker Text Size:", 
                                                                        value = 12, step = 1, min=1, max=15)
                                                    ),
                                                    column(2,
                                                           numericInput("M_xlab_size", "x Label Size:", 
                                                                        value = 2, step = 1, min=1, max=10)
                                                    ),
                                                    column(2,
                                                           numericInput("M_legendTextSize", "Legend Size:", 
                                                                        value = 1, step = 0.5, min=1, max=10)
                                                    ),
                                                    column(2,
                                                           numericInput("M_legendRow", "Legend Row:", 
                                                                        value = 2, step = 1, min=1, max=10)
                                                    )
                                                ),
                                                uiOutput("M_markerSelect"),
                                                hr(),
                                                actionButton("M_updateDensityPlot", "Update Plot", icon = icon("hand-pointer-o")),
                                                plotOutput("M_stackDensityPlot", width = "100%")
                                                
                                                ),
                               conditionalPanel("input.m_plotType == 'Label Clusters'",
                                                uiOutput("M_labelCluster"),
                                                hr(),
                                                lapply(1:100, function(i) {
                                                    uiOutput(paste0('Cluster', i))
                                                }),
                                                hr(),
                                                actionButton("updatelabel", "Submit Cluster Label", icon = icon("hand-o-right")),
                                                hr()
                                                )
                           )),
                           
                           tabPanel("Subset Progression", fluidPage(
                               hr(),
                               wellPanel(
                                   uiOutput("P_plotType")
                               ),
                               hr(),
                               conditionalPanel("input.p_plotType == 'Subset Relationship'",
                                                fluidRow(
                                                    column(3,
                                                           uiOutput("P_xlab")
                                                    ),
                                                    column(3,
                                                           uiOutput("P_ylab")
                                                    ),
                                                    column(3,
                                                           numericInput("P_PointSize", "Point Size:", value = 3)
                                                    ),
                                                    column(3, 
                                                           numericInput("P_LabelSize", "Label Size:", value = 12)
                                                    )
                                                ),
                                                plotOutput("P_scatterPlot", width = "80%")
                                                
                               ),
                               conditionalPanel("input.p_plotType == 'Marker Expression Profile'",
                                                fluidRow(
                                                    column(3,
                                                           uiOutput("P_orderBy")
                                                    ),
                                                    column(2,
                                                           numericInput("P_LabelSize2", "Label Size:", value = 5)
                                                    ),
                                                    column(7,
                                                           uiOutput("P_clusterSelect")
                                                    )
                                                ),
                                                hr(),
                                                uiOutput("P_markerSelect"),
                                                hr(),
                                                fluidRow(
                                                    column(2,
                                                           checkboxInput("P_reverseOrder", label = "Reverse Order", value = FALSE)
                                                    ),
                                                    column(3,
                                                           checkboxInput("P_combineTrends", label = "Combine Trend Lines", value = FALSE)
                                                    ),
                                                    column(3,
                                                           actionButton("P_updateRegressionPlot", "Update Plot", icon = icon("hand-pointer-o"))
                                                    ),
                                                    column(4)
                                                    
                                                ),
                                                plotOutput("P_markerPlot", width = "100%")
                                                ),
                               
                               conditionalPanel("input.p_plotType == 'Run Diffusionmap'",
                                                h4("Diffusionmap Setup:"),
                                                
                                                wellPanel(
                                                    h5("Cluster-based down-sampling to remove subset aboundance heterogeneity"),
                                                    
                                                    fluidRow(
                                                        column(4,
                                                               uiOutput("P_clusterMethod")
                                                        ),
                                                        column(4,
                                                               numericInput("P_clusterSampleSize", "Cluster Sample Size", value = 500, 
                                                                            min = 10, max = 1000, step = 5, width = "100%")
                                                        ),
                                                        column(4,
                                                               selectInput('P_sampleMethod', 'Downsample Method:', choices = c("ceil", "all", "fixed", "min"), 
                                                                           selected = "ceil", width = "100%")
                                                        )
                                                    ),
                                                    
                                                    tableOutput('P_clusterTable'),
                                                    
                                                    uiOutput("P_clusterFilter"),
                                                    hr(),
                                                
                                                    h5("Diffusionmap Parameters"),
                                                    fluidRow(
                                                        column(6,
                                                               selectInput('P_distMethod', 'Distance calculation Method:', choices = c("euclidean"), 
                                                                           selected = "euclidean", width = "100%")
                                                        ),
                                                        column(6,
                                                               numericInput("P_outDim", "Output Dimensionality:", value = 4, 
                                                                            min = 1, max = 6, step = 1, width = "100%")
                                                        )
                                                    )
                                                ),
                                                hr(),
                                                actionButton("P_runDiffusionmap", "Run Diffusionmap", icon = icon("hand-pointer-o"))
                        )
                    )
                ) 
            )
        )
    )
))