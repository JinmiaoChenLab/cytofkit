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
                   checkboxInput("labelRepel", label = "Repel Cluster Labels", value = FALSE),
                   checkboxInput("sampleLabel", label = "Label Samples by Shapes", value = FALSE),
                   checkboxInput("facetPlot", label = "Seperate Plot by Samples", value = FALSE)
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
                           tabPanel("Scatter Plot", fluidPage(
                               hr(),
                               fluidRow(
                                   column(2,
                                          uiOutput("S_PlotMethod")
                                   ),
                                   column(2, 
                                          uiOutput("S_PlotFunction")
                                   ),
                                   column(2,
                                          numericInput("S_PointSize", "Point Size:", value = 3)
                                   ),
                                   column(2, 
                                          numericInput("S_LabelSize", "Label Size:", value = 12)
                                   ),
                                   column(2,
                                          selectInput('colorPalette', label = "Color Palette:", 
                                                      choices = c("bluered", "topo", "heat", "terrain", "cm"), 
                                                      selected = "bluered", width = "100%")
                                   ),
                                   column(2,
                                          uiOutput("S_ifFlowSOM")
                                   )
                               ),
                               hr(),
                               conditionalPanel("input.s_ifFlowSOM == 'No'",
                                                plotOutput("S_ScatterPlot", width = "80%")
                                                ),
                               conditionalPanel("input.s_ifFlowSOM == 'Yes'",
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
                           
                           tabPanel("Heat Map", fluidPage(
                               hr(),
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
                                                checkboxInput("P_reverseOrder", label = "Reverse Order", value = FALSE),
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
                                                
                                                    h5("Diffusionmap Parameters"),
                                                    fluidRow(
                                                        column(6,
                                                               selectInput('P_distMethod', 'Distance calculation Method:', choices = c("euclidean"), 
                                                                           selected = "euclidean", width = "100%")
                                                        ),
                                                        column(6,
                                                               numericInput("P_outDim", "Output Dimensionality:", value = 3, 
                                                                            min = 1, max = 5, step = 1, width = "100%")
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