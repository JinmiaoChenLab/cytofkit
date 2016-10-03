library(shiny)

shinyUI(fluidPage(
    titlePanel("Interactive Visualization of cytofkit Results"),
    hr(),
    fluidRow(
        column(3,
               h4('Load cytofkit RData:'),
               wellPanel(
                   fileInput(inputId = 'cytofkitObj',
                             label = NULL,
                             multiple = FALSE,
                             accept = c('text/RData', '.RData')),
                   actionButton("goButton", "Submit", icon = icon("hand-o-right"))
               ),
               
               hr(),
               h4("Plot Control:"),
               wellPanel(
                   checkboxInput("addLabel", label = "Add Cluster Labels", value = TRUE),
                   checkboxInput("labelRepel", label = "Repel Cluster Labels", value = FALSE),
                   checkboxInput("facetPlot", label = "Separate Plot by Samples", value = FALSE)
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
                           tabPanel("Cluster Panel", fluidPage(
                               hr(),
                               tabsetPanel(id="C_clusterTabs",
                                           tabPanel(title="Cluster Plot", value="C_panel1", 
                                                    br(),
                                                    fluidRow(
                                                        column(3,
                                                               uiOutput("C_PlotMethod")
                                                        ),
                                                        column(3, 
                                                               uiOutput("C_PlotFunction")
                                                        ),
                                                        column(3,
                                                               numericInput("S_PointSize", "Point Size:", value = 1)
                                                        ),
                                                        column(3, 
                                                               numericInput("S_LabelSize", "Label Size:", value = 12)
                                                        )
                                                    ),
                                                    hr(),
                                                    plotOutput("C_ScatterPlot", width = "80%")
                                                    ),
                                           tabPanel(title="Annotate Clusters", value="C_panel2",
                                                    br(),
                                                    uiOutput("C_labelCluster"),
                                                    hr(),
                                                    lapply(1:100, function(i) {
                                                        uiOutput(paste0('Cluster', i))
                                                    }),
                                                    hr(),
                                                    actionButton("updatelabel", "Submit Cluster Label", icon = icon("hand-o-right")),
                                                    hr()),
                                           tabPanel(title="Run FlowSOM", value="C_panel3",
                                                    br(),
                                                    h4("FlowSOM Clustering Setup:"),
                                                    hr(),
                                                    wellPanel(
                                                        numericInput("S_FlowSOM_k", "Cluster k", value = 10, width = "30%"),
                                                        uiOutput("C_markerSelect")
                                                    ),
                                                    hr(),
                                                    actionButton("C_runFlowSOM", "Run FlowSOM", icon = icon("hand-pointer-o")))
                                           )
                           )),
                           
                           tabPanel("Marker Panel", fluidPage(
                               hr(),
                               tabsetPanel(id="M_markerTabs",
                                   tabPanel(title="Level Plot", value="M_panel1",
                                            br(),
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
                                            plotOutput("M_markerExpressionPlot", width = "100%")), 
                                   tabPanel(title="Histogram", value="M_panel2", 
                                            br(),
                                            fluidRow(
                                                column(4,
                                                       uiOutput("M_stackFactor")
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
                                            plotOutput("M_stackDensityPlot", width = "100%")),
                                   tabPanel(title="Heat Map", value="M_panel3", 
                                            br(),
                                            fluidRow(
                                                column(4, 
                                                       uiOutput("H_plotCluster")
                                                ),
                                                column(2,
                                                       selectInput('H_plotMethod', strong('Heatmap Type:'), 
                                                                   choices = c("mean", "median"), 
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
                                            plotOutput("H_heatmapPlot", width = "100%"))
                                   )
                           )),
                           
                           tabPanel("Sample Panel", fluidPage(
                               hr(),
                               tabsetPanel(id="S_sampleTabs",
                                           tabPanel(title="Cell Counts", value="S_panel1", 
                                                    br(),
                                                    fluidRow(
                                                        column(6,
                                                               uiOutput("S_clusterMethod")
                                                        ),
                                                        column(6, 
                                                               selectInput('S_viewOption', strong('View Option:'), 
                                                                           choices = c("Count Table", "Percentage Heat Map", 
                                                                                       "Percentage Change Plot"), 
                                                                           selected = "Count Table", width = "100%")
                                                        )
                                                    ),
                                                    hr(),
                                                    conditionalPanel(" input.S_viewOption == 'Count Table' ",
                                                                     tableOutput('S_clusterTable')),
                                                    conditionalPanel(" input.S_viewOption == 'Percentage Heat Map' ",
                                                                     plotOutput("S_heatmapPlot", width = "100%")),
                                                    conditionalPanel(" input.S_viewOption == 'Percentage Change Plot' ",
                                                                     plotOutput("S_rateChangePlot", width = "100%"))
                                                    ),
                                           
                                           tabPanel(title="Group Samples", value="S_panel2",
                                                    br(),
                                                    h4("Type in the Group Name for Each Sample:"),
                                                    lapply(1:100, function(i) {
                                                        uiOutput(paste0('S_sample', i))
                                                    }),
                                                    hr(),
                                                    textInput("sampleGroupLevels", "Group Name Levels: (to order the group names)", 
                                                              value = "", width = "100%",
                                                              placeholder = "Type in group names in order, seperated by semicolon(;)"),
                                                    hr(),
                                                    fluidRow(
                                                        column(3,
                                                               actionButton("updateSampleGroups", "Submit New Sample Groups", icon = icon("hand-o-right"))
                                                        ),
                                                        column(3, 
                                                               actionButton("revertSampleNames", "Revert to Old Sample Names", icon = icon("hand-o-right"))
                                                        ),
                                                        column(6)
                                                    ),
                                                    hr())
                               )
                           )),
                           
                           tabPanel("Progression Panel", fluidPage(
                               hr(),
                               tabsetPanel(id="P_progressionTabs",
                                   tabPanel(title="Subset Relationship Plot", value="P_panel1", 
                                            br(),
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
                                                       numericInput("P_LabelSize", "Label Size:", value = 8)
                                                )
                                            ),
                                            plotOutput("P_scatterPlot", width = "80%")), 
                                   tabPanel(title="Marker Expression Profile", value="P_panel2", 
                                            br(),
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
                                                column(3,
                                                       actionButton("P_updateRegressionPlot", "Update Plot", icon = icon("hand-pointer-o"))
                                                ),
                                                column(2,
                                                       checkboxInput("P_reverseOrder", label = "Reverse Order", value = FALSE)
                                                ),
                                                column(3,
                                                       checkboxInput("P_combineTrends", label = "Combine Trend Lines", value = FALSE)
                                                ),
                                                column(4)
                                                
                                            ),
                                            plotOutput("P_markerPlot", width = "100%")), 
                                   tabPanel(title="Run Diffusion Map", value="P_panel3",
                                            br(),
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
                                            actionButton("P_runDiffusionmap", "Run Diffusionmap", icon = icon("hand-pointer-o")))
                                   )
                          )) 
              )
        )
    )
))