library(shiny)
library(shinyFiles)
#currently unused

shinyUI(fluidPage(
    titlePanel("Interactive Exploration of cytofkit Analysis Results"),
    hr(),
    fluidRow(
        ## side panel--take 1/4 space
        column(3,
               h4('Load cytofkit RData:'),
               wellPanel(
                   fileInput(inputId = 'cytofkitObj',
                             label = NULL,
                             multiple = FALSE,
                             accept = c('text/RData', '.RData')),
                   shinyFilesButton('serverObj', label = "Server File Select", title = "Please select your RData", multiple = FALSE),
                   textOutput("queryText"),
                   textOutput("rdata_desc"),
                   fluidRow(
                       column(6,
                              actionButton("goButton", "Submit", icon = icon("hand-o-right"))
                       ),
                       column(6,
                              actionButton("reset", "Reset Data", icon = icon("repeat"))
                       )
                   )
               ),
               
               hr(),
               
               conditionalPanel(" input.main_panel == 'C_panel' && input.C_clusterTabs == 'C_tab1' ",
                                h4("Plot Control:"),
                                wellPanel(
                                    checkboxInput(inputId = "C_addLabel", label = "Add Cluster Labels", value = TRUE),
                                    checkboxInput(inputId = "C_labelRepel", label = "Repel Cluster Labels", value = FALSE),
                                    checkboxInput(inputId = "C_facetPlot", label = "Separate Plot by Samples", value = FALSE)
                                ),
                                actionButton("PDFClusterPlot", "Download Cluster Plot in PDF", icon = icon("download"))
               ),
               conditionalPanel(" input.main_panel == 'M_panel' && input.M_markerTabs == 'M_tab1' ",
                                h4("Plot Control:"),
                                wellPanel(
                                    selectInput('M_heatmap_dendrogram', strong('Heatmap Dendrogram:'), 
                                                choices = c("both","row","column","none"), 
                                                selected = "both", width = "100%"),  
                                    selectInput('M_heatmap_colorPalette', strong('Color Palette:'), 
                                                choices = c("bluered", "greenred", "spectral1", "spectral2"), 
                                                selected = "bluered", width = "100%")
                                ),
                                actionButton("PDFHeatmap", "Download Marker Heatmap in PDF", icon = icon("download"))
               ),
               conditionalPanel(" input.main_panel == 'M_panel' && input.M_markerTabs == 'M_tab2' ",
                                h4("Plot Control:"),
                                wellPanel(
                                    actionButton("PDFExpPlot", "Download Exp Plot in PDF", icon = icon("download"))
                                )),
               conditionalPanel(" input.main_panel == 'M_panel' && input.M_markerTabs == 'M_tab3' ",
                                h4("Plot Control:"),
                                wellPanel(
                                    actionButton("PDFHistogram", "Download Histogram in PDF", icon = icon("download"))
                                )),
               conditionalPanel(" input.main_panel == 'S_panel' && input.S_sampleTabs == 'S_tab1' ",
                                h4("Plot Control:"),
                                wellPanel(
                                    selectInput('S_heatmap_dendrogram', strong('Heatmap Dendrogram:'), 
                                                choices = c("both","row","column","none"), 
                                                selected = "both", width = "100%"),  
                                    selectInput('S_heatmap_colorPalette', strong('Color Palette:'), 
                                                choices = c("bluered", "greenred", "spectral1", "spectral2"), 
                                                selected = "bluered", width = "100%")
                                ),
                                actionButton("PDFSamHeat", "Download Sample Heatmap in PDF", icon = icon("download"))
               ),
               conditionalPanel(" input.main_panel == 'S_panel' && input.S_sampleTabs == 'S_tab2' ",
                                h4("Plot Control:"),
                                actionButton("PDFrateChange", "Download Rate Change Plot in PDF", icon = icon("download"))
               ),
               conditionalPanel(" input.main_panel == 'P_panel' && input.P_progressionTabs == 'P_tab1' ",
                                h4("Plot Control:"),
                                wellPanel(
                                    checkboxInput(inputId = "P_addLabel", label = "Add Cluster Labels", value = TRUE),
                                    checkboxInput(inputId = "P_labelRepel", label = "Repel Cluster Labels", value = FALSE),
                                    checkboxInput(inputId = "P_facetPlot", label = "Separate Plot by Samples", value = FALSE)
                                ),
                                actionButton("PDFScatter", "Download Scatterplot in PDF", icon = icon("download"))
               ),
               conditionalPanel(" input.main_panel == 'P_panel' && input.P_progressionTabs == 'P_tab2' ",
                                h4("Plot Control:"),
                                wellPanel(
                                    checkboxInput(inputId = "P_addLabel2", label = "Add Cluster Labels", value = TRUE)
                                ),
                                actionButton("PDFmarkerPlot", "Download Marker Plot in PDF", icon = icon("download"))
               ),
               br(),
               fluidRow(
                   column(6,
                          sliderInput(inputId="tab_w", label = "PDF width(in):", 
                                      min=3, max=20, value=8, width=100, ticks=FALSE)
                   ),
                   column(6, 
                          sliderInput(inputId="tab_h", label = "PDF height(in):", 
                                      min=3, max=20, value=8, width=100, ticks=FALSE)
                   )),
               
               actionButton("OpenDir", "Open download folder", icon = icon("folder")),
               
               hr(),
               h4("Sample Filter:"),
               wellPanel(uiOutput("selectAll"),
                         uiOutput("sampleSelect")),
               
               hr(),
               h4("Data Summary:"),
               wellPanel(
                   h5("Expression Data:"),
                   textOutput("summaryText1"),
                   h5("Markers used for dimension reduction and clustering:"),
                   textOutput("summaryText5"),
                   h5("Cluster Method(s):"),
                   textOutput("summaryText2"),
                   h5("Visualization Method(s):"),
                   textOutput("summaryText3"),
                   h5("Progression Method(s):"),
                   textOutput("summaryText4")
               ),
               
               hr(),
               h4("Save results:"),
               h5("Outputs to save"),
               fluidRow(
                   column(4,
                          checkboxInput(inputId = "saveFCS", label = "FCS", value = TRUE)
                   ),
                   column(4,
                          checkboxInput(inputId = "saveRData", label = "RData", value = TRUE)
                   ),
                   column(4,
                          checkboxInput(inputId = "saveCsv", label = "csv", value = FALSE)
                   )
               ),
               actionButton("saveButton", "Save Data", icon = icon("download")),
               
               hr(),
               h4(tags$a(href="mailto:jinmiao@gmail.com,a0124008@u.nus.edu?subject=[cytofkit-question]", 
                         "Contact Us")),
               imageOutput("logo", height = "60px")
        ),
        ## main panel--take 3/4 space
        column(9,
               tabsetPanel(id="main_panel", type = "pills",
                           tabPanel(title="Cluster Panel", value="C_panel", fluidPage(
                               hr(),
                               tabsetPanel(id="C_clusterTabs", type = "tabs",
                                           tabPanel(title="Cluster Plot", value="C_tab1", 
                                                    br(),
                                                    fluidRow(
                                                        column(3,
                                                               uiOutput("C_PlotMethod")
                                                        ),
                                                        column(3, 
                                                               uiOutput("C_PlotFunction")
                                                        ),
                                                        column(3,
                                                               numericInput("C_PointSize", "Point Size:", value = 1)
                                                        ),
                                                        column(3, 
                                                               numericInput("C_LabelSize", "Label Size:", value = 12)
                                                        )
                                                    ),
                                                    uiOutput("C_clusterSelect"),
                                                    hr(),
                                                    plotOutput("C_ScatterPlot", width = "100%")
                                           ),
                                           tabPanel(title="Change Cluster Color", value="C_tab2",
                                                    br(),
                                                    wellPanel(
                                                        uiOutput("C_colourCluster")
                                                    ),
                                                    hr(),
                                                    lapply(1:100, function(i) {
                                                        uiOutput(paste0('Cluster_', i, '_col'))
                                                    }),
                                                    hr(),
                                                    fluidRow(
                                                        column(3,
                                                               actionButton("C_updateClusterColor", "Update Cluster Color", 
                                                                            icon = icon("hand-o-right"), width = "100%")
                                                        ),
                                                        column(3, 
                                                               actionButton("C_revertClusterColor", "Revert to default", 
                                                                            icon = icon("hand-o-right"), width = "100%")
                                                        ),
                                                        column(6)
                                                    ),
                                                    hr()),
                                           tabPanel(title="Annotate Clusters", value="C_tab3",
                                                    br(),
                                                    wellPanel(
                                                        uiOutput("C_labelCluster"),
                                                        uiOutput("C_labelCluster_name")
                                                    ),
                                                    hr(),
                                                    lapply(1:100, function(i) {
                                                        uiOutput(paste0('Cluster', i))
                                                    }),
                                                    hr(),
                                                    actionButton("updatelabel", "Submit Cluster Label", icon = icon("hand-o-right")),
                                                    hr()),
                                           tabPanel(title="Run FlowSOM", value="C_tab4",
                                                    br(),
                                                    h4("FlowSOM Clustering Setup:"),
                                                    hr(),
                                                    wellPanel(
                                                        numericInput("C_FlowSOM_k", "Cluster k", value = 10, width = "30%"),
                                                        uiOutput("C_markerSelect")
                                                    ),
                                                    hr(),
                                                    actionButton("C_runFlowSOM", "Run FlowSOM", icon = icon("hand-pointer-o")))
                               )
                           )),
                           
                           tabPanel(title = "Marker Panel", value = "M_panel", fluidPage(
                               hr(),
                               
                               tabsetPanel(id="M_markerTabs", type = "tabs",
                                           tabPanel(title="Expression Heat Map", value="M_tab1", 
                                                    br(),
                                                    fluidRow(
                                                        column(4, 
                                                               uiOutput("M_plotCluster")
                                                        ),
                                                        column(2,
                                                               selectInput('M_plotMethod', strong('Heatmap Type:'), 
                                                                           choices = c("mean", "median"), 
                                                                           selected = "mean", width = "100%")
                                                        ),
                                                        column(2,
                                                               selectInput('M_scaleMethod', strong('Scale Data:'), 
                                                                           choices = c("none", "row", "column"), 
                                                                           selected = "none", width = "100%")
                                                        ),
                                                        column(2,
                                                               numericInput("M_rowLabelSize", "Row Label Size:", value = 1, step = 0.5)
                                                        ),
                                                        column(2, 
                                                               numericInput("M_colLabelSize", "Col Label Size:", value = 1, step = 0.5)
                                                        )
                                                    ),
                                                    fluidRow(
                                                      column(10,
                                                             uiOutput("M_heatmapmarkerSelect")
                                                      ),
                                                      column(2,
                                                             actionButton("M_heatmapSelectAll", "All Markers"),
                                                             actionButton("M_updateHeatmap", "Update Plot")
                                                      )
                                                    ),
                                                    hr(),
                                                    plotOutput("M_heatmapPlot", width = "100%")),
                                           tabPanel(title="Expression Level Plot", value="M_tab2",
                                                    br(),
                                                    fluidRow(
                                                      column(3,
                                                             uiOutput("M_PlotMethod")
                                                      ),
                                                      column(3,
                                                             numericInput("M_PointSize", "Point Size:", value = 1),
                                                             sliderInput("M_Alpha", "Transparency:", value = 1, min = 0, max = 1, step = 0.1)
                                                      ),
                                                      column(3,
                                                             selectInput('M_colorPalette', label = "Color Palette:", 
                                                                         choices = c("bluered", "spectral1", "spectral2", "heat"), 
                                                                         selected = "bluered", width = "100%")
                                                      ),
                                                      column(3,
                                                             selectInput('M_ScaleOptions', label = "Scaling Range:", 
                                                                         choices = c("Local", "Global"), 
                                                                         selected = "Local", width = "100%"),
                                                             selectInput('M_scaledData', label = "Centering:", 
                                                                         choices = c("Un-centered", "Centered"), 
                                                                         selected = "Un-centered", width = "100%")
                                                      )
                                                    ),
                                                    fluidRow(
                                                      column(10,
                                                             uiOutput("M_PlotMarker")
                                                      ),
                                                      column(2,
                                                             actionButton("M_chooseAllMarker", "All Markers"),
                                                             actionButton("M_updateExPlot", "Update Plot")
                                                      )
                                                    ),
                                                    hr(),
                                                    plotOutput("M_markerExpressionPlot", width = "100%")), 
                                           tabPanel(title="Expression Histogram", value="M_tab3", 
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
                                                    fluidRow(
                                                      column(10,
                                                             uiOutput("M_markerSelect")
                                                      ),
                                                      column(2,
                                                             actionButton("M_histSelectAll", "All Markers")
                                                      )
                                                    ),
                                                    hr(),
                                                    actionButton("M_updateDensityPlot", "Update Plot", icon = icon("hand-pointer-o")),
                                                    plotOutput("M_stackDensityPlot", width = "100%")),
                                           tabPanel(title="Update Marker Names", value="M_tab4", 
                                                    h5('Type in Your New Name for Each Marker:'),
                                                    hr(),
                                                    lapply(1:100, function(i) {
                                                        uiOutput(paste0('Marker_', i, "_name"))
                                                    }),
                                                    hr(),
                                                    actionButton("C_updateMarkerNames", "Update Marker Name", icon = icon("hand-pointer-o")))
                               )
                           )),
                           
                           tabPanel(title = "Sample Panel", value = "S_panel", fluidPage(
                               hr(),
                               tabsetPanel(id="S_sampleTabs", type = "tabs",
                                           tabPanel(title="Cell Percentage Heatmap", value="S_tab1", 
                                                    br(),
                                                    fluidRow(
                                                        column(4, 
                                                               uiOutput("S_plotCluster")
                                                        ),
                                                        column(2,
                                                               selectInput('S_plotMethod', strong('Heatmap Type:'), 
                                                                           choices = c("percentage"), 
                                                                           selected = "percentage", width = "100%")
                                                        ),
                                                        column(2,
                                                               selectInput('S_scaleMethod', strong('Scale Data:'), 
                                                                           choices = c("none", "row", "column"), 
                                                                           selected = "none", width = "100%")
                                                        ),
                                                        column(2,
                                                               numericInput("S_rowLabelSize", "Row Label Size:", value = 1, step = 0.5)
                                                        ),
                                                        column(2, 
                                                               numericInput("S_colLabelSize", "Col Label Size:", value = 1, step = 0.5)
                                                        )
                                                    ),
                                                    hr(),
                                                    plotOutput("S_heatmapPlot", width = "100%")
                                           ),
                                           tabPanel(title="Cell Percentage Line Chart", value="S_tab2", 
                                                    br(),
                                                    uiOutput("S_clusterMethod2"),
                                                    uiOutput("S_clusterFilter"),
                                                    hr(),
                                                    plotOutput("S_rateChangePlot", width = "100%")
                                           ),
                                           tabPanel(title="Regroup Samples", value="S_tab3",
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
                           
                           tabPanel(title="Progression Panel", value = "P_panel", fluidPage(
                               hr(),
                               tabsetPanel(id="P_progressionTabs", type = "tabs",
                                           tabPanel(title="Subset Relationship Plot", value="P_tab1", 
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
                                           tabPanel(title="Marker Expression Profile", value="P_tab2", 
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
                                           tabPanel(title="Run Diffusion Map", value="P_tab3",
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