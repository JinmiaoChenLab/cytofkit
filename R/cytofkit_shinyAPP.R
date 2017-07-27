#' A Shiny APP to interactively visualize the analysis results 
#' 
#' Take the the RData object file saved by cytofkit as input, automatically load the data and allow exploration of the analysis results with interactive control
#'
#'
#' @param RData Either the RData object file or data object, if missing, RData file need to be loaded on the ShinyAPP
#' 
#' @return Opens shinyApp session for data visualisation
#' @import shiny
#' @import shinyFiles
#' @importFrom grDevices dev.copy2pdf
#' @importFrom graphics plot
#' @author Hao Chen
#' @export
#' @examples 
#' d <- system.file('extdata', package = 'cytofkit')
#' Rdata <- list.files(d, pattern = '.RData$', full.names = TRUE)
#' #only for interactive sessions, remove hash to run
#' #cytofkitShinyAPP(Rdata)
cytofkitShinyAPP <- function(RData = NULL) {
    
    source(system.file('shiny', "global.R", package = 'cytofkit'))
  
    analysis_results <- NULL
    sampleInformation <- NULL
    progCluster <- NULL
    serverObj <- NULL
    roots <- c(wd=getwd())
    
    if(!missing(RData)){
        if(class(RData) == "character"){
            if(file.exists(RData)){
              if(tools::file_ext(RData) == "RData"){
                load(RData)
                direct_analysis_results <- analysis_results
                message(".RData loaded!")
              }else{
                stop("Argument is not .RData file!")
              }
            }else{
                stop("RData file doesn't exist! Please check your obj file")
            }
        }else{
            analysis_results <- RData
        }
        
        if(is.null(analysis_results$projectName)){
            analysis_results$projectName <- "cytofkit_shinyAPP_output"
        }
        
        if(!is.null(analysis_results$progressionRes)){
            ## default the first cluster results are used for progression analysis
            progCluster <- names(analysis_results$clusterRes)[1]
        }
        
        sampleInformation <- data.frame(cellID = row.names(analysis_results$expressionData),
                                        cellSample = factor(sub("_[0-9]*$", "", row.names(analysis_results$expressionData))),
                                        stringsAsFactors = FALSE)
        analysis_results$sampleInfo <- sampleInformation
    }
    
    #shiny::runApp(system.file('shiny', package = 'cytofkit'))
    options(launch.browser = TRUE)
    shinyApp(
        ui = fluidPage(
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
                                    wellPanel(
                                      actionButton("PDFClusterPlot", "Download Cluster Plot in PDF", icon = icon("download")),
                                      hr(),
                                      fluidRow(
                                        column(6,
                                               sliderInput(inputId="C_tab1_w", label = "PDF width(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ),
                                        column(6, 
                                               sliderInput(inputId="C_tab1_h", label = "PDF height(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ))
                                    )),
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
                                    wellPanel(
                                      actionButton("PDFHeatmap", "Download Marker Heatmap in PDF", icon = icon("download")),
                                      hr(),
                                      fluidRow(
                                        column(6,
                                               sliderInput(inputId="M_tab3_w", label = "PDF width(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ),
                                        column(6, 
                                               sliderInput(inputId="M_tab3_h", label = "PDF height(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ))
                                    )),
                   conditionalPanel(" input.main_panel == 'M_panel' && input.M_markerTabs == 'M_tab2' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      actionButton("PDFExpPlot", "Download Exp Plot in PDF", icon = icon("download")),
                                      hr(),
                                      fluidRow(
                                        column(6,
                                               sliderInput(inputId="M_tab1_w", label = "PDF width(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ),
                                        column(6, 
                                               sliderInput(inputId="M_tab1_h", label = "PDF height(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ))
                                    )),
                   conditionalPanel(" input.main_panel == 'M_panel' && input.M_markerTabs == 'M_tab3' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      actionButton("PDFHistogram", "Download Histogram in PDF", icon = icon("download")),
                                      hr(),
                                      fluidRow(
                                        column(6,
                                               sliderInput(inputId="M_tab2_w", label = "PDF width(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ),
                                        column(6, 
                                               sliderInput(inputId="M_tab2_h", label = "PDF height(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ))
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
                                    wellPanel(
                                      actionButton("PDFSamHeat", "Download Sample Heatmap in PDF", icon = icon("download")),
                                      hr(),
                                      fluidRow(
                                        column(6,
                                               sliderInput(inputId="S_tab1_w", label = "PDF width(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ),
                                        column(6, 
                                               sliderInput(inputId="S_tab1_h", label = "PDF height(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ))
                                    )),
                   conditionalPanel(" input.main_panel == 'S_panel' && input.S_sampleTabs == 'S_tab2' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      actionButton("PDFrateChange", "Download Rate Change Plot in PDF", icon = icon("download")),
                                      hr(),
                                      fluidRow(
                                        column(6,
                                               sliderInput(inputId="S_tab2_w", label = "PDF width(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ),
                                        column(6, 
                                               sliderInput(inputId="S_tab2_h", label = "PDF height(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ))
                                    )),
                   conditionalPanel(" input.main_panel == 'P_panel' && input.P_progressionTabs == 'P_tab1' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      checkboxInput(inputId = "P_addLabel", label = "Add Cluster Labels", value = TRUE),
                                      checkboxInput(inputId = "P_labelRepel", label = "Repel Cluster Labels", value = FALSE),
                                      checkboxInput(inputId = "P_facetPlot", label = "Separate Plot by Samples", value = FALSE)
                                    ),
                                    wellPanel(
                                      actionButton("PDFScatter", "Download Scatterplot in PDF", icon = icon("download")),
                                      hr(),
                                      fluidRow(
                                        column(6,
                                               sliderInput(inputId="P_tab1_w", label = "PDF width(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ),
                                        column(6, 
                                               sliderInput(inputId="P_tab1_h", label = "PDF height(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ))
                                    )),
                   conditionalPanel(" input.main_panel == 'P_panel' && input.P_progressionTabs == 'P_tab2' ",
                                    h4("Plot Control:"),
                                    wellPanel(
                                      checkboxInput(inputId = "P_addLabel2", label = "Add Cluster Labels", value = TRUE)
                                    ),
                                    wellPanel(
                                      actionButton("PDFmarkerPlot", "Download Marker Plot in PDF", icon = icon("download")),
                                      hr(),
                                      fluidRow(
                                        column(6,
                                               sliderInput(inputId="P_tab2_w", label = "PDF width(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ),
                                        column(6, 
                                               sliderInput(inputId="P_tab2_h", label = "PDF height(in):", 
                                                           min=3, max=20, value=8, width=100, ticks=FALSE)
                                        ))
                                    )),
                   
                   actionButton("OpenDir", "Open download folder", icon = icon("folder")),
                   
                   hr(),
                   h4("Sample Filter:"),
                   wellPanel(uiOutput("sampleSelect")),
                   
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
                   actionButton("saveButton", "Save Data", icon = icon("download")),
                   
                   hr(),
                   h4(tags$a(href="mailto:chen_hao@immunol.a-star.edu.sg,Chen_Jinmiao@immunol.a-star.edu.sg?subject=[cytofkit-question]", 
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
                                                      uiOutput("M_heatmapmarkerSelect"),
                                                      hr(),
                                                      plotOutput("M_heatmapPlot", width = "100%")),
                                             tabPanel(title="Expression Level Plot", value="M_tab2",
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
                                                               selectInput('M_colorPalette', label = "Color Palette:", 
                                                                           choices = c("bluered", "spectral1", "spectral2", "heat"), 
                                                                           selected = "bluered", width = "100%")
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
                                                      uiOutput("M_markerSelect"),
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
        ),
        
        server = function(input, output, session) {
          
          ##------------------Reactive Values and Reactive Objects-------------------
          
          #if?
          v <- reactiveValues(data = NULL, sampleInfo = NULL)
          c <- reactiveValues(clusterCol = list())
          p <- reactiveValues(progressionCluster = NULL)
          
          if(!is.null(analysis_results)) {
            v$data <- analysis_results
            v$sampleInfo <- data.frame(cellID = row.names(analysis_results$expressionData),
                                       cellSample = factor(sub("_[0-9]*$", "", row.names(analysis_results$expressionData))),
                                       stringsAsFactors = FALSE)
            p$progressionCluster <- names(analysis_results$clusterRes)[1]
          }
          
          ## Scatter plot methods
          visualizationMethods <- reactive({
            if(is.null(v$data) || is.null(v$data$visualizationMethods)){
              return(NULL)
            }else{
              return(v$data$visualizationMethods)
            }
          })
          
          ## Scatter plot functions
          visualizationFunctions <- reactive({
            if(is.null(v$data) || is.null(v$data$clusterRes)){
              return(NULL)
            }else{
              return(c(names(v$data$clusterRes), 
                       "Sample",
                       "Density",
                       "None"))
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
              
              if(!is.null(v$data$progressionRes)){
                ## default the first cluster results are used for progression analysis
                p$progressionCluster <- names(v$data$clusterRes)[1]
              }
              
              
              # Need modification later
              # currently doesn't update sampleInfo with v$data$sampleInfo
              v$sampleInfo <- data.frame(cellID = row.names(v$data$expressionData),
                                         cellSample = factor(sub("_[0-9]*$", "", row.names(v$data$expressionData))),
                                         stringsAsFactors = FALSE)
              v$data$sampleInfo <- v$sampleInfo
            }
          })
          
          ## For user, set roots option to your server directory 
          shinyFileChoose(input, 'serverObj', session = session, roots = roots, filetypes = "RData")
          
          observeEvent(input$serverObj, {
            inServer <- parseFilePaths(roots= roots, input$serverObj)
            print(inServer$datapath)
            load(as.character(inServer$datapath))
            v$data <- analysis_results
              if(is.null(v$data$projectName)){
                v$data$projectName <- "cytofkit_shinyAPP_output"
              }
              if(!is.null(v$data$progressionRes)){
                ## default the first cluster results are used for progression analysis
                p$progressionCluster <- names(v$data$clusterRes)[1]
              }
              # Need modification later
              # currently doesn't update sampleInfo with v$data$sampleInfo
              v$sampleInfo <- data.frame(cellID = row.names(v$data$expressionData),
                                         cellSample = factor(sub("_[0-9]*$", "", row.names(v$data$expressionData))),
                                         stringsAsFactors = FALSE)
              v$data$sampleInfo <- v$sampleInfo
            })
          
          output$rdata_desc <- renderText({
            if(is.null(v$data)){
              paste0("No .RData loaded yet")
            }else{
              paste0("Loaded: ", v$data$resultDir, v$data$projectName, ".RData")
            }
          })
          
          observeEvent(input$reset, {
            analysis_results <- NULL
            session$reload()
            print("Reset done")
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
            paste0("-- ", nrow(v$data[[1]]), " cells x ", ncol(v$data[[1]]), " markers")
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
          
          output$summaryText5 <- renderText({
            if(is.null(v$data))
              return(NULL)
            paste0("-- ", paste(v$data$dimRedMarkers, collapse =  " | "))
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
                  warning(paste0("Path for original FCS files doesn't exist, 
                                 data cannnot be saved to new copies of FCS files.", 
                                 "Please check path: ", v$data$rawFCSdir))
                }
                
                ## NOTE: if samples are regrouped, then new FCS file cannot be saved
                incProgress(1/2, message = paste0("To ", v$data$resultDir))
                v$data$sampleInfo <- v$sampleInfo
                analysis_results <<- v$data
                cytof_writeResults(analysis_results,
                                   saveToRData = TRUE,
                                   saveToFCS = saveToFCS,
                                   saveToFiles = FALSE)
                incProgress(1/2)
                ## open the results directory
                opendir(v$data$resultDir)
              })
              }
          })
          
          observeEvent(input$OpenDir, {
            pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
            if(dir.exists(pdfDir)){
              opendir(pdfDir)
            }else{
              stop("PDF not created yet!")
            }
          })
          
          output$logo <- renderImage({
            return(list(
              src = "vignettes/logo.png",
              contentType = "image/png",
              alt = "Singapore Immunology Network"
            ))
          }, deleteFile = FALSE)
          
          ##------------------------------Cluster Panel------------------------------
          
          ##-----cluster plot-----
          output$C_PlotMethod <- renderUI({
            if(is.null(v$data) || is.null(visualizationMethods())){
              return(NULL)
            }else{
              selectInput('c_PlotMethod', 'Visualization Method:', choices = visualizationMethods(), 
                          selected = visualizationMethods()[1], width = "100%")
            }   
          })
          
          output$C_PlotFunction <- renderUI({
            if(is.null(v$data) || is.null(visualizationFunctions())){
              return(NULL)
            }else{
              selectInput('c_PlotFunction', 'Cluster By:', choices = visualizationFunctions(), 
                          selected = visualizationFunctions()[1], width = "100%")
            }   
          })
          
          output$C_markerSelect <- renderUI({
            if(is.null(v$data)){
              return(NULL)
            }else{
              markerNames <- colnames(v$data$expressionData)
              markerNames <- markerNames[order(markerNames)]
              checkboxGroupInput('c_markerSelect', strong('Select Markers:'),
                                 markerNames, selected = markerNames, inline = TRUE)
            }   
          })
          
          output$C_clusterSelect <- renderUI({
            if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_PlotFunction))
              return(NULL)
            if(input$c_PlotFunction %in% c("Sample", "Density","None")){
              return(NULL)
            }else{
              clusterMethod <- input$c_PlotFunction
              clusterIDs <- sort(unique(v$data$clusterRes[[clusterMethod]]))
              selectizeInput('c_clusterSelect', 'Clusters Filter:', 
                             choices = clusterIDs, selected = clusterIDs, 
                             multiple = TRUE, width = "100%")
              # checkboxGroupInput('p_clusterSelect', strong('Select Clusters:'), 
              #                    clusterIDs, selected = clusterIDs, inline = TRUE)
            }   
          })
          
          ## Complex dependencies here: --> (depends on)
          ## C_ScatterPlotInput --> c_PlotMethod + c_clusterSelect 
          ## c_clusterSelect --> c_PlotMethod
          ## carefull checkings are applied to solve concurrency conflicts
          C_ScatterPlotInput <- function(){
            if(is.null(v$data) || is.null(input$c_PlotMethod) || 
               is.null(input$c_PlotFunction) || is.null(input$c_clusterSelect)){
              return(NULL)
            }else if(!all(input$c_clusterSelect %in% v$data$clusterRes[[input$c_PlotFunction]]) &&
                     !(input$c_PlotFunction %in% c("Sample", "Density","None"))){
              return(NULL)
            }else{
              
              withProgress(message="Generating Cluster Scatter Plot", value=0, {
                if(input$c_PlotFunction %in% c("Sample", "Density", "None")){
                  clusterSelect <- NULL
                  clusterColor <- NULL
                }else{
                  clusterSelect <- input$c_clusterSelect
                  clusterMethod <- input$c_PlotFunction
                  if(!is.null(c$clusterCol[[clusterMethod]])){
                    clusterColor <- c$clusterCol[[clusterMethod]]
                  }else{
                    cluster_num <- length(unique(v$data$clusterRes[[clusterMethod]]))
                    clusterColor <- rainbow(cluster_num)
                  }
                }
                gp <- scatterPlot(obj = v$data,
                                  plotMethod = input$c_PlotMethod,
                                  plotFunction = input$c_PlotFunction,
                                  pointSize = input$C_PointSize,
                                  addLabel = input$C_addLabel,
                                  labelSize = input$C_LabelSize,
                                  sampleLabel = FALSE,
                                  FlowSOM_k = input$C_FlowSOM_k, 
                                  selectCluster = clusterSelect,
                                  selectSamples = input$samples, 
                                  facetPlot = input$C_facetPlot,
                                  labelRepel = input$C_labelRepel,
                                  removeOutlier = TRUE,
                                  clusterColor = clusterColor)
                incProgress(1/2)
                plot(gp)
                incProgress(1/2)
              })
            }
          }
          
          output$C_ScatterPlot <- renderPlot({
            C_ScatterPlotInput()
          }, height = 900, width = 950)
          
          observeEvent(input$PDFClusterPlot, {
            if(!is.null(v$data)){
              withProgress(message="Downloading Clusterplot PDF files...", value=0, {
                print(getwd())
                dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
                filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Clusterplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                  filename1 <- paste0(pdfDir, .Platform$file.sep,
                                      "cytofkit_shinyAPP_Clusterplot_",
                                      Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                  i = i + 1;
                }
                pdf(filename1, 
                    width=as.numeric(input$H_tab1_w), 
                    height=as.numeric(input$H_tab1_h))
                C_ScatterPlotInput()
                dev.off()
              })
            }
          })
          

          ##----- change cluster colour -----
          output$C_colourCluster <- renderUI({
            if(is.null(v$data) || is.null(v$data$clusterRes)){
              return(NULL)
            }else{
              clusterMethods <- c(names(v$data$clusterRes)) 
              #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
              selectInput('c_colourCluster', 'Choose Cluster to Change the Colour :', 
                          choices = clusterMethods, 
                          selected = clusterMethods[1], width = "50%")
            }   
          })
          
          ## currently use 100 as a limit for cluster numbers 
          ## --- TODO: use reactiveValues to automatically retrive cluster numbers --- ## 
          lapply(1:100, function(i) {
            output[[paste0('Cluster_', i, "_col")]] <- renderUI({
              if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_colourCluster)){
                return(NULL)
              }
              
              clusters <- v$data$clusterRes[[input$c_colourCluster]]
              clusterLabel <- levels(as.factor(clusters))
              if(is.null(c$clusterCol[[input$c_colourCluster]])){
                clusterColor <- rainbow(length(unique(clusters)))
              }else{
                clusterColor <- c$clusterCol[[input$c_colourCluster]]
              }
              
              if (i <= length(clusterLabel)){
                x <- clusterLabel[i]
                colourInput(inputId=paste0('cluster_', i, '_col'), 
                            label=paste0('Cluster ', x," Colour :"), 
                            value = clusterColor[i], showColour = "both", 
                            palette = "square")
              }
            })
          })
          
          ## update cluster color
          observeEvent(input$C_updateClusterColor, {
            if(!is.null(v$data) && !is.null(input$c_colourCluster)){
              clusterMethod <- input$c_colourCluster
              clusterVec<- v$data$clusterRes[[clusterMethod]]
              clusters <- levels(as.factor(clusterVec))
              clusterCols <- NULL
              for (i in 1:length(clusters)){
                clusteri <- clusters[i]
                iCol <- input[[paste0('cluster_', i, '_col')]]
                clusterCols <- c(clusterCols, iCol)
              }
              
              ## update new cluster colours
              c$clusterCol[[clusterMethod]] <- clusterCols
              
              ## jump to C_tab1
              updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
            }
          })
          
          ## revert default cluster colors
          observeEvent(input$C_revertClusterColor, {
            if(!is.null(v$data) && !is.null(input$c_colourCluster)){
              clusterMethod <- input$c_colourCluster
              c$clusterCol[[clusterMethod]] <- NULL
              
              ## jump to C_tab1
              updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
            }
          })
          
          
          ## ------annotate clusters-----
          output$C_labelCluster <- renderUI({
            if(is.null(v$data) || is.null(v$data$clusterRes)){
              return(NULL)
            }else{
              clusterMethods <- c(names(v$data$clusterRes)) 
              #clusterMethods <- clusterMethods[!grepl("Subset", clusterMethods)]
              selectInput('c_labelCluster', 'Choose Cluster Results to Annotate:', 
                          choices = clusterMethods, 
                          selected = clusterMethods[1], width = "50%")
            }   
          })
          
          output$C_labelCluster_name <- renderUI({
            if(is.null(v$data) || is.null(v$data$clusterRes) || is.null(input$c_labelCluster)){
              return(NULL)
            }else{
              textInput("c_labelCluster_name", label = "Type In Your Name for Annotated Cluster", 
                        value = paste0("Annotated_", input$c_labelCluster), width = "50%")
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
            if(!is.null(v$data) && !is.null(input$c_labelCluster) && !is.null(input$c_labelCluster_name)){
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
              
              ## update new cluster results
              labelName <- input$c_labelCluster_name
              obj$clusterRes[[labelName]] <- clusterLabels
              
              ## update the project name
              obj$projectName <- paste0(obj$projectName, "_annotated_", clusterMethod)
              
              v$data <- obj
              
              ## jump to C_tab1
              updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
            }
          })
          
          
          
          ##-----RUN flowSOM-----
          ## result object which will be updated by C_runFlowSOM
          observeEvent(input$C_runFlowSOM, {
            if(!is.null(v$data) && !is.null(input$c_markerSelect)){
              obj <- v$data
              withProgress(message=paste0('Running FlowSOM using k=', input$C_FlowSOM_k), value=0, {
                FlowSOM_cluster <- cytof_cluster(xdata = obj$expressionData[ ,input$c_markerSelect],
                                                 method = "FlowSOM",
                                                 FlowSOM_k = input$C_FlowSOM_k)
                incProgress(1/2)
                ## update FlowSOM cluster results
                obj$clusterRes[["FlowSOM"]] <- FlowSOM_cluster
                ## update the project name
                obj$projectName <- paste0(obj$projectName, "_add_FlowSOM")
                v$data <- obj
                incProgress(1/2)
              })
              
              ## jump to C_tab1
              updateTabsetPanel(session, "C_clusterTabs", selected = "C_tab1")
            }
          })
          
          
          ##------------------------------Marker Panel-------------------------------
          
          ##-----heat map plot-----
          output$M_plotCluster <- renderUI({
            if(is.null(v$data) || is.null(clusterMethods())){
              return(NULL)
            }else{
              selectInput('m_plotCluster', 'Cluster Method:', choices = clusterMethods(), 
                          selected = clusterMethods()[1], width = "100%")
            }   
          })
          
          output$M_heatmapmarkerSelect <- renderUI({
            if(is.null(v$data)){
              return(NULL)
            }else{
              sorted_markerNames <- colnames(v$data$expressionData)
              markerNames <- sorted_markerNames[order(sorted_markerNames)]
              initNum <- ifelse(length(markerNames) >=4, 4, 1)
              selectizeInput('m_heatmapmarkerSelect', 'Select Markers:', 
                             choices = markerNames, selected = markerNames[1:initNum], 
                             multiple = TRUE, width = "100%")
            }   
          })
          
          M_heatmapPlotInput <- reactive({
            if(is.null(v$data) || is.null(input$m_plotCluster))
              return(NULL)
            heatMap(data = v$data, 
                    clusterMethod = input$m_plotCluster, 
                    type = input$M_plotMethod, 
                    dendrogram = input$M_heatmap_dendrogram,
                    colPalette = input$M_heatmap_colorPalette,
                    selectSamples = input$samples,
                    selectMarkers = input$m_heatmapmarkerSelect,
                    cex_row_label= input$M_rowLabelSize, 
                    cex_col_label= input$M_colLabelSize, 
                    scaleMethod = input$M_scaleMethod)
          })
          
          output$M_heatmapPlot <- renderPlot({
            M_heatmapPlotInput()
          }, height = 900, width = 950)
          
          observeEvent(input$PDFHeatmap, {
            if(!is.null(v$data)){
              withProgress(message="Downloading Marker Heatmap PDF files...", value=0, {
                print(getwd())
                dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
                filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Marker_Heatmap_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                  filename1 <- paste0(pdfDir, .Platform$file.sep,
                                      "cytofkit_shinyAPP_Marker_Heatmap_",
                                      Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                  i = i + 1;
                }
                file.copy("cytofkit_shinyAPP_marker_heatmap_plot.pdf", filename1)
              })
            }
          })
          
          
          ##-----level plot-----
          output$M_PlotMethod <- renderUI({
            if(is.null(v$data) || is.null(visualizationMethods())){
              return(NULL)
            }else{
              selectInput('m_PlotMethod', 'Visualization Method:', choices = visualizationMethods(), 
                          selected = visualizationMethods()[1], width = "100%")
            }   
          })
          
          output$M_PlotMarker <- renderUI({
            if(is.null(v$data)){
              return(NULL)
            }else{
              sorted_markers <- colnames(v$data$expressionData)
              sorted_markers <- sorted_markers[order(sorted_markers)]
              markers <- c(sorted_markers, "All Markers", "All Markers(scaled)")
              selectInput('m_PlotMarker', 'Plot Marker:', choices = markers, 
                          selected = markers[1], width = "100%")
            }   
          })
          
          M_markerExpressionPlotInput <- function(){
            if(is.null(v$data) || is.null(input$m_PlotMethod) || is.null(input$m_PlotMarker)){
              return(NULL)
            }else{
              withProgress(message="Generating Marker Expression Plot", value=0, {
                gp <- scatterPlot(obj = v$data,
                                  plotMethod = input$m_PlotMethod,
                                  plotFunction = input$m_PlotMarker,
                                  pointSize = input$M_PointSize,
                                  addLabel = FALSE,
                                  labelSize = input$S_LabelSize,
                                  sampleLabel = FALSE,
                                  FlowSOM_k = input$C_FlowSOM_k, 
                                  selectSamples = input$samples, 
                                  facetPlot = FALSE,
                                  colorPalette = input$M_colorPalette,
                                  labelRepel = FALSE,
                                  removeOutlier = TRUE)
                incProgress(1/2)
                plot(gp)
                incProgress(1/2)
              })
            }
          }
          
          output$M_markerExpressionPlot <- renderPlot({
            M_markerExpressionPlotInput()
          }, height = 900, width = 950)
          
          observeEvent(input$PDFExpPlot, {
            if(!is.null(v$data)){
              withProgress(message="Downloading Marker Expression Plot PDF files...", value=0, {
                print(getwd())
                dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
                filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Marker_Expression_Plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                  filename1 <- paste0(pdfDir, .Platform$file.sep,
                                      "cytofkit_shinyAPP_Marker_Expression_Plot_",
                                      Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                  i = i + 1;
                }
                pdf(filename1, 
                    width=as.numeric(input$H_tab1_w), 
                    height=as.numeric(input$H_tab1_h))
                M_markerExpressionPlotInput()
                dev.off()
              })
            }
          })
          

          ##-----histogram plot-----
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
              sorted_markerNames <- colnames(v$data$expressionData)
              markerNames <- sorted_markerNames[order(sorted_markerNames)]
              initNum <- ifelse(length(markerNames) >=4, 4, 1)
              selectizeInput('m_markerSelect', 'Select Markers:', 
                             choices = markerNames, selected = markerNames[1:initNum], 
                             multiple = TRUE, width = "100%")
            }   
          })
          
          M_stackDensityPlotInput <- function(){
            m_markerSelect <- isolate(input$m_markerSelect)
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
                
                if(input$m_stackFactor == "sample"){
                  stackFactorColours <- NULL
                }else{
                  clusterMethod <- input$m_stackFactor
                  clusterVec <- v$data$clusterRes[[clusterMethod]]
                  cluster_num <- length(unique(clusterVec))
                  selectColors <- match(levels(as.factor(stackFactor)), levels(as.factor(clusterVec)))
                  if(!is.null(c$clusterCol[[clusterMethod]])){
                    stackFactorColours <- c$clusterCol[[clusterMethod]][selectColors]
                  }else{
                    stackFactorColours <- rainbow(cluster_num)[selectColors]
                  }
                }
                
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
                                       legend_title = input$m_stackFactor,
                                       stackFactorColours = stackFactorColours)
                incProgress(1/3)
                plot(gp)
                incProgress(1/3)
              })
            }
          }
          
          observeEvent(input$M_updateDensityPlot, {
            output$M_stackDensityPlot <- renderPlot({
              M_stackDensityPlotInput()
            }, height = 900, width = 950)
          })
          
          observeEvent(input$PDFHistogram, {
            if(!is.null(v$data)){
              withProgress(message="Downloading Stack Density Plot PDF files...", value=0, {
                print(getwd())
                dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
                filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Stack_Density_Plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                  filename1 <- paste0(pdfDir, .Platform$file.sep,
                                      "cytofkit_shinyAPP_Stack_Density_Plot_",
                                      Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                  i = i + 1;
                }
                pdf(filename1, 
                    width=as.numeric(input$H_tab1_w), 
                    height=as.numeric(input$H_tab1_h))
                M_stackDensityPlotInput()
                dev.off()
              })
            }
          })
          

          ##----- update marker names -----
          
          ## currently use 100 as a limit for marker number
          ## --- TODO: use reactiveValues to automatically retrive marker numbers --- ## 
          lapply(1:100, function(i) {
            output[[paste0('Marker_', i, "_name")]] <- renderUI({
              if(is.null(v$data)){
                return(NULL)
              }
              sorted_markerNames <- colnames(v$data$expressionData)
              markerNames <- sorted_markerNames[order(sorted_markerNames)]
              
              if (i <= length(markerNames)){
                markeri <- markerNames[i]
                textInput(inputId = paste0('marker_', i, "_name"), 
                          label = markeri, value = markeri, width = "30%", 
                          placeholder = "Type in your new name for this marker")
              }
            })
          })
          
          
          ## update cluster labels
          observeEvent(input$C_updateMarkerNames, {
            if(!is.null(v$data)){
              markerNames <- colnames(v$data$expressionData)
              newMarkerNames <- NULL
              for (i in 1:length(markerNames)){
                iName <- input[[paste0('marker_', i, '_name')]]
                newMarkerNames <- c(newMarkerNames, iName)
              }
              ## update new cluster colours
              colnames(v$data$expressionData) <- newMarkerNames
              ## jump to C_tab1
              updateTabsetPanel(session, "M_markerTabs", selected = "M_tab1")
            }
          })
          
          
          ##------------------------------Sample Panel-------------------------------
          
          ##-----cell percentage heatmap-----
          output$S_plotCluster <- renderUI({
            if(is.null(v$data) || is.null(clusterMethods())){
              return(NULL)
            }else{
              selectInput('s_plotCluster', 'Cluster Method:', choices = clusterMethods(), 
                          selected = clusterMethods()[1], width = "100%")
            }   
          })
          
          S_heatmapPlotInput <- reactive({
            if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_plotCluster))
              return(NULL)
            
            heatMap(data = v$data, 
                    clusterMethod = input$s_plotCluster, 
                    type = input$S_plotMethod, 
                    dendrogram = input$S_heatmap_dendrogram,
                    colPalette = input$S_heatmap_colorPalette,
                    selectSamples = input$samples,
                    cex_row_label= input$S_rowLabelSize, 
                    cex_col_label= input$S_colLabelSize, 
                    scaleMethod = input$S_scaleMethod)
            
            dev.copy2pdf(file = "cytofkit_shinyAPP_cells_heatmap_plot_plot.pdf",
                         width=as.numeric(input$S_tab1_w), 
                         height=as.numeric(input$S_tab1_h))
          })
          
          output$S_heatmapPlot <- renderPlot({
            S_heatmapPlotInput()
          }, height = 900, width = 950)
          
          observeEvent(input$PDFSamHeat, {
            if(!is.null(v$data)){
              withProgress(message="Downloading Sample Heatmap PDF files...", value=0, {
                print(getwd())
                dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
                filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Sample_Heatmap_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                  filename1 <- paste0(pdfDir, .Platform$file.sep,
                                      "cytofkit_shinyAPP_Sample_Heatmap_",
                                      Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                  i = i + 1;
                }
                file.copy("cytofkit_shinyAPP_cells_heatmap_plot_plot.pdf", filename1)
              })
            }
          })
          
          
          ##-----cell percentage line chart-----
          output$S_clusterMethod2 <- renderUI({
            if(is.null(v$data) || is.null(clusterMethods())){
              return(NULL)
            }else{
              selectInput('s_clusterMethod2', 'Cluster Method:', choices = clusterMethods(), 
                          selected = clusterMethods()[1], width = "100%")
            }   
          })
          
          output$S_clusterFilter <- renderUI({
            if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod2)){
              return(NULL)
            }else{
              clusterIDs <- sort(unique(v$data$clusterRes[[input$s_clusterMethod2]]))
              selectizeInput('s_clusterFilter', 'Filter Clusters:', 
                             choices = clusterIDs, selected = clusterIDs, 
                             multiple = TRUE, width = "100%")
            }   
          })
          
          S_rateChangePlotInput <- function(){
            if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod2) || is.null(input$s_clusterFilter))
              return(NULL)
            withProgress(message="Generating Rate Change Plot", value=0, {
              ## percentage stat
              data <- data.frame(sample = v$sampleInfo$cellSample,
                                 cluster = as.factor(v$data$clusterRes[[input$s_clusterMethod2]]),
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
          }
          
          output$S_rateChangePlot <- renderPlot({
            S_rateChangePlotInput()
          }, height = 500, width = 950)
          
          observeEvent(input$PDFrateChange, {
            if(!is.null(v$data)){
              withProgress(message="Downloading Rate Change Plot PDF files...", value=0, {
                print(getwd())
                dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
                filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Rate_Change_Plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                  filename1 <- paste0(pdfDir, .Platform$file.sep,
                                      "cytofkit_shinyAPP_Rate_Change_Plot_",
                                      Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                  i = i + 1;
                }
                pdf(filename1, 
                    width=as.numeric(input$H_tab1_w), 
                    height=as.numeric(input$H_tab1_h))
                S_rateChangePlotInput()
                dev.off()
              })
            }
          })

          
          # output$S_clusterTable <- renderTable({
          #     if(is.null(v$data) || is.null(clusterMethods()) || is.null(input$s_clusterMethod2)){
          #         return(NULL)
          #     }else{
          #         data <- data.frame(sample = v$sampleInfo$cellSample,
          #                            cluster = as.factor(v$data$clusterRes[[input$s_clusterMethod2]]),
          #                            counts = 1)
          #         
          #         statData1 <- aggregate(counts ~ ., data = data, sum)
          #         statData2 <- aggregate(counts ~ sample, data = data, sum)
          #         statData <- merge(statData1, statData2, by="sample", suffixes = c("InAll","InSample"))
          #         if(is.numeric(statData$cluster)) statData$cluster <- as.integer(statData$cluster)
          #         statData$counts <- as.integer(statData$countsInAll)
          #         statData$percentageInAll <- round(statData$countsInAll/nrow(data), 4)
          #         statData$percentageInSample <- round(statData$countsInAll/statData$countsInSample, 2)
          #         statData[, c("sample", "cluster", "counts", "percentageInSample", "percentageInAll")]
          #     }   
          # }) 
          
          
          ##-----Regroup samples-----
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
              
              cellID_number <- do.call(base::c, regmatches(v$sampleInfo$cellID,
                                                           gregexpr("_[0-9]*$", v$sampleInfo$cellID, perl=TRUE)))
              
              ## update reactive object v$sampleInfo
              ## newCellID = "sampleGroup" + "_cellID" + "globalID" to avoid dumplicates
              v$sampleInfo$newCellID <- paste0(as.character(v$sampleInfo$cellSample), 
                                               "_",
                                               1:length(cellID_number))
              
              
              ## update reactive object v$data
              expressionData <- v$data$expressionData
              row.names(expressionData) <- v$sampleInfo$newCellID
              v$data$expressionData <- expressionData
              
              ## update the project name
              v$data$projectName <- paste0(v$data$projectName, "_grouped_samples")
              
              ## update v$data$progressionRes
              if(!is.null(v$data$progressionRes)){
                sampleExpressData <- v$data$progressionRes$sampleData
                row.names(sampleExpressData) <- v$sampleInfo$newCellID[match(row.names(sampleExpressData),
                                                                             v$sampleInfo$cellID)]
                v$data$progressionRes$sampleData <- sampleExpressData
              }
              
              ## jump to S_tab1
              updateTabsetPanel(session, "S_sampleTabs", selected = "S_tab1")
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
              ## jump to S_tab1
              updateTabsetPanel(session, "S_sampleTabs", selected = "S_tab1")
            }
          })
          
          
          
          ##---------------------------Progression Panel------------------------------
          
          ##-----subset relationship plot-----
          
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
          
          P_scatterPlotInput <- function(){
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
                
                clusterMethod <- p$progressionCluster
                clusterVec <- v$data$clusterRes[[clusterMethod]]
                cluster_num <- length(unique(clusterVec))
                selectColors <- match(levels(as.factor(data$cluster)), levels(as.factor(clusterVec)))
                
                if(!is.null(c$clusterCol[[clusterMethod]])){
                  clusterColor <- c$clusterCol[[clusterMethod]][selectColors]
                }else{
                  clusterColor <- rainbow(cluster_num)[selectColors]
                }
                
                gp <- cytof_clusterPlot(data = data, 
                                        xlab = input$p_xlab, 
                                        ylab = input$p_ylab, 
                                        cluster = "cluster", 
                                        sample = "sample",
                                        title = "Subset Relationship", 
                                        type = ifelse(input$P_facetPlot, 2, 1),
                                        point_size = input$P_PointSize, 
                                        addLabel = input$P_addLabel, 
                                        labelSize = input$P_LabelSize, 
                                        sampleLabel = FALSE, 
                                        labelRepel = input$P_labelRepel,
                                        fixCoord = FALSE,
                                        clusterColor = clusterColor)
                incProgress(1/3)
                plot(gp)
                incProgress(1/3)
              })
            }
          }
          
          output$P_scatterPlot <- renderPlot({
            P_scatterPlotInput()
          }, height = 900, width = 950)
          
          observeEvent(input$PDFScatter, {
            if(!is.null(v$data)){
              withProgress(message="Downloading Progression Scatterplot PDF files...", value=0, {
                print(getwd())
                dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
                filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Scatterplot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                  filename1 <- paste0(pdfDir, .Platform$file.sep,
                                      "cytofkit_shinyAPP_Scatterplot_",
                                      Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                  i = i + 1;
                }
                pdf(filename1, 
                    width=as.numeric(input$H_tab1_w), 
                    height=as.numeric(input$H_tab1_h))
                P_scatterPlotInput()
                dev.off()
              })
            }
          })
          
          
          ##-----marker expression profile-----
          
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
              sorted_markerNames <- colnames(v$data$progressionRes$sampleData)  
              markerNames <- sorted_markerNames[order(sorted_markerNames)]
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
          
          P_markerPlotInput <- function(){
            p_markerSelect <- isolate(input$p_markerSelect)
            p_clusterSelect <- isolate(input$p_clusterSelect)
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
                                             addClusterLabel = input$P_addLabel2,
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
                                            addClusterLabel = input$P_addLabel2,
                                            clusterLabelSize = input$P_LabelSize2,
                                            segmentSize = 0.5,
                                            min_expr = NULL) 
              }
              incProgress(1/3)
              plot(pp)
              incProgress(1/3)
            })
          }
          
          observeEvent(input$P_updateRegressionPlot, {
            output$P_markerPlot <- renderPlot({
              P_markerPlotInput()
            }, height = 900, width = 950)  
          })
          
          observeEvent(input$PDFmarkerPlot, {
            if(!is.null(v$data)){
              withProgress(message="Downloading Marker Plot PDF files...", value=0, {
                print(getwd())
                dir.create(paste0("cytofkit_PDF_Plots_", Sys.Date()))
                pdfDir <- paste0(getwd(), .Platform$file.sep, "cytofkit_PDF_Plots_", Sys.Date())
                filename1 <- paste0(pdfDir, .Platform$file.sep, "cytofkit_shinyAPP_Marker_Plot_", Sys.Date(), ".pdf")
                i = 0
                while(file.exists(filename1)){
                  filename1 <- paste0(pdfDir, .Platform$file.sep,
                                      "cytofkit_shinyAPP_Marker_Plot_",
                                      Sys.Date(), "_", sprintf("%03d", i + 1), ".pdf");
                  i = i + 1;
                }
                pdf(filename1, 
                    width=as.numeric(input$H_tab1_w), 
                    height=as.numeric(input$H_tab1_h))
                P_markerPlotInput()
                dev.off()
              })
            }
          })
          
          
          ##-----Run Diffusionmap-----
          
          output$P_clusterMethod <- renderUI({
            if(is.null(v$data) || is.null(clusterMethods())){
              return(NULL)
            }else{
              selectInput('p_clusterMethod', 'Cluster Method:', choices = clusterMethods(), 
                          selected = clusterMethods()[1], width = "100%")
            }   
          })
          
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
          
          
          ## result object which will be updated by P_runDiffusionmap
          observeEvent(input$P_runDiffusionmap, {
            
            if(!is.null(v$data)){
              obj <- v$data
              usedClusters <- input$p_clusterFilter
              clusterCheck <- obj$clusterRes[[input$p_clusterMethod]] %in% usedClusters
              mdata <- obj$expressionData[clusterCheck, ]
              mcluster <- obj$clusterRes[[input$p_clusterMethod]][clusterCheck]
              withProgress(message="Running Diffusionmap", value=0, {
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
              p$progressionCluster <- input$p_clusterMethod
              ## jump to P_tab1
              updateTabsetPanel(session, "P_progressionTabs", selected = "P_tab1")
            }
          })
        }
    )
    
    
    
}

