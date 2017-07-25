\name{NEWS}
\title{News for Package \pkg{cytofkit}}

\section{Changes in version 1.8.2 (2017-07-24)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item cytofkitShinyAPP function now takes RData as argument to skip reuploading RData
      \item While choosing selected markers for dimension reduction and clustering, all markers can be visualised in the shiny app
      \item Added a "reset" button to ShinyApp to clear the session and start over
      \item Added Server file select button to allow browsing server files
    }
  }
  \subsection{MODIFICATIONS}{
    \itemize{
      \item Combined both cytofkitShinyAPP functions into one
      \item ShinyApp uses actionButton for download instead of a downloadHandler
      \item ShinyApp displays what .RData is loaded into its reactive data
      \item ShinyApp lists markers in alphabetical order
      \item Marker selection done at dimensionality reduction stage rather than at raw data transformation, to allow all expression data to be visualised at later stages
      \item Updated maintainer email address
    }
  }
}

\section{Changes in version 1.8.1 (2017-04-26)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item fixed documentation warning for function cytofkitShinyAPP2
      \item updated my maintainer email address
    }
  }
}

\section{Changes in version 1.6.6 (2017-04-11)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item add fixedLogicle transformation option on the GUI, with window popup to allow specifing the w, t, m, a parameters for logicle transformation.
      \item add openShinyAPP (boolean parameter) option in cytofkit main function, which can open shinyAPP once the analysis was done and automatically load the result object into the shinyAPP for exploration.
      \item add cytofkitShinyAPP2 function which can take cytofkit analysis_results (either file name of R object) as input and automatically load to shinyAPP once launched.
    }
  }
}


\section{Changes in version 1.6.5 (2017-03-27)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item debug that FlowSOM_k doesn't work in cytofkit main function
    }
  }
}


\section{Changes in version 1.6.4 (2017-03-17)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item in the shiny server code change function call of c in do.call to base::c
    }
  }
}

\section{Changes in version 1.6.3 (2017-03-16)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item debug the FSC|SSC channel processing error
    }
  }
}


\section{Changes in version 1.6.2 (2017-03-08)}{
  \subsection{NEW FEATURES}{
    \itemize{
      \item add default linear transformation to FSC and SSC channels
      \item add support for PDF figure download on shinyAPP, update the side panel to be tab dependent
      \item add new color palatte in heatmap (greenred and spectral) and level plot (spectral)
      \item Add cluster filter in cluster plot on shinyAPP
      \item Allow multiple annotation for same cluster (specify cluster_annotation name) on shinyAPP
      \item Allow color selection for each cluster on shinyAPP
      \item Allow modification of the marker name on shinyAPP
    }
  }
}


\section{Changes in version 1.6.1 (2016-10-27)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item updated colorPalette options in cytof_colorPlot function, also the Shiny APP, added spectral.
      \item debugged rowname conflication when regroup the samples in shinyAPP, now only use global ID, discarded the local cell ID, which avoid the dumplicate rownames conflication but results in failure in saving new FCS files.
    }
  }
}


\section{Changes in version 1.4.10 (2016-10-05)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item corrected the citation title
      \item added cluster filter in rateChange line plot in ShinyAPP
      \item debugged the error of loading back exported RData file
      \item updated the saving button, make it more robust. If cannot find FCS path, then doesn't save new FCS files
      \item added progression indicator in shinyAPP
    }
  }
}


\section{Changes in version 1.4.9 (2016-09-27)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item debug when w is negative in autoLgcl transformation, convert to logicle transformation in this case
      \item modify the cytof_writeResults function to make it more robust, with several more saving options.
      \item set autoLgcl as the default transformation method.
    }
  }
  
  \subsection{NEW FEATURES}{
    \itemize{
      \item big updates on the layout of shinyAPP, functions categorized into four panels ("cluster", "marker", "sample", "progression")
      \item added "group samples" function to relabel and group samples in sample panel
      \item added subset precentage change plot in sample panel
      \item added case checking of Nan of w in autoLgcl function
      \item Added cluster filtering and cluster table on shinyAPP diffusionmap set up page
      \item Added combined view of marker expression patten on scatter plot and marker expression trend on subset progression
      \item Added stack density plot in Marker Plot
      \item Added cluster annotation (Label Clusters) in Marker Plot
      \item add none option for transformation method to support FCS files with data already transformed
      \item add "Seperate Plot by Samples" on side panel, remove option "Label Samples by Shapes"
      \item make sample filter works on subset progression panel, support plot labeld on grid plot 
      \item add cytofkit plos computational biology paper in CITATION
      }
  }
}


\section{Changes in version 1.4.8 (2016-07-08)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item Adjusted the windown width and height of getParameters_GUI() to fit long name of FCM data
      \item debugged, set full.names = TRUE in getParameters_GUI() when fcsFile == NULL
    }
  }
}


\section{Changes in version 1.4.7 (2016-07-04)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item Concised the autoLgcl function
    }
  }
}

\section{Changes in version 1.4.6 (2016-06-14)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item concised the title in vignette, tiny modification
      \item added missed halo variable in ClusterX
      \item added transformation=FALSE in calling read.FCS function to avoid unexpected transformation for flow cytometry data.
      \item changed cast from reshap to dcast from reshap2 in function cytof_writeResults
      \item added projectName, rawFCSdir, resultDir entries in analysis_results object
      \item modified parameters for cytof_writeResults, only need analysis_results object
      \item replaced parameter uniformClusterSize in cytofkit to clusterSampleSize
      
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item redesigned shiny APP, big updates on subset progression tab, support FlowSOM and Diffusion map running. Added save buttion.
      \item rewrited most of the codes in function cytof_progression
      \item added diffusionmap in cytof_progression, updated on GUI
      \item added reverseOrder option in cytof_progressionPlot function
      \item added clusterLabelSize option in cytof_progressionPlot function
      \item added segmentSize option in cytof_progressionPlot function
      \item added cluster filetering and addClusterLabel option in cytof_progressionPlot
      \item added fixCoord option to function cytof_clusterPlot
      \item added distMethod option in cytof_progression function
      \item added distance calculation options in cytof_dimReduction
      \item added tsneSeed in cytof_dimReduction for reproducible t-SNE results
      \item added cytof_clusterStat and cytof_colorPlot function in cytof_postProcess
      \item added a button on GUI to open the resultDir once the cytof_writeResults was done, cross platform.
    }
  }
}

\section{Changes in version 1.4.5 (2016-06-12)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item Solved examples runing bug on 1.4.4
    }
  }
}

\section{Changes in version 1.4.4 (2016-05-25)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item use match.arg together with match to concise the codes
      \item renamed function auto_lgcl to autoLgcl
      \item renamed function densVM to DensVM
      \item renamed function distMethod to isomap_distMethod, more informative
      \item changed usage of apply to rowSums in ClusterX, much faster
      \item ClusterX use RANN::nn2 to speedup core/halo detection
      \item apply binary search to speedup estimateDc in ClusterX, sample size is set to 5000 by default
      \item Modified the vignettes, applied BiocStyle format
      \item debugged the dumplicate name confliction in cytof_addToFCS by modify to_add_i[[addColName]] the to_add_i[[j]]
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item added function cytofAsinh for CyTOF data transformation
      \item added transformMethod on the GUI, two options include cytofAsinh and fcmLgcl
      \item added FlowSOM clustering method
      \item added visualizaiton option of runFlowSOM, DensityPlot, DotPlot, ColorBySample in shinyAPP
      \item added sampleLabel option in function cytof_clusterPlot
      \item added cex_row_label and cex_col_label options in function cytof_heatmap
    }
  }
}

\section{Changes in version 1.4.3 (2016-05-16)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item Dubugged dimReduction selection error
    }
  }
}

\section{Changes in version 1.2.3 (2016-04-04)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item updated the shiny APP, heatmap scale, plot with selection of all markers
      \item output with a copy of expression with all markers
    }
  }
}

\section{Changes in version 1.2.3 (2016-04-04)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item updated the shiny APP, heatmap scale, plot with selection of all markers
      \item output with a copy of expression with all markers
    }
  }
}

\section{Changes in version 1.2.2 (2016-03-04)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item removed comp parameter in cytofkit-function
      \item changed parameter comp to ifCompensation in cytofkit-function
      \item adjusted the order of parameters in cytofkit-funciton
    }
  }
}

\section{Changes in version 1.2.1 (2015-12-30)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item change main function cytof_tsne_densvm to cytofkit, correspondingly change cytof_tsne_densvm_GUI to cytofkit_GUI
      \item change fcs_lgcl to cytof_exprsExtract
      \item change fcs_lgcl_merge to cytof_exprsMerge
      \item change densVM_cluster to cytof_cluster
      \item change sign_auto to auto_lgcl
      \item updated the GUI 
      \item updated vignettes 
    }
  }
  \subsection{NEW FEATURES}{
    \itemize{
      \item redesigned the pipeline with a general framework
      \item added new clustering methods ClusterX and Rphenograph
      \item added progression estimation using ISOMAP
      \item added a shiny APP to access the results interactively
    }
  }
}

\section{Changes in version 1.0.0 (2015-2-10)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item renamed package name from sidap to cytofkit
    }
  }
}




