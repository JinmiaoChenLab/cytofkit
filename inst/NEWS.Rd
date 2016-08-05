\name{NEWS}
\title{News for Package \pkg{cytofkit}}

\section{Changes in version 1.5.7 (2016-07-13)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item Smart handle of non-transformation markers(FSC,SSC), and exclude of channels(like Event, Time) when markers==NULL specified in cytof_exprsExtract function
    }
  }
  
  \subsection{NEW FEATURES}{
    \itemize{
      \item Redesigned the shinyAPP tab panels: Cluster Plot, Marker Plot, Subset Progression 
      \item Added cluster filtering and cluster table on shinyAPP diffusionmap set up page
      \item Added combined view of marker expression patten on scatter plot and marker expression trend on subset progression
      \item Added stack density plot in Marker Plot
      \item Added cluster annotation (Label Clusters) in Marker Plot
    }
}


\section{Changes in version 1.5.6 (2016-07-08)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item Adjusted the windown width and height of getParameters_GUI() to fit long name of FCM data
      \item debugged, set full.names = TRUE in getParameters_GUI() when fcsFile == NULL
    }
  }
}


\section{Changes in version 1.5.5 (2016-07-04)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item Concised the autoLgcl function
    }
  }
}

\section{Changes in version 1.5.4 (2016-06-14)}{
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
