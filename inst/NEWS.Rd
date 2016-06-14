\name{NEWS}
\title{News for Package \pkg{cytofkit}}

\section{Changes in version 1.4.4 (2016-05-25)}{
  \subsection{MODIFICATION}{
    \itemize{
      \item use match.arg together with match to concise the codes
      \item rename function auto_lgcl to autoLgcl
      \item rename function densVM to DensVM
      \item rename function distMethod to isomap_distMethod, more informative
      \item change usage of apply to rowSums in ClusterX, much faster
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




