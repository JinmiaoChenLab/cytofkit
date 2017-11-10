#' RphenoGraph clustering
#' 
#' R implementation of the phenograph algorithm
#' 
#' A simple R implementation of the phenograph [PhenoGraph](http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm, 
#' which is a clustering method designed for high-dimensional single-cell data analysis. It works by creating a graph ("network") representing 
#' phenotypic similarities between cells by calculating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities 
#' using the well known [Louvain method](https://sites.google.com/site/findcommunities/) in this graph. 
#' 
#' @param data Input data matrix.
#' @param k Number of nearest neighbours, default is 30.
#'
#' @return a communities object, the operations of this class contains:
#' \item{print}{returns the communities object itself, invisibly.}
#' \item{length}{returns an integer scalar.}
#' \item{sizes}{returns a numeric vector.}
#' \item{membership}{returns a numeric vector, one number for each vertex in the graph that was the input of the community detection.}
#' \item{modularity}{returns a numeric scalar.}
#' \item{algorithm}{returns a character scalar.}
#' \item{crossing}{returns a logical vector.}
#' \item{is_hierarchical}{returns a logical scalar.}
#' \item{merges}{returns a two-column numeric matrix.}
#' \item{cut_at}{returns a numeric vector, the membership vector of the vertices.}
#' \item{as.dendrogram}{returns a dendrogram object.}
#' \item{show_trace}{returns a character vector.}
#' \item{code_len}{returns a numeric scalar for communities found with the InfoMAP method and NULL for other methods.}
#' \item{plot}{for communities objects returns NULL, invisibly.}
#'  
#' @importFrom igraph graph.data.frame cluster_louvain modularity membership
#' @importFrom Rcpp sourceCpp
#' 
#' @export
#' 
#' @author Chen Hao
#' @references Jacob H. Levine and et.al. Data-Driven Phenotypic Dissection of AML Reveals Progenitor-like Cells that Correlate with Prognosis. Cell, 2015.                
#' @examples
#' iris_unique <- unique(iris) # Remove duplicates
#' data <- as.matrix(iris_unique[,1:4])
#' Rphenograph_out <- Rphenograph(data, k = 45)
Rphenograph <- function(data, k=30){
    if(is.data.frame(data))
        data <- as.matrix(data)
    
    if(!is.matrix(data))
        stop("Wrong input data, should be a data frame or matrix!")
    
    if(k<1){
        stop("k must be a positive integer!")
    }else if (k > nrow(data)-2){
        stop("k must be smaller than the total number of points!")
    }
    
    message("Run Rphenograph starts:","\n", 
        "  -Input data of ", nrow(data)," rows and ", ncol(data), " columns","\n",
        "  -k is set to ", k)
    
    cat("  Finding nearest neighbors...")
    t1 <- system.time(neighborMatrix <- find_neighbors(data, k=k+1)[,-1])
    cat("DONE ~",t1[3],"s\n", " Compute jaccard coefficient between nearest-neighbor sets...")
    t2 <- system.time(links <- jaccard_coeff(neighborMatrix))

    cat("DONE ~",t2[3],"s\n", " Build undirected graph from the weighted links...")
    links <- links[links[,1]>0, ]
    relations <- as.data.frame(links)
    colnames(relations)<- c("from","to","weight")
    t3 <- system.time(g <- graph.data.frame(relations, directed=FALSE))
    
    # Other community detection algorithms: 
    #    cluster_walktrap, cluster_spinglass, 
    #    cluster_leading_eigen, cluster_edge_betweenness, 
    #    cluster_fast_greedy, cluster_label_prop  
    cat("DONE ~",t3[3],"s\n", " Run louvain clustering on the graph ...")
    t4 <- system.time(community <- cluster_louvain(g))
    cat("DONE ~",t4[3],"s\n")
    
    message("Run Rphenograph DONE, took a total of ", sum(c(t1[3],t2[3],t3[3],t4[3])), "s.")
    cat("  Return a community class\n  -Modularity value:", modularity(community),"\n")
    cat("  -Number of clusters:", length(unique(membership(community))))
    
    return(community)
}


#' K Nearest Neighbour Search
#'
#' Uses a kd-tree to find the p number of near neighbours for each point in an input/output dataset.
#' 
#' Use the nn2 function from the RANN package, utilizes the Approximate Near Neighbor (ANN) C++ library, 
#' which can give the exact near neighbours or (as the name suggests) approximate near neighbours 
#' to within a specified error bound. For more information on the ANN library please 
#' visit http://www.cs.umd.edu/~mount/ANN/.
#' 
#' @param data Input data matrix.
#' @param k Number of nearest neighbours.
#' 
#' @return a n-by-k matrix of neighbor indices
#' 
#' @noRd
#' 
#' @importFrom RANN nn2
find_neighbors <- function(data, k){
    nearest <- nn2(data, data, k, treetype = "kd", searchtype = "standard")
    return(nearest[[1]])
}