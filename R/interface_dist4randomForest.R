#' Compute Distance for Query Data and a 'randomForest' Object
#' 
#' Given an object of S3 class \code{randomForest} and a query data for prediction, 
#' compute the distance among samples in a query induced by the given random forest. 
#' Currently, no \code{NA} value is allowed in a query dataset. 
#' 
#' @param forest a S3 class \code{randomForest} object.
#' @param query a \code{data.frame} of \eqn{m} observations whose distance will be computed.
#' @param dtype name of the distance; either \code{"proximity"} or \code{"geodesic"} (default: \code{"proximity"}).
#' 
#' @return an object of \code{\link{dist}} class for \eqn{m} observations in a query.
#' 
#' @examples 
#' \donttest{
#' ## load the 'iris' data for classification
#' data("iris")
#' iris.df = as.data.frame(iris)
#' 
#' ## set seed for replicability
#' set.seed(496)
#' 
#' ## fit the randomForest classification model with default settings
#' iris.rf = randomForest::randomForest(Species ~ ., data=iris.df)
#' 
#' ## compute two types of distances
#' dP = dist4randomForest(iris.rf, iris.df, dtype="proximity")
#' dG = dist4randomForest(iris.rf, iris.df, dtype="geodesic")
#' 
#' ## visualize two distance matrices
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(as.matrix(dP), xaxt='n', yaxt='n', main="dtype='proximity'")
#' image(as.matrix(dG), xaxt='n', yaxt='n', main="dtype='geodesic'")
#' par(opar)
#' }
#' 
#' 
#' @concept interface
#' @export
dist4randomForest <- function(forest, query, dtype=c("proximity","geodesic")){
  # ---------------------------------------------------------------
  # PREP
  # input : forest
  if (!inherits(forest, "randomForest")){
    stop("* dist4randomForest : 'forest' is not a valid object.")
  }
  # input : query
  if (!is.data.frame(query)){
    stop("* dist4randomForest : 'query' should be a dataframe.")
  }
  if (any(is.na(query))){
    stop("* dist4randomForest : 'query' does not allow NA/NaN values.")
  }
  # input : dtype
  mydtype = match.arg(dtype)
  
  # ---------------------------------------------------------------
  # COMPUTE
  # case branching by dtype
  if (identical(mydtype, "proximity")){
    # extract the distance by (1.0-proximity)
    mat_dist = 1.0-as.matrix(predict(forest, query, proximity=TRUE)$proximity)
    rownames(mat_dist) = NULL
    colnames(mat_dist) = NULL
    
    # return a distance object
    return(stats::as.dist(mat_dist))
  } else {
    # extract the distance by geodesic average
    # compute the node matrix
    nodemat = attr(predict(forest, query, nodes=TRUE), "nodes")
    
    # dimensions
    nobj  = base::nrow(nodemat)
    ntree = base::ncol(nodemat) 
    
    # iterate over all trees
    gdmat = array(0, c(nobj, nobj))
    for (it in 1:ntree){
      # extract a tree & membership
      tree_it = randomForest::getTree(forest, k=it, labelVar = TRUE)
      node_it = as.vector(nodemat[,it])
      
      # compute the geodesic distance & update
      gdmat = gdmat + (rftreeSingleGeodesicTree(tree_it, node_it)/ntree)
    }
    
    # return a distance object
    return(stats::as.dist(gdmat))
  }
}



# auxiliary for 'dist4randomForest' ---------------------------------------
#' @keywords internal
#' @noRd
rftreeSingleGeodesicTree <- function(tree_now, node_now){
  # extract some relevant information
  nnode = base::nrow(tree_now)
  vec_l = round(tree_now$`left daughter`)
  vec_r = round(tree_now$`right daughter`)
  
  # create an adjacency matrix
  mat_adj = array(0, c(nnode, nnode))
  for (i in 1:length(vec_l)){
    id_l = vec_l[i]
    id_r = vec_r[i]
    
    if (id_l > 0){
      mat_adj[i,id_l] <- mat_adj[id_l,i] <- 1
      mat_adj[i,id_r] <- mat_adj[id_r,i] <- 1 
    }
  }
  
  # igraph manipulation + distance computation
  obj_igraph = igraph::graph_from_adjacency_matrix(mat_adj, 
                                                   mode="undirected", 
                                                   diag=FALSE, 
                                                   weighted=NULL)
  obj_distmat = igraph::distances(obj_igraph)
  
  # now arrange into returnable object
  nleafs = length(node_now)
  output = array(0,c(nleafs, nleafs))
  for (i in 1:(nleafs-1)){
    id_i = node_now[i]
    for (j in (i+1):nleafs){
      id_j = node_now[j]
      
      output[i,j] <- output[j,i] <- obj_distmat[id_i, id_j]
    }
  }
  
  # return
  return(output)
}

