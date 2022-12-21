#' Affinity by Adaptive Kernel by Zelnik-Manor and Perona (2005)
#' 
#' \insertCite{zelnik-manor_self-tuning_2005;textual}{dglearn} proposed to define data-driven bandwidth parameters 
#' using the nearest-neighbor distances. Let \eqn{\sigma_i} be the distance from a point \eqn{x_i} to its \code{nnbd}-th 
#' nearest neighbor. Then the affinity matrix is defined as
#' \deqn{A_{ij} = \exp(-\|x_i - x_j\|^2 / \sigma_i \sigma_j)} 
#' so that it reflects density of the data. 
#' 
#' @param X an \eqn{(m\times p)} matrix of row-stacked observations of S3 \code{dist} object of \eqn{m} observations.
#' @param nbdk neighborhood size to define data-driven bandwidth parameter (default: 5).
#' @param ... extra parameters including \describe{
#' \item{alpha}{normalization constant; one of \code{c(0, 0.5, 1)} (default: \code{0}).}
#' \item{zero.diag}{a logical; set the diagonal entries as zeros (default: \code{FALSE}).}
#' }
#' 
#' @return an \eqn{(m\times m)} affinity matrix of \code{"affinity"} class.
#' 
#' @examples 
#' \donttest{
#' ## load 'iris' data
#' data("iris")
#' 
#' ## extract the numerical part
#' iris_mat = as.matrix(iris[,1:4])
#' 
#' ## compute affinity using different neighbor sizes
#' aff1 = aff05Z(iris_mat, nbdk=5)
#' aff2 = aff05Z(iris_mat, nbdk=10)
#' aff3 = aff05Z(iris_mat, nbdk=20)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(aff1, xaxt="n", yaxt="n", xlab="", ylab="", main='nbdk=5')
#' image(aff2, xaxt="n", yaxt="n", xlab="", ylab="", main='nbdk=10')
#' image(aff3, xaxt="n", yaxt="n", xlab="", ylab="", main='nbdk=20')
#' par(opar)
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept affinity
#' @export
aff05Z <- function(X, nbdk=5, ...){
  # ---------------------------------------------------------------
  # PREP
  # explicit
  par_nbdk = max(1, round(nbdk))
  if ((par_nbdk<1)||(is.infinite(par_nbdk))||(!is.numeric(par_nbdk))){
    stop("* aff05Z : 'nbdk' should be a positive integer.")
  }
  
  # implicit
  params = list(...)
  pnames = names(params)
  
  if ("alpha"%in%pnames){
    par_alpha = as.double(params$alpha)
    if (!(par_alpha %in% c(0, 0.5, 1))){
      stop("* affgaussian : 'alpha' should be one of 0, 0.5, or 1.")
    }
  } else {
    par_alpha = 0
  }
  if ("zero.diag"%in%pnames){
    par_zerodiag = as.logical(params$zero.diag)
  } else {
    par_zerodiag = FALSE
  }
  
  # ---------------------------------------------------------------
  # CASE BRANCHING BY TYPE OF X : 'dist' or 'matrix
  if (inherits(X, "dist")){
    output = aff05Z_dist(X, par_nbdk)
  } else if (is.matrix(X)){
    output = aff05Z_mat(X, par_nbdk)
  } else {
    stop("* aff05Z : 'X' should be either a 'dist' class object or a 'matrix'.")
  }
  
  # ---------------------------------------------------------------
  # POST PROCESSING
  # zero out the diagonal entries
  if (par_zerodiag){
    diag(output) = 0
  }
  # diffusion maps applies normalization : alpha={0,0.5,1}.
  if (par_alpha > 0){
    vecd = base::rowSums(output)
    if (par_alpha < 1){
      output = output/base::outer(sqrt(vecd), sqrt(vecd))
    } else {
      output = output/base::outer(vecd, vecd)
    }
  }
  
  # ---------------------------------------------------------------
  # RETURN
  class(output) <- c("matrix","array","affinity")
  return(output)
}



# auxiliary for branching -------------------------------------------------
#' @keywords internal
#' @noRd
aff05Z_dist <- function(distobj, nbdk){
  x = as.matrix(distobj)
  n = base::nrow(x)
  par_k = max(1, round(nbdk))
  
  # k-th nearest neighbor distance
  nndist = rep(0, n)
  for (i in 1:n){
    tgt = as.vector(x[i,])
    nndist[i] = tgt[order(tgt)[nbdk+1]]
  }
  
  # build kernel & affinity
  matK = array(1, c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      matK[i,j] <- matK[j,i] <- exp(-((x[i,j])^2)/(nndist[i]*nndist[j]))
    }
  }
  
  # return
  return(matK)
}
#' @keywords internal
#' @noRd
aff05Z_mat <- function(X, nbdk){
  return(aff05Z_dist(stats::dist(X), nbdk))
}
