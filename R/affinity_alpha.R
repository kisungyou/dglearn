#' Affinity by Alpha-Decay Kernel
#' 
#' \insertCite{moon_visualizing_2019}{dglearn} proposed a data-driven kernel called alpha-decay kernel 
#' that uses a user-specified decay parameter as well as the nearest-neighbor distance. 
#' Let \eqn{\sigma_i} be the distance from a point \eqn{x_i} to its \code{nnbd}-th 
#' nearest neighbor. Then the affinity matrix using the alpha-decay kernel is defined as
#' \deqn{A_{ij} = \frac{1}{2}
#' \exp\left(-\left(\frac{\|x_i - x_j\|}{\sigma_i}\right)^\alpha\right) + 
#' \frac{1}{2}\exp\left(-\left(\frac{\|x_i - x_j\|}{\sigma_j}\right)^\alpha\right).}
#' 
#' @param X an \eqn{(m\times p)} matrix of row-stacked observations of S3 \code{dist} object of \eqn{m} observations.
#' @param nbdk neighborhood size to define data-driven bandwidth parameter (default: 5).
#' @param alpha decay parameter (default: 2).
#' @param ... extra parameters including \describe{
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
#' aff1 = affalpha(iris_mat, nbdk=5)
#' aff2 = affalpha(iris_mat, nbdk=10)
#' aff3 = affalpha(iris_mat, nbdk=20)
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
affalpha <- function(X, nbdk=5, alpha=2, ...){
  # ---------------------------------------------------------------
  # PREP
  # explicit
  par_nbdk = max(1, round(nbdk))
  if ((par_nbdk<1)||(is.infinite(par_nbdk))||(!is.numeric(par_nbdk))){
    stop("* affalpha : 'nbdk' should be a positive integer.")
  }
  par_alpha = as.double(alpha)
  if ((par_alpha < .Machine$double.eps)||(is.infinite(par_alpha))){
    stop("* affalpha : 'alpha' should be a nonnegative real number.")
  }
  
  # implicit
  params = list(...)
  pnames = names(params)
  
  if ("zero.diag"%in%pnames){
    par_zerodiag = as.logical(params$zero.diag)
  } else {
    par_zerodiag = FALSE
  }
  
  # ---------------------------------------------------------------
  # CASE BRANCHING BY TYPE OF X : 'dist' or 'matrix
  if (inherits(X, "dist")){
    output = affalpha_dist(X, par_nbdk, par_alpha)
  } else if (is.matrix(X)){
    output = affalpha_mat(X, par_nbdk, par_alpha)
  } else {
    stop("* affafinity : 'X' should be either a 'dist' class object or a 'matrix'.")
  }
  
  # ---------------------------------------------------------------
  # POST PROCESSING
  # zero out the diagonal entries
  if (par_zerodiag){
    diag(output) = 0
  }
  
  # ---------------------------------------------------------------
  # RETURN
  class(output) <- c("matrix","array","affinity")
  return(output)
}


# auxiliary ---------------------------------------------------------------
#' @keywords internal
#' @noRd
affalpha_dist <- function(dx, nbdk, alpha){
  x = as.matrix(dx)
  n = base::nrow(x)
  nbdk  = max(1, round(nbdk))
  alpha = max(sqrt(.Machine$double.eps), as.double(alpha))
  
  # k-th nearest distance
  nndist = rep(0,n)
  for (i in 1:n){
    tgt = as.vector(x[i,])
    nndist[i] = tgt[order(tgt)[nbdk+1]]
  }
  
  # Build Kernel Matrix
  matK = array(1,c(n,n))
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      term1 = exp(-((x[i,j]/nndist[i])^(alpha)))
      term2 = exp(-((x[i,j]/nndist[j])^(alpha)))
      matK[i,j] <- matK[j,i] <- 0.5*(term1+term2)
    }
  }
  
  # return
  return(matK)
}
#' @keywords internal
#' @noRd
affalpha_mat <- function(x, nbdk, alpha){
  return(affalpha_dist(stats::dist(x), nbdk, alpha))
}