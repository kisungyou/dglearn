#' Affinity by Gaussian Kernel
#' 
#' The standard method for constructing the affinity matrix is to use Gaussian kernel. 
#' When two points \eqn{x_i} and \eqn{x_j} are given, the affinity is computed as 
#' \deqn{A_{ij} = \exp(-d(x_i, d_j)^2 / 2\epsilon^2)} 
#' for a bandwidth parameter \eqn{\epsilon}. When the option \code{bandwidth="auto"} is used, 
#' \eqn{\epsilon} is chosen in an ad hoc manner.
#' 
#' @param X an \eqn{(m\times p)} matrix of row-stacked observations of S3 \code{dist} object of \eqn{m} observations.
#' @param bandwidth \eqn{\epsilon} value or \code{"auto"}.
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
#' ## compute affinity using different scales
#' aff1 = affgaussian(iris_mat, bandwidth=0.1)
#' aff2 = affgaussian(iris_mat, bandwidth=1)
#' aff3 = affgaussian(iris_mat, bandwidth="auto")
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(aff1, xaxt="n", yaxt="n", xlab="", ylab="", main='bandwidth=0.1')
#' image(aff2, xaxt="n", yaxt="n", xlab="", ylab="", main='bandwidth=1')
#' image(aff3, xaxt="n", yaxt="n", xlab="", ylab="", main='bandwidth="auto"')
#' par(opar)
#' }
#' 
#' @concept affinity
#' @export
affgaussian <- function(X, bandwidth="auto"){
  # ---------------------------------------------------------------
  # PREP
  if (is.character(bandwidth)){
    if (identical(bandwidth,"auto")){
      par_auto = TRUE
      par_bandwidth = 1
    } else {
      stop("* affgaussian : 'bandwidth' is invalid. Use 'auto' or numeric option.")
    }
  } else {
    par_auto = FALSE
    par_bandwidth = as.double(bandwidth)
    if ((length(par_bandwidth)>1)||(!is.finite(par_bandwidth))||(par_bandwidth<.Machine$double.eps)){
      stop("* affgaussian : 'bandwidth' should be a nonnegative real number.")
    }
  }
  
  # ---------------------------------------------------------------
  # CASE BRANCHING BY TYPE OF X : 'dist' or 'matrix
  if (inherits(X, "dist")){
    output = affgaussian_dist(X, par_auto, par_bandwidth)
  } else if (is.matrix(X)){
    output = src_affgaussian(X, par_auto, par_bandwidth)
  } else {
    stop("* affgaussian : 'X' should be either a 'dist' class object or a 'matrix'.")
  }
  
  # ---------------------------------------------------------------
  # RETURN
  class(output) <- c("matrix","array","affinity")
  return(output)
}



# use R when 'dist' comes in ----------------------------------------------
#' @keywords internal
#' @noRd
affgaussian_dist <- function(matD, use_auto, bandwidth){
  if (use_auto){ # automatic case
    mval = as.double(stats::median(matD))
    output = exp(-((as.matrix(matD)^2)/(2*(mval^2)))) 
  } 
  else {         # otherwise
    output = exp(-((as.matrix(matD)^2)/(2*(bandwidth^2)))) 
  }
  return(output)
}