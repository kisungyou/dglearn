#' Embedding by Diffusion Maps
#' 
#' Diffusion Maps \insertCite{coifman_geometric_2005,coifman_diffusion_2006}{dglearn} is a nonlinear manifold learning technique based on an idea 
#' to approximate the diffusion process on data manifold. Unlike the full implementation 
#' available in other packages like \pkg{Rdimtools} or \pkg{diffusionMap}, this code 
#' uses a precomputed affinity matrix by one of \code{aff*} functions in the \pkg{dglearn} package. 
#' 
#' @param A an \eqn{(m\times m)} affinity matrix of \code{"affinity"} class.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param tsteps an integer-valued timescale for smoothing (default: 1).
#' 
#' @return a named list containing \describe{
#' \item{Y}{an \eqn{(m\times ndim)} matrix whose rows are embedded observations.}
#' }
#' 
#' @examples 
#' \donttest{
#' ## load 'iris' data
#' data("iris")
#' 
#' ## extract the numerical part and label
#' iris_mat = as.matrix(iris[,1:4])
#' iris_lab = as.factor(iris[,5])
#' 
#' ## compute affinity with automatic Gaussian kernel
#' iris_aff = dglearn::affgaussian(iris_mat)
#' 
#' ## run diffusion maps at different step sizes
#' dm1 = embedDM(iris_aff, tsteps=1)$Y
#' dm2 = embedDM(iris_aff, tsteps=2)$Y
#' dm3 = embedDM(iris_aff, tsteps=5)$Y
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(dm1, col=iris_lab, cex=0.5, pch=19, asp=1, main="tsteps=1")
#' plot(dm2, col=iris_lab, cex=0.5, pch=19, asp=1, main="tsteps=2")
#' plot(dm3, col=iris_lab, cex=0.5, pch=19, asp=1, main="tsteps=5")
#' par(opar)
#' }
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept embedding
#' @export
embedDM <- function(A, ndim=2, tsteps=1){
  # ---------------------------------------------------------------
  # PREP
  if (!inherits(A, "affinity")){
    stop("* embedDM : 'A' should be an 'affinity' class object. See the documentation for more details.")
  }
  par_ndim  = round(ndim)
  par_tstep = round(tsteps)
  
  if ((par_ndim < 1)||(par_ndim >= base::ncol(A))){
    stop(paste0("* embedDM : 'ndim' should be an integer in [1,",base::ncol(A),")."))
  }
  if ((length(par_tstep)>1)||(!is.finite(par_tstep))||(par_tstep<1)){
    stop("* embedDM : 'tsteps' should be a nonnegative integer.")
  }
  
  # ---------------------------------------------------------------
  # COMPUTE
  # transition matrix
  P = A/base::rowSums(A)
  
  # eigendecomposition
  eigP = RSpectra::eigs(P, k = round(par_ndim+1), which="LM")
  rm(P)

  # compute an embedding
  embedding = base::Re(eigP$vectors[,2:(par_ndim+1)])%*%diag(base::Re(eigP$values[2:(par_ndim+1)])^par_tstep)

  # ---------------------------------------------------------------
  # RETURN
  output = list()
  output$Y = embedding
  return(output)
}