#' Embedding by PHATE
#' 
#' Potential of Heat Diffusion for Affinity-based Transition Embedding (PHATE) 
#' is a variant of Diffusion Maps proposed by \insertCite{moon_visualizing_2019;textual}{dglearn}. 
#' Given an arbitrary affinity matrix computed by \code{aff*} functions in the \pkg{dglearn} package, 
#' it computes a low-dimensional embedding as described in the paper.
#' 
#' @param A an \eqn{(m\times m)} affinity matrix of \code{"affinity"} class.
#' @param ndim an integer-valued target dimension (default: 2).
#' @param tsteps an integer-valued timescale for smoothing or \code{"auto"} to automatically detect an optimal balancing timescale (default: \code{"auto"}).
#' @param potential type of potential distance transformation; \code{"log"} or \code{"sqrt"} (default: \code{"log"}). 
#' 
#' @return a named list containing \describe{
#' \item{Y}{an \eqn{(m\times ndim)} matrix whose rows are embedded observations.}
#' \item{tsteps}{when \code{tsteps="auto"}, it returns an optimal balancing timescale. Otherwise, it is identical to the input parameter of same name.}
#' }
#' 
#' @seealso \code{\link{findsteps}}.
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
#' ## compute multiple affinities for data-driven kernels
#' aff1 = dglearn::aff05Z(iris_mat, nbdk=5)
#' aff2 = dglearn::aff05Z(iris_mat, nbdk=10)
#' aff3 = dglearn::aff05Z(iris_mat, nbdk=20)
#' 
#' ## run PHATE with automatic stepsize rule
#' ph1 = dglearn::embedPHATE(aff1)
#' ph2 = dglearn::embedPHATE(aff2)
#' ph3 = dglearn::embedPHATE(aff3)
#' 
#' ## step size information
#' mm1 = paste0("nbdk=5:stepsize=", ph1$tsteps)
#' mm2 = paste0("nbdk=10:stepsize=",ph2$tsteps)
#' mm3 = paste0("nbdk=20:stepsize=",ph3$tsteps)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(ph1$Y, col=iris_lab, cex=0.5, pch=19, asp=1, main=mm1)
#' plot(ph2$Y, col=iris_lab, cex=0.5, pch=19, asp=1, main=mm2)
#' plot(ph3$Y, col=iris_lab, cex=0.5, pch=19, asp=1, main=mm3)
#' par(opar)
#' }
#' 
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @concept embedding
#' @export
embedPHATE <- function(A, ndim=2, tsteps="auto", potential=c("log","sqrt")){
  # ---------------------------------------------------------------
  # PREP
  # A, ndim, potential
  if (!inherits(A, "affinity")){
    stop("* embedPHATE : 'A' should be an 'affinity' class object. See the documentation for more details.")
  }
  par_ndim = round(ndim)
  par_trf  = match.arg(potential)
  if ((par_ndim < 1)||(par_ndim >= base::ncol(A))){
    stop(paste0("* embedPHATE : 'ndim' should be an integer in [1,",base::ncol(A),")."))
  }
  
  # tsteps
  if (is.character(tsteps)){
    if (identical(tsteps,"auto")){
      par_auto  = TRUE
      par_tstep = NA
    } else {
      stop("* embedPHATE : 'tsteps' is invalid. Use 'auto' or numeric value.")
    }
  } else {
    par_auto  = FALSE
    par_tstep = round(tsteps)
    if ((length(par_tstep)>1)||(!is.finite(par_tstep))||(par_tstep<1)){
      stop("* embedPHATE : 'tsteps' should be a nonnegative integer.")
    }
  }
  
  
  # ---------------------------------------------------------------
  # COMPUTE
  # degree vector
  d = base::rowSums(A)
  n = base::nrow(A)
  
  # if auto (default), compute the time step
  if (par_auto){
    # eigenvalue computation
    matL = A/base::outer(base::sqrt(d),base::sqrt(d))
    eigL = base::eigen(matL, only.values = TRUE)
    
    # compute the optimal stepsize
    par_tstep = find_steps_entmse(eigL$values)
    rm(matL)
  }
  
  # compute the diffusion operator
  P1 = diag(1/d)%*%A
  Pt = dglearn::matpower(P1, round(par_tstep))
  rm(P1)
  
  # potential distance
  if (identical(par_trf,"sqrt")){
    Pt <- base::sqrt(Pt)
  } else {
    Pt <- Pt + (1e-7)
    Pt <- Pt/base::rowSums(Pt)
    Pt <- base::log(Pt)
  }
  distPt <- stats::dist(Pt)
  rm(Pt)
  
  # compute embedding
  fun_smacof = utils::getFromNamespace("hidden_mmds","maotai")
  embedding  = fun_smacof(distPt, ndim=par_ndim)
  rm(distPt)
  
  # ---------------------------------------------------------------
  # RETURN
  output = list()
  output$Y = embedding
  output$tsteps = par_tstep
  return(output)
}

# X = as.matrix(iris[,1:4])
# y = as.factor(iris[,5])
# A = aff05Z(X, nbdk=10)
# 
# Yc = embedPHATE(A)$Yc
# Ym = embedPHATE(A)$Ym
# 
# par(mfrow=c(1,2))
# plot(Yc, pch=19, col=y, main="CMDS")
# plot(Ym, pch=19, col=y, main="MMDS")