#' Estimate the Diffusion Step Size
#' 
#' The core idea of diffusion maps is to power the Markov kernel matrix 
#' a specified number of times in order to see multi-scale aspect of the 
#' data manifold. While the step size is a job left for users, there exist 
#' a number of methods to \emph{estimate} the desired step sizes for a diffusion 
#' operator to take.
#' 
#' @param P an \eqn{(m\times m)} row-stochastic matrix.
#' @param method name of the method to be used, including\describe{
#' \item{\code{"entmse"}}{changepoint of von Neumann entropy by two linear regression lines (default).}
#' }
#' 
#' @return the estimated step size for \code{P}.
#' 
#' @seealso \code{\link{matpower}}.
#' 
#' @examples 
#' \donttest{
#' ## load 'iris' data
#' data("iris")
#' 
#' ## scale the numerical part
#' i_dat = as.matrix(scale(iris[,1:4]))
#' 
#' ## simple affinity and Markov kernels
#' i_aff = exp(-as.matrix(stats::dist(i_dat))^2)
#' i_ker = i_aff/base::rowSums(i_aff)
#' 
#' ## compute the optimal step size
#' t_entmse = dglearn::findsteps(i_ker, method="entmse")
#' 
#' ## power the transition matrices at different scales
#' Pa = dglearn::matpower(i_ker, 1)
#' Pb = dglearn::matpower(i_ker, 5)
#' Pc = dglearn::matpower(i_ker, 10)
#' Pentmse = dglearn::matpower(i_ker, t_entmse)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3), pty="s")
#' image(Pa, xaxt="n", yaxt="n", xlab="", ylab="", main="step size=1")
#' image(Pb, xaxt="n", yaxt="n", xlab="", ylab="", main="step size=5")
#' image(Pc, xaxt="n", yaxt="n", xlab="", ylab="", main="step size=10")
#' image(Pentmse, xaxt="n", yaxt="n", xlab="", ylab="", 
#'       main=paste0("entmse:",t_entmse," steps"))
#' par(opar)
#' }
#' 
#' @concept operation
#' @export
findsteps <- function(P, method=c("entmse")){
  # ---------------------------------------------------------------
  # PREP
  # explicit
  if (!check_markov(P)){
    stop("* export : 'P' is not a proper row-stochastic matrix.")
  }
  
  # implicit : if uses 'eigenvalue'-based method, add to the list!
  par_method   = match.arg(method)
  par_eigbased = c("entmse")
  
  # ---------------------------------------------------------------
  # COMPUTE
  # eigenvalue-based methods
  if (par_method %in% par_eigbased){
    # EVD
    eigP    = base::eigen(P)
    # compute
    output = switch(par_method,
                    "entmse" = find_steps_entmse(eigP$values))
  }
  
  
  # ---------------------------------------------------------------
  # RETURN
  return(output)
}

# internal use ------------------------------------------------------------
#' @keywords internal
#' @noRd
find_steps_entmse <- function(eigPval){
  # select the strictly positive ones
  eigA = base::Re(eigPval)
  eigA = eigA[(eigA>100*.Machine$double.eps)]
  
  # compute von Neumann entropy
  n_max = 1000
  vec.t = 1:n_max
  vec.H = rep(0, n_max)
  for (it in 1:n_max){
    eig.t = eigA^it
    eig.t = eig.t/base::sum(eig.t)
    # eig.t = eig.t[(eig.t > .Machine$double.eps)]
    vec.H[it] = (-base::sum(eig.t*base::log(eig.t)))
  }
  
  # select the non-nan ones
  id_min   = is.finite(vec.H)
  finite_t = vec.t[id_min]
  finite_H = vec.H[id_min]
  
  # optimal 't' by clamped least-square fits
  opt.t = aux_knee_mse(finite_t, finite_H)
  return(opt.t)
}