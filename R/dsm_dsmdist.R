#' Diffusion-based Similarity Metric between Multiple Stochastic Matrices
#' 
#' Given multiple row-stochastic matrices with nonnegative entries, 
#' compute all pairwise diffusion-based similarity metric (DSM) distances.
#' 
#' @param Plist a length-\eqn{N} list of \eqn{(m\times m)} row-stochastic matrices.
#' 
#' @return an \eqn{(N\times N)} symmetric matrix of pairwise DSM values.
#' 
#' @examples 
#' \donttest{
#' ## load 'iris' and 'USArrests' data
#' data("iris")
#' data("USArrests")
#' 
#' ## extract the numerical part and select 25 random objects with scaling
#' dat_i = as.matrix(scale(iris[sample(1:150, 25, replace=F),1:4]))
#' dat_u = as.matrix(scale(USArrests[sample(1:50, 25, replace=F),]))
#' 
#' ## generate 10 replicates of each with perturbation and create kernels
#' list20 = vector("list", length=20L)
#' sdlev  = 0.3
#' for (i in 1:10){
#'   perturb_i   = dat_i + matrix(rnorm(25*4, sd=sdlev), ncol=4)
#'   aff_i       = exp(-as.matrix(stats::dist(perturb_i))^2)
#'   list20[[i]] = aff_i/base::rowSums(aff_i)
#' }
#' for (u in 11:20){
#'   perturb_u   = dat_u + matrix(rnorm(25*4, sd=sdlev), ncol=4)
#'   aff_u       = exp(-as.matrix(stats::dist(perturb_u))^2)
#'   list20[[u]] = aff_u/base::rowSums(aff_u)
#' }
#' 
#' ## compute the distance
#' dist20 = dsmdist(list20)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(pty="s")
#' image(x=1:20, y=1:20, z=dist20, zlim=c(0,1), 
#'       xaxt="n", yaxt="n", xlab="", ylab="", 
#'       main="All DSMs for 2 Classes")
#' par(opar)
#' }
#' 
#' @concept dsm
#' @export
dsmdist <- function(Plist){
  # ---------------------------------------------------------------
  # PREP
  if (!is.list(Plist)){
    stop("* dsmdist : 'Plist' should be a list.")
  }
  if (any(unlist(lapply(Plist, check_markov))==FALSE)){
    stop("* dsmdist : at least one of the input matrices is not a row-stochastic matrix.")
  }
  if (length(unique(unlist(lapply(Plist, nrow))))!=1){
    stop("* dsmdist : at least one of the input matrices has a non-matching size.")
  }
  
  # ---------------------------------------------------------------
  # COMPUTE AND RETURN
  # total number of Markov kernels
  N = length(Plist)
  
  # iterate
  output = array(0,c(N,N))
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      output[i,j] <- output[j,i] <- src_dsm_dist2diff(Plist[[i]], Plist[[j]])
    }
  }
  
  # return
  return(output)
}
