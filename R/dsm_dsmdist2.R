#' Diffusion-based Similarity Metric between Two Stochastic Matrices
#' 
#' Given two \eqn{(m\times m)} row-stochastic matrices with nonnegative entries, 
#' compute the diffusion-based similarity metric (DSM). 
#' 
#' @param P1 an \eqn{(m\times m)} row-stochastic matrix.
#' @param P2 an \eqn{(m\times m)} row-stochastic matrix.
#' 
#' @return the DSM value. 
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
#' ## simple affinity and Markov kernels
#' aff_i = exp(-as.matrix(stats::dist(dat_i))^2)
#' aff_u = exp(-as.matrix(stats::dist(dat_u))^2)
#' 
#' ker_i = aff_i/base::rowSums(aff_i)
#' ker_u = aff_u/base::rowSums(aff_u)
#' 
#' ## compute the distance
#' print(paste0("DSM value is ",round(dsmdist2(ker_i, ker_u),4)))
#' }
#' 
#' @concept dsm
#' @export
dsmdist2 <- function(P1, P2){
  # ---------------------------------------------------------------
  # PREP
  if (!check_markov(P1)){
    stop("* dsmdist2 : 'P1' is not a proper row-stochastic matrix.")
  }
  if (!check_markov(P2)){
    stop("* dsmdist2 : 'P2' is not a proper row-stochastic matrix.")
  }
  if (base::nrow(P1)!=base::nrow(P2)){
    stop("* dsmdist2 : 'P1' and 'P2' are not of same size.")
  }
  
  # ---------------------------------------------------------------
  # COMPUTE AND RETURN
  return(as.double(src_dsm_dist2diff(P1,P2)))
}


# direct access -----------------------------------------------------------
#' @keywords internal
#' @noRd
nocheck_dsmdist2 <- function(P1, P2){
  return(as.double(src_dsm_dist2diff(P1,P2)))
}
