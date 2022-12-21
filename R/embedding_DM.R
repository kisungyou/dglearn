#' Embedding by Diffusion Maps
#' 
#' 
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
}