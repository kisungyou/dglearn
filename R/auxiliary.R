# BOOLEAN TYPE
# check_sqmat  : check whether input is a square matrix
# check_markov : markov matrix





# BOOLEAN TYPE ------------------------------------------------------------
#' @keywords internal
#' @noRd
check_sqmat <- function(X){
  if (!is.matrix(X)){
    return(FALSE)
  }
  if (nrow(X)!=ncol(X)){
    return(FALSE)
  }
  return(TRUE)
}
#' @keywords internal
#' @noRd
check_markov <- function(P){
  if (!check_sqmat(P)){
    # square matrix
    return(FALSE)
  } else {
    # all nonnegative
    if (any(P<0)){
      return(FALSE)
    }
    # rowsum
    if (any(abs(rowSums(P)-rep(1,base::nrow(P))) > 1e-10)){
      return(FALSE)
    }
  }
  return(TRUE)
}

