# BOOLEAN TYPE
# check_sqmat  : check whether input is a square matrix
# check_markov : markov matrix
#
# DOUBLE TYPE
# aux_knee_mse : least-squares knee point detection given as an index




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



# DOUBLE TYPE -------------------------------------------------------------
#' @keywords internal
#' @noRd
hidden_knee_clamped_basic <- function(x, y){
  m = length(x)
  c = x[1]
  d = y[1]
  a = x[m]
  b = y[m]
  
  y2 = (((b-d)/(a-c))*(x-c))+d
  return(sum((y-y2)^2))
}
#' @keywords internal
#' @noRd
aux_knee_mse <- function(x, y){
  x = as.vector(x)
  y = as.vector(y)
  n = length(x)
  if (n < 3){
    stop("* aux_knee_mse : length must be larger than 2.")
  }
  scores = rep(Inf, n)
  for (i in 2:(n-1)){
    x.left = x[1:i]
    y.left = y[1:i]
    
    x.right = x[i:n]
    y.right = y[i:n]
    
    term1 = hidden_knee_clamped_basic(x.left, y.left)
    term2 = hidden_knee_clamped_basic(x.right, y.right)
    scores[i] = term1+term2
  }
  return(which.min(scores)) # return the index of the minimal SSE's
}
