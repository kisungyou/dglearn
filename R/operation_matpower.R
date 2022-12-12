#' Matrix Power
#' 
#' Given a square matrix \eqn{X}, compute its matrix power
#' \deqn{X^n = X\times \cdots \times X,}
#' where the multiplication is applied \eqn{n} times. This is an equivalent to 
#' \code{matrix_power} in Python's numpy linear algebra module.
#' 
#' @param X an \eqn{(m\times m)} square matrix.
#' @param n the number of times a matrix is multiplied. It should be a finite nonnegative integer.
#' 
#' @return an \eqn{(m\times m)} square matrix for \eqn{X^n}.
#' 
#' @examples 
#' ## generate a small row-stochastic matrix
#' X0 = abs(matrix(stats::rnorm(100), ncol=10))
#' X1 = X0/base::rowSums(X0)
#' 
#' ## take a few different powers
#' X2 = matpower(X1, 2)
#' X4 = matpower(X1, 4)
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' image(X1, xaxt='n', yaxt='n', main="exponent=1")
#' image(X2, xaxt='n', yaxt='n', main="exponent=2")
#' image(X4, xaxt='n', yaxt='n', main="exponent=4")
#' par(opar)
#' 
#' @concept operation
#' @export
matpower <- function(X, n){
  # ---------------------------------------------------------------
  # PREP
  if (!check_sqmat(X)){
    stop("* matpower : input 'X' should be a square matrix.")
  }
  par_n = round(n)
  if ((!is.finite(par_n))||(par_n < 1)){
    stop("* matpower : input 'n' should be a finite nonnegative integer.")
  }

  # ---------------------------------------------------------------
  # COMPUTE AND RETURN
  # initialize
  output = X
  
  # matrix power
  if (par_n > 1){
    for (i in 1:(par_n-1)){
      output = output%*%X
    }
  }
  
  # return
  return(output)
}