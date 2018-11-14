
#' QR Decomposition
#'
#' Decompose a matrix into an orthogonal matrix Q and Upper triangular matrix R
#' @param mat a matrix of real numbers
#' @return list containing Q and R matrices
#' @export
#' @examples
#' QRdecompose(rbind(c(2,-2,18),c(2,1,0),c(1,2,0)))
QRdecompose <- function(mat) {

  num_cols <- ncol(mat)
  num_rows <- nrow(mat)

  R <- mat

  num_iter <- min(num_cols, num_rows)

  H <- lapply(1:num_iter, function(i){
    x <- R[i:num_rows, i]

    len_x <- length(x)

    x_ext <- matrix(c(rep(0, (i-1)),
                      x), ncol = 1)

    e <- matrix(c(rep(0, (i-1)),
                  1,
                  rep(0, (len_x - 1))), ncol = 1)
    u <- sqrt(sum(x^2))
    v <- x_ext + sign(x[1]) * u * e

    v_inner <- (t(v) %*% v)[1]

    Hi <- diag(num_rows) - 2 * (v %*% t(v)) / v_inner

    R <<- Hi %*% R

    return(Hi)
  })

  Q <- Reduce("%*%", H)
  res <- list('Q'=Q,'R'=R)

  return(res)
}


#' LQ Decomposition
#'
#' Decompose a matrix into a Lower triangular matrix L and an orthogonal matrix Q
#' @param mat a matrix of real numbers
#' @return list containing L and Q matrices
#' @export
#' @examples
#' LQdecompose(rbind(c(2,-2,18),c(2,1,0),c(1,2,0)))
LQdecompose <- function(mat){

  mat_transpose <- t(mat)

  qr.mat_transpose <- QRdecompose(mat_transpose)

  res <- list('L' = t(qr.mat_transpose$R), 'Q' = t(qr.mat_transpose$Q))

  return(res)
}


#' Solve Rational Matrix Equation
#'
#' Given a symmetrix positive definite matrix Q and a non-singular matrix L, Find symmetric positive definite solution X such that X = Q + L (X inv) L^T
#' @param Q a symmetrix postive definite matrix of real numbers
#' @param L a non-singular matrix of real numbers
#' @param num_iterations Number of iterations to run for convergence
#' @return X : solution to the equation X = Q + L (X inv) L^T
#' @export
#' @examples
#' sol.rationalmatrix.euqation(matrix(c(2,-1,-1,2), 2, 2), rbind(c(2,3),c(2,1)))
sol.rationalmatrix.euqation <- function(Q, L, num_iterations = 50){

  if(!is.matrix(Q) || !is.matrix(L)){
    stop("Q and L must be matrices")
  }

  if(any(is.complex(Q)) || any(is.complex(L))){
    stop("Q and L must be real matrices")
  }

  if(nrow(Q) != nrow(L)){
    stop("Q and L must have same number of rows")
  }

  num_rows <- nrow(Q)

  C <- chol(Q)

  Cinv <- solve(C)

  chol0 <- cbind(C, L %*% Cinv)

  lq0 <- LQdecompose(chol0)

  L_curr <- lq0$L[1:num_rows, 1:num_rows]

  sapply(1:num_iterations, function(i){
    Y.i <- L_curr
    Y.i.inv.trans <- t(solve(Y.i))
    chol.i <- cbind(C, L %*% Y.i.inv.trans)

    lq.i <- LQdecompose(chol.i)
    L_curr <<- lq.i$L[1:num_rows, 1:num_rows]
  })

  X <- L_curr %*% t(L_curr)

  return(X)

}
