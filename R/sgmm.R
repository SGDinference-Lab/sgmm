#' SGMM and Inference via Random Scaling
#'
#' sgmm computes the SGMM estimator and the confidence interval via random scaling method.
#'
#' @param x numeric. (n x p) matrix of regressors. Should not include 1 (the intercept)
#' @param y numeric
#' @param z numeric. (n x q) matrix of instruments. Should not include 1 (the intercept)
#' @param gamma_0 numeric
#' @param alpha numeric
#' @param bt_start numeric
#' @param inference character specifying the inference method. Default is "rs" (random scaling)
#' @param n0 numeric
#' @param Phi_start (p x q) matrix of starting values for Phi
#' @param w_start (q x q) matrix of starting values for W 
#'
#' @return
#' #' An object of class \code{"sgdi"}, which is a list containing the following
#' components:
#'
#' @export
#'
#' @examples
#' n = 1e05
#' p = 5
#' bt0 = rep(5,p)
#' x = matrix(rnorm(n*(p-1)), n, (p-1))
#' y = cbind(1,x) %*% bt0 + rnorm(n)

sgmm = function(x, y, z, gamma_0=1, alpha=0.501, bt_start = NULL,
                  inference = "rs", n0 = 0,
                  Phi_start,
                  w_start
  ){
  x = as.matrix(x)
  z = as.matrix(z)

  # Get the dimension of x and the sample size: p and n
  p = ncol(as.matrix(x))
  n = length(y)

  # Initialize the bt_t, A_t, b_t, c_t
  if (is.null(bt_start)){
    bt_t = bar_bt_t = bt_start = matrix(0, nrow=p, ncol=1)
  } else {
    bt_t = bar_bt_t = matrix(bt_start, nrow=p, ncol=1)
  }
  A_t = matrix(0, p, p)
  b_t = matrix(0, p, 1)
  c_t = 0
  V_t = NULL

  #----------------------------------------------
  # Linear Instrumental Variable Mean Regression
  #----------------------------------------------
  out = sgmm_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start)
  beta_hat = out$beta_hat
  V_hat = out$V_hat

  return(list(coefficient=beta_hat, V_hat = V_hat))

}
