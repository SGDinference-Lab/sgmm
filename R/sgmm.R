#' SGMM: Estimation and Inference via Random Scaling
#'
#' sgmm computes the SGMM estimator and the confidence interval via random scaling method.
#'
#' @param x (n x p) matrix of regressors 
#' @param y (n x 1) vector of the dependent variable
#' @param z (n x q) matrix of instruments 
#' @param gamma_0 a constant for the learning rate (default: 1)
#' @param alpha an exponent for the learning rate (default: 0.501)
#' @param bt_start starting values for the parameters (default: zero vector)
#' @param inference character specifying the inference method (default: "rs",i.e. random scaling)
#' @param weight character specifying the type of the weighting matrix (default: "2sls")
#' @param n0 sample size for the initial sample (default: 0)
#' @param n1 sample size for the s2sls sample when "gmm" option is chosen (default: 0)
#' @param Phi_start (p x q) matrix of starting values for Phi
#' @param w_start (q x q) matrix of starting values for W 
#' @param n_perm number of permutations (default: 1). If n_perm=0, there is no permutation. 
#' @param w_option character specifying the option on beta_{i-1} when computing the weight. Default is "frequent"
#' @param path_index index number for computing the whole path of estimates and variances (default: -1. No path output)
#' @param n_eq number of equations for the System GMM estimator (default: 1. a signle equation model)
#' @note all x, y, z variables should be normalized to have sample mean zero. 
#' No intercept should be included in x and z for which p <= q.
#' The learning rate has the form: (gamma_0) times k^(alpha) for k=1,2,...  
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

sgmm = function(x=x, y=y, z=x, gamma_0=1, alpha=0.501, bt_start = NULL, 
                inference="rs", weight="2sls", n0=0, n1=0, Phi_start=Phi_start, 
                #w_start=w_start, n_perm=1, w_option="frequent"){
                w_start=w_start, n_perm=1, w_option, path_index=-1, n_eq=1){
  
  x = as.matrix(x)
  z = as.matrix(z)
  
  # Get the dimension of x and the sample size: p and n
  p = ncol(as.matrix(x))
  n = length(y)
  
  beta_hat_all = {}
  V_hat_all = {}
  
  if (n_perm == 0) {
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
    if (inference == "rs"){
      if (weight=="2sls"){
        out = s2sls_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start)
      } else if (weight=="2sls_so"){
        out = s2sls_so_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start)
      } else if (weight=="gmm"){
        out = sgmm_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start)
      } else if (weight=="gmm_new"){ 
        out = sgmm_new_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option=w_option)
      } else if (weight=="gmm_so"){
        out = sgmm_so_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option=w_option, path_index=path_index)      
      } else if (weight=="gmm_sys"){
        out = sgmm_sys_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option=w_option, path_index=path_index, n_eq=n_eq)            
      }      
      beta_hat = out$beta_hat
      V_hat = out$V_hat
      beta_hat_path = out$beta_hat_path
      V_hat_path = out$V_hat_path
      
      beta_hat_all = cbind(beta_hat_all, beta_hat)
      V_hat_all = abind::abind(V_hat_all, V_hat, along=3)
    } else if (inference == "plugin"){
      if (weight=="gmm_so"){
        out = sgmm_so_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option=w_option, path_index=path_index)            
      } else if (weight=="gmm_sys"){
        out = sgmm_sys_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option=w_option, path_index=path_index, n_eq=n_eq)            
      }      
      beta_hat = out$beta_hat
      V_hat = out$V_hat
      beta_hat_path = out$beta_hat_path
      V_hat_path = out$V_hat_path
      
      beta_hat_all = cbind(beta_hat_all, beta_hat)
      V_hat_all = abind::abind(V_hat_all, V_hat, along=3)
    }
    
  } else if (n_perm > 0) {
    for (i_p in 1:n_perm){
      
      ind = sample.int(n, n)
      x = x[ind,]
      z = z[ind,]
      y = y[ind]
      
      x = as.matrix(x)
      z = as.matrix(z)
      
      
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
      
      if (weight=="2sls"){
        out = s2sls_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start)
      } else if (weight=="2sls_so"){
        out = s2sls_so_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start)
      } else if (weight=="gmm"){
        out = sgmm_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start)
      } else if (weight=="gmm_new"){ 
        out = sgmm_new_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option=w_option)
      } else if (weight=="gmm_so"){
        out = sgmm_so_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option=w_option, path_index=path_index)      
      }
      
      beta_hat = out$beta_hat
      V_hat = out$V_hat
      beta_hat_path = out$beta_hat_path
      V_hat_path = out$V_hat_path
      
      beta_hat_all = cbind(beta_hat_all, beta_hat)
      V_hat_all = abind::abind(V_hat_all, V_hat, along=3)
      
    }
  }
  return(list(coefficient=beta_hat_all, V_hat=V_hat_all, coef_path=beta_hat_path, V_hat_path=V_hat_path))
}