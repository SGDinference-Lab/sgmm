#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List s2sls_so_cpp(const arma::mat& x, 
              const arma::colvec& y, 
              const arma::mat& z, 
              const double& gamma_0, 
              const double& alpha,
              const arma::colvec& bt_start, 
              const std::string inference,
              const int& n0,
              const arma::mat& Phi_start,
              const arma::mat& w_start){

  int n = y.n_elem;
  int p = bt_start.n_elem;
  int q = z.n_cols;

  arma::mat A_i;
  arma::vec b_i;
  double c_i = 0.0;
  arma::mat V_n;

  A_i.zeros(p,p);
  b_i.zeros(p);
  V_n.zeros(p,p);

  double A_i1 = 0.0;
  double b_i1 = 0.0;
  double V_n1 = 0.0;

  double gamma_i;
  arma::mat Phi_lag = Phi_start;
  arma::mat w_i = w_start;
  arma::colvec bt_i = bt_start;
  arma::colvec bar_bt_i;
  bar_bt_i.zeros(p);
  arma::mat i_mat = eye(q,q);
  double m_i;
  arma::vec z_i;
  arma::mat G_i;
  arma::vec H_i;
  arma::vec small_g_i;
  arma::mat so_mul;
  arma::mat so_inv;
  
  for (int obs = 1; obs < (n+1); obs++){

    gamma_i = gamma_0 * std::pow(obs, -alpha);
    z_i = trans(z.row(obs-1));
    G_i = z_i * x.row(obs-1);
    H_i = - z_i * y(obs-1);
    small_g_i = G_i * bt_i + H_i;
    so_mul = trans(Phi_lag) * w_i * Phi_lag;
    so_inv = inv_sympd(so_mul);
    //so_inv = inv(so_mul);
    bt_i = bt_i - gamma_i * so_inv * trans(Phi_lag) * w_i * small_g_i;
    Phi_lag = (n0 + obs - 1) * Phi_lag /(n0 + obs)  + (1)* G_i/(n0+obs);
    m_i = ( (n0 + obs - 1) + trans(z_i) * w_i * z_i ).eval()(0,0);
    w_i = ((n0 + obs) * w_i) / (n0 + obs - 1) * ( i_mat - z_i * trans(z_i) * w_i / m_i );
    bar_bt_i = ( bar_bt_i*(obs - 1) + bt_i ) / (obs);

    if ( inference == "rs") {
      A_i = A_i + std::pow(obs, 2.0) * bar_bt_i * trans(bar_bt_i);
      b_i = b_i + std::pow(obs, 2.0) * bar_bt_i;
      c_i = c_i + std::pow(obs, 2.0);
    } else if ( inference == "rs1") {
      A_i1 = A_i1 + std::pow(obs, 2.0) * bar_bt_i[1] * bar_bt_i[1];
      b_i1 = b_i1 + std::pow(obs, 2.0) * bar_bt_i[1];
      c_i = c_i + std::pow(obs, 2.0);
    }
  }

  if ( inference == "rs") {
    V_n = ( A_i - b_i * trans(bar_bt_i) - bar_bt_i * trans(b_i) + c_i * bar_bt_i * trans(bar_bt_i) ) / (std::pow(n, 2.0));
  } else if ( inference == "rs1") {
    V_n1 = ( A_i1 - b_i1 * bar_bt_i[1] - bar_bt_i[1] * b_i1 + c_i * bar_bt_i[1] * bar_bt_i[1] ) / (std::pow(n, 2.0));
  }
  
  
  return List::create(Named("beta_hat") = bar_bt_i,
                      Named("V_hat") = V_n);
}

