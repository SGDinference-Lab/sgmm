#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List sgmm_so_cpp(const arma::mat& x, 
                  const arma::colvec& y, 
                  const arma::mat& z, 
                  const double& gamma_0, 
                  const double& alpha,
                  const arma::colvec& bt_start, 
                  const std::string inference,
                  const int& n0,
                  const int& n1,
                  const arma::mat& Phi_start,
                  const arma::mat& w_start,
                  const std::string w_option,
                  const int& path_index){

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
  arma::colvec bar_bt_i = bt_start;
  arma::mat i_mat = eye(q,q);
  double m_i;
  arma::vec z_i = trans(z.row(0));
  arma::mat G_i = z_i * x.row(0);
  arma::vec H_i = - z_i * y(0);
  arma::vec small_g_i = G_i * bt_i + H_i;
  arma::vec bar_small_g_i;
  bar_small_g_i.zeros(q);
  arma::vec small_g_i_at_bt_bar;
  small_g_i_at_bt_bar.zeros(q);
  double norm_g_i;
  arma::vec bar_bt_i_fix;
  arma::mat so_mul;
  arma::mat so_inv;
  
  // Declare matrices for the path
  arma::vec bar_bt_i_path;
  bar_bt_i_path.zeros(n);
  arma::vec so_inv_path;
  so_inv_path.zeros(n);
  
  // S2SLS procedure for the first 'n1' observations. 
  for (int obs = 1; obs < (n1+1); obs++){
    
    gamma_i = gamma_0 * std::pow(obs, -alpha);
    z_i = trans(z.row(obs-1));
    G_i = z_i * x.row(obs-1);
    H_i = - z_i * y(obs-1);
    small_g_i = G_i * bt_i + H_i;
    //bar_small_g_i = ( bar_small_g_i*(obs - 1) + small_g_i ) / obs;
    so_mul = trans(Phi_lag) * w_i * Phi_lag;
    //so_inv = inv_sympd(so_mul);
    so_inv = pinv(so_mul);
    bt_i = bt_i - gamma_i * so_inv * trans(Phi_lag) * w_i * small_g_i;
    Phi_lag = (n0 + obs - 1) * Phi_lag /(n0 + obs)  + (1)* G_i/(n0+obs);
    m_i = ( (n0 + obs - 1) + trans(z_i) * w_i * z_i ).eval()(0,0);
    w_i = ((n0 + obs) * w_i) / (n0 + obs - 1) * ( i_mat - z_i * trans(z_i) * w_i / m_i );
    bar_bt_i = ( bar_bt_i*(obs - 1) + bt_i ) / (obs);

    if ( inference == "rs" || inference == "plugin") {
      A_i = A_i + std::pow(obs, 2.0) * bar_bt_i * trans(bar_bt_i);
      b_i = b_i + std::pow(obs, 2.0) * bar_bt_i;
      c_i = c_i + std::pow(obs, 2.0);
    }
    
    if ( path_index >= 0 ) {
      bar_bt_i_path(obs-1) = bar_bt_i(path_index-1);
      so_inv_path(obs-1) = so_inv(path_index-1,path_index-1);
    }
  }
  
  // bar_bt_i_fix = {1.0,1.0,1.0,1.0,1.0};
  bar_bt_i_fix = bar_bt_i;
  
  // SGMM procedure for the remaining 'n-n1' observations. 
  for (int obs = (n1+1); obs < (n+1); obs++){

    gamma_i = gamma_0 * std::pow(obs, -alpha);
    z_i = trans(z.row(obs-1));
    G_i = z_i * x.row(obs-1);
    H_i = - z_i * y(obs-1);
    small_g_i = G_i * bt_i + H_i;
    bar_small_g_i = ( bar_small_g_i*(obs - 1) + small_g_i ) / obs;
    so_mul = trans(Phi_lag) * w_i * Phi_lag;
    //so_inv = inv_sympd(so_mul);
    so_inv = pinv(so_mul);
    bt_i = bt_i - gamma_i * so_inv * trans(Phi_lag) * w_i * small_g_i;
    Phi_lag = (n0 + obs - 1) * Phi_lag /(n0 + obs)  + (1)* G_i/(n0+obs);
    
    if (w_option == "frequent"){
      m_i = ( (n0 + obs - 1) + trans(small_g_i) * w_i * small_g_i ).eval()(0,0);
      w_i = ((n0 + obs) * w_i) / (n0 + obs - 1) * ( i_mat - (small_g_i * trans(small_g_i) ) * w_i / m_i );
    } else if (w_option == "single"){
      small_g_i = G_i * bar_bt_i_fix + H_i;
      m_i = ( (n0 + obs - 1) + trans(small_g_i) * w_i * small_g_i ).eval()(0,0);
      w_i = ((n0 + obs) * w_i) / (n0 + obs - 1) * ( i_mat - (small_g_i * trans(small_g_i) ) * w_i / m_i );
    } else if (w_option == "infrequent"){
      small_g_i = G_i * bar_bt_i_fix + H_i;
      m_i = ( (n0 + obs - 1) + trans(small_g_i) * w_i * small_g_i ).eval()(0,0);
      w_i = ((n0 + obs) * w_i) / (n0 + obs - 1) * ( i_mat - (small_g_i * trans(small_g_i) ) * w_i / m_i );
      if (obs*1.0/n == 0.4) {
       bar_bt_i_fix = bar_bt_i;
      }
      if (obs*1.0/n == 0.6) {
       bar_bt_i_fix = bar_bt_i;
      }
      if (obs*1.0/n == 0.8) {
       bar_bt_i_fix = bar_bt_i;
      }
    } else {
      Rcerr << "w_option is chosen incorrectly. \n";
    }
    
    // Alternative 1: Demean bar_small_g_i. We apply the Woodbury formular twice after getting w_i
    // w_i = w_i + (w_i * bar_small_g_i * trans(bar_small_g_i) * w_i) / (1-trans(bar_small_g_i)*w_i*bar_small_g_i).eval()(0,0);
    
    // Alternative 2: weight is evaluated at bar_bt_i
    // small_g_i_at_bt_bar = G_i * bar_bt_i + H_i;
    // m_i = ( (n0 + obs - 1) + trans(small_g_i_at_bt_bar) * w_i * small_g_i_at_bt_bar ).eval()(0,0);    
    // w_i = ((n0 + obs) * w_i) / (n0 + obs - 1) * ( i_mat - (small_g_i_at_bt_bar * trans(small_g_i_at_bt_bar) ) * w_i / m_i );    
    
    // Alternative 3: normalize the small_g_i
    // norm_g_i = norm(small_g_i);
    // small_g_i = small_g_i/norm_g_i;
    // m_i = ( (n0 + obs - 1) + trans(small_g_i) * w_i * small_g_i ).eval()(0,0);
    // w_i = ((n0 + obs) * w_i) / (n0 + obs - 1) * ( i_mat - (small_g_i * trans(small_g_i) ) * w_i / m_i );
    
    // Alternative 4: bar_bt_i fix and compute weight
    // small_g_i = G_i * bar_bt_i_fix + H_i;
    // m_i = ( (n0 + obs - 1) + trans(small_g_i) * w_i * small_g_i ).eval()(0,0);
    // w_i = ((n0 + obs) * w_i) / (n0 + obs - 1) * ( i_mat - (small_g_i * trans(small_g_i) ) * w_i / m_i );
    // if (obs*10/n == 4) {
    //   bar_bt_i_fix = bar_bt_i;
    // }
    // if (obs*10/n == 6) {
    //   bar_bt_i_fix = bar_bt_i;
    // }
    // if (obs*10/n == 8) {
    //   bar_bt_i_fix = bar_bt_i;
    // }

    // Alternative 5: bar_bt_i fix at the true beta_0
    // small_g_i = G_i * bar_bt_i_fix + H_i;
    // m_i = ( (n0 + obs - 1) + trans(small_g_i) * w_i * small_g_i ).eval()(0,0);
    // w_i = ((n0 + obs) * w_i) / (n0 + obs - 1) * ( i_mat - (small_g_i * trans(small_g_i) ) * w_i / m_i );
        
    
    bar_bt_i = ( bar_bt_i*(obs - 1) + bt_i ) / (obs);

    if ( inference == "rs" ) {
      A_i = A_i + std::pow(obs, 2.0) * bar_bt_i * trans(bar_bt_i);
      b_i = b_i + std::pow(obs, 2.0) * bar_bt_i;
      c_i = c_i + std::pow(obs, 2.0);
    } else if ( inference == "rs1") {
      A_i1 = A_i1 + std::pow(obs, 2.0) * bar_bt_i[1] * bar_bt_i[1];
      b_i1 = b_i1 + std::pow(obs, 2.0) * bar_bt_i[1];
      c_i = c_i + std::pow(obs, 2.0);
    }
    if ( path_index >= 0 ) {
      bar_bt_i_path(obs-1) = bar_bt_i(path_index-1);
      so_inv_path(obs-1) = so_inv(path_index-1,path_index-1);
    }
    
  }

  if ( inference == "rs") {
    V_n = ( A_i - b_i * trans(bar_bt_i) - bar_bt_i * trans(b_i) + c_i * bar_bt_i * trans(bar_bt_i) ) / (std::pow(n, 2.0));
  } else if ( inference == "plugin") {
    V_n = so_inv;
  } else if ( inference == "rs1") {
    V_n1 = ( A_i1 - b_i1 * bar_bt_i[1] - bar_bt_i[1] * b_i1 + c_i * bar_bt_i[1] * bar_bt_i[1] ) / (std::pow(n, 2.0));
  }
  
  return List::create(Named("beta_hat") = bar_bt_i,
                      Named("V_hat") = V_n,
                      Named("beta_hat_path") = bar_bt_i_path,
                      Named("V_hat_path") = so_inv_path);
}

