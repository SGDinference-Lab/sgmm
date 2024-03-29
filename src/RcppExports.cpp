// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// s2sls_cpp
List s2sls_cpp(const arma::mat& x, const arma::colvec& y, const arma::mat& z, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const std::string inference, const int& n0, const arma::mat& Phi_start, const arma::mat& w_start);
RcppExport SEXP _sgmm_s2sls_cpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP inferenceSEXP, SEXP n0SEXP, SEXP Phi_startSEXP, SEXP w_startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type inference(inferenceSEXP);
    Rcpp::traits::input_parameter< const int& >::type n0(n0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Phi_start(Phi_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_start(w_startSEXP);
    rcpp_result_gen = Rcpp::wrap(s2sls_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start));
    return rcpp_result_gen;
END_RCPP
}
// s2sls_so_cpp
List s2sls_so_cpp(const arma::mat& x, const arma::colvec& y, const arma::mat& z, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const std::string inference, const int& n0, const arma::mat& Phi_start, const arma::mat& w_start);
RcppExport SEXP _sgmm_s2sls_so_cpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP inferenceSEXP, SEXP n0SEXP, SEXP Phi_startSEXP, SEXP w_startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type inference(inferenceSEXP);
    Rcpp::traits::input_parameter< const int& >::type n0(n0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Phi_start(Phi_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_start(w_startSEXP);
    rcpp_result_gen = Rcpp::wrap(s2sls_so_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start));
    return rcpp_result_gen;
END_RCPP
}
// sgmm_cpp
List sgmm_cpp(const arma::mat& x, const arma::colvec& y, const arma::mat& z, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const std::string inference, const int& n0, const arma::mat& Phi_start, const arma::mat& w_start);
RcppExport SEXP _sgmm_sgmm_cpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP inferenceSEXP, SEXP n0SEXP, SEXP Phi_startSEXP, SEXP w_startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type inference(inferenceSEXP);
    Rcpp::traits::input_parameter< const int& >::type n0(n0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Phi_start(Phi_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_start(w_startSEXP);
    rcpp_result_gen = Rcpp::wrap(sgmm_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, Phi_start, w_start));
    return rcpp_result_gen;
END_RCPP
}
// sgmm_new_cpp
List sgmm_new_cpp(const arma::mat& x, const arma::colvec& y, const arma::mat& z, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const std::string inference, const int& n0, const int& n1, const arma::mat& Phi_start, const arma::mat& w_start, const std::string w_option);
RcppExport SEXP _sgmm_sgmm_new_cpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP inferenceSEXP, SEXP n0SEXP, SEXP n1SEXP, SEXP Phi_startSEXP, SEXP w_startSEXP, SEXP w_optionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type inference(inferenceSEXP);
    Rcpp::traits::input_parameter< const int& >::type n0(n0SEXP);
    Rcpp::traits::input_parameter< const int& >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Phi_start(Phi_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_start(w_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type w_option(w_optionSEXP);
    rcpp_result_gen = Rcpp::wrap(sgmm_new_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option));
    return rcpp_result_gen;
END_RCPP
}
// sgmm_so_cpp
List sgmm_so_cpp(const arma::mat& x, const arma::colvec& y, const arma::mat& z, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const std::string inference, const int& n0, const int& n1, const arma::mat& Phi_start, const arma::mat& w_start, const std::string w_option, const int& path_index);
RcppExport SEXP _sgmm_sgmm_so_cpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP inferenceSEXP, SEXP n0SEXP, SEXP n1SEXP, SEXP Phi_startSEXP, SEXP w_startSEXP, SEXP w_optionSEXP, SEXP path_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type inference(inferenceSEXP);
    Rcpp::traits::input_parameter< const int& >::type n0(n0SEXP);
    Rcpp::traits::input_parameter< const int& >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Phi_start(Phi_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_start(w_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type w_option(w_optionSEXP);
    Rcpp::traits::input_parameter< const int& >::type path_index(path_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(sgmm_so_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option, path_index));
    return rcpp_result_gen;
END_RCPP
}
// sgmm_sys_cpp
List sgmm_sys_cpp(const arma::mat& x, const arma::colvec& y, const arma::mat& z, const double& gamma_0, const double& alpha, const arma::colvec& bt_start, const std::string inference, const int& n0, const int& n1, const arma::mat& Phi_start, const arma::mat& w_start, const std::string w_option, const int& path_index, const int& n_eq);
RcppExport SEXP _sgmm_sgmm_sys_cpp(SEXP xSEXP, SEXP ySEXP, SEXP zSEXP, SEXP gamma_0SEXP, SEXP alphaSEXP, SEXP bt_startSEXP, SEXP inferenceSEXP, SEXP n0SEXP, SEXP n1SEXP, SEXP Phi_startSEXP, SEXP w_startSEXP, SEXP w_optionSEXP, SEXP path_indexSEXP, SEXP n_eqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type z(zSEXP);
    Rcpp::traits::input_parameter< const double& >::type gamma_0(gamma_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type bt_start(bt_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type inference(inferenceSEXP);
    Rcpp::traits::input_parameter< const int& >::type n0(n0SEXP);
    Rcpp::traits::input_parameter< const int& >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Phi_start(Phi_startSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_start(w_startSEXP);
    Rcpp::traits::input_parameter< const std::string >::type w_option(w_optionSEXP);
    Rcpp::traits::input_parameter< const int& >::type path_index(path_indexSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_eq(n_eqSEXP);
    rcpp_result_gen = Rcpp::wrap(sgmm_sys_cpp(x, y, z, gamma_0, alpha, bt_start, inference, n0, n1, Phi_start, w_start, w_option, path_index, n_eq));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sgmm_s2sls_cpp", (DL_FUNC) &_sgmm_s2sls_cpp, 10},
    {"_sgmm_s2sls_so_cpp", (DL_FUNC) &_sgmm_s2sls_so_cpp, 10},
    {"_sgmm_sgmm_cpp", (DL_FUNC) &_sgmm_sgmm_cpp, 10},
    {"_sgmm_sgmm_new_cpp", (DL_FUNC) &_sgmm_sgmm_new_cpp, 12},
    {"_sgmm_sgmm_so_cpp", (DL_FUNC) &_sgmm_sgmm_so_cpp, 13},
    {"_sgmm_sgmm_sys_cpp", (DL_FUNC) &_sgmm_sgmm_sys_cpp, 14},
    {NULL, NULL, 0}
};

RcppExport void R_init_sgmm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
