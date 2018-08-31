// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// RCKCpp
Rcpp::List RCKCpp(const arma::colvec& I, const arma::colvec& IC, const arma::colvec& IF, const arma::colvec& V, const arma::mat& R0, const std::vector<arma::mat>& Rk, const std::vector<arma::mat>& Ck, const arma::mat& Cap, const arma::mat& OCV, const arma::colvec SOCList, const arma::colvec IList, const double& dt, const double& SOCStart, const bool& trace, const unsigned int& traceLimit);
RcppExport SEXP _BatteryPrognostics_RCKCpp(SEXP ISEXP, SEXP ICSEXP, SEXP IFSEXP, SEXP VSEXP, SEXP R0SEXP, SEXP RkSEXP, SEXP CkSEXP, SEXP CapSEXP, SEXP OCVSEXP, SEXP SOCListSEXP, SEXP IListSEXP, SEXP dtSEXP, SEXP SOCStartSEXP, SEXP traceSEXP, SEXP traceLimitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type IC(ICSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type IF(IFSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R0(R0SEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type Rk(RkSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type Ck(CkSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Cap(CapSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type OCV(OCVSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type SOCList(SOCListSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type IList(IListSEXP);
    Rcpp::traits::input_parameter< const double& >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< const double& >::type SOCStart(SOCStartSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type traceLimit(traceLimitSEXP);
    rcpp_result_gen = Rcpp::wrap(RCKCpp(I, IC, IF, V, R0, Rk, Ck, Cap, OCV, SOCList, IList, dt, SOCStart, trace, traceLimit));
    return rcpp_result_gen;
END_RCPP
}
// SOC
Rcpp::List SOC(const arma::colvec& I, const arma::colvec& V, const bool& trace);
RcppExport SEXP _BatteryPrognostics_SOC(SEXP ISEXP, SEXP VSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(SOC(I, V, trace));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BatteryPrognostics_RCKCpp", (DL_FUNC) &_BatteryPrognostics_RCKCpp, 15},
    {"_BatteryPrognostics_SOC", (DL_FUNC) &_BatteryPrognostics_SOC, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_BatteryPrognostics(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
