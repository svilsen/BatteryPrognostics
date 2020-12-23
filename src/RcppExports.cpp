// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// SOCADEFKFilterCpp
Rcpp::List SOCADEFKFilterCpp(const arma::colvec& I, const arma::colvec& V, const arma::colvec& Temp, const arma::colvec& Time, const arma::colvec& theta_0, const arma::colvec& ocv_0, const double& SOC_0, const std::vector<arma::mat>& P, const std::vector<arma::mat>& Q, const std::vector<arma::mat>& R, const unsigned int& K, const bool& dual, const bool& trace, const unsigned int& trace_limit);
RcppExport SEXP _BatteryPrognostics_SOCADEFKFilterCpp(SEXP ISEXP, SEXP VSEXP, SEXP TempSEXP, SEXP TimeSEXP, SEXP theta_0SEXP, SEXP ocv_0SEXP, SEXP SOC_0SEXP, SEXP PSEXP, SEXP QSEXP, SEXP RSEXP, SEXP KSEXP, SEXP dualSEXP, SEXP traceSEXP, SEXP trace_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Temp(TempSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Time(TimeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta_0(theta_0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type ocv_0(ocv_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type SOC_0(SOC_0SEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const bool& >::type dual(dualSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type trace_limit(trace_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(SOCADEFKFilterCpp(I, V, Temp, Time, theta_0, ocv_0, SOC_0, P, Q, R, K, dual, trace, trace_limit));
    return rcpp_result_gen;
END_RCPP
}
// ADUKF_Cpp
Rcpp::List ADUKF_Cpp(const arma::colvec& I, const arma::colvec& V, const arma::colvec& Temp, const arma::colvec& Time, const arma::colvec& OCV_0, const arma::colvec& Theta_0, const arma::colvec& sigma_0, const arma::mat& R_X, const arma::mat& R_V, const double& SOC_0, const unsigned int& K, const bool& trace, const unsigned int& trace_limit);
RcppExport SEXP _BatteryPrognostics_ADUKF_Cpp(SEXP ISEXP, SEXP VSEXP, SEXP TempSEXP, SEXP TimeSEXP, SEXP OCV_0SEXP, SEXP Theta_0SEXP, SEXP sigma_0SEXP, SEXP R_XSEXP, SEXP R_VSEXP, SEXP SOC_0SEXP, SEXP KSEXP, SEXP traceSEXP, SEXP trace_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Temp(TempSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Time(TimeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type OCV_0(OCV_0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Theta_0(Theta_0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type sigma_0(sigma_0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R_X(R_XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R_V(R_VSEXP);
    Rcpp::traits::input_parameter< const double& >::type SOC_0(SOC_0SEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type trace_limit(trace_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(ADUKF_Cpp(I, V, Temp, Time, OCV_0, Theta_0, sigma_0, R_X, R_V, SOC_0, K, trace, trace_limit));
    return rcpp_result_gen;
END_RCPP
}
// SOCDEFKFilterCpp
Rcpp::List SOCDEFKFilterCpp(const arma::colvec& I, const arma::colvec& V, const arma::colvec& theta_0, const arma::colvec& ocv_0, const double& SOC_0, const double& C_max, const std::vector<arma::mat>& P, const std::vector<arma::mat>& Q, const std::vector<arma::mat>& R, const double& dt, const unsigned int K, const bool& trace, const unsigned int& trace_limit);
RcppExport SEXP _BatteryPrognostics_SOCDEFKFilterCpp(SEXP ISEXP, SEXP VSEXP, SEXP theta_0SEXP, SEXP ocv_0SEXP, SEXP SOC_0SEXP, SEXP C_maxSEXP, SEXP PSEXP, SEXP QSEXP, SEXP RSEXP, SEXP dtSEXP, SEXP KSEXP, SEXP traceSEXP, SEXP trace_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type theta_0(theta_0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type ocv_0(ocv_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type SOC_0(SOC_0SEXP);
    Rcpp::traits::input_parameter< const double& >::type C_max(C_maxSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const std::vector<arma::mat>& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const double& >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type trace_limit(trace_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(SOCDEFKFilterCpp(I, V, theta_0, ocv_0, SOC_0, C_max, P, Q, R, dt, K, trace, trace_limit));
    return rcpp_result_gen;
END_RCPP
}
// RCKCpp
Rcpp::List RCKCpp(const arma::colvec& I, const arma::colvec& IC, const arma::colvec& IF, const arma::mat& R0, const std::vector<arma::mat>& Rk, const std::vector<arma::mat>& Ck, const arma::mat& Cap, const arma::mat& OCV, const arma::colvec SOCList, const arma::colvec IList, const double& dt, const double& SOCStart, const bool& trace, const unsigned int& trace_limit);
RcppExport SEXP _BatteryPrognostics_RCKCpp(SEXP ISEXP, SEXP ICSEXP, SEXP IFSEXP, SEXP R0SEXP, SEXP RkSEXP, SEXP CkSEXP, SEXP CapSEXP, SEXP OCVSEXP, SEXP SOCListSEXP, SEXP IListSEXP, SEXP dtSEXP, SEXP SOCStartSEXP, SEXP traceSEXP, SEXP trace_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type IC(ICSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type IF(IFSEXP);
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
    Rcpp::traits::input_parameter< const unsigned int& >::type trace_limit(trace_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(RCKCpp(I, IC, IF, R0, Rk, Ck, Cap, OCV, SOCList, IList, dt, SOCStart, trace, trace_limit));
    return rcpp_result_gen;
END_RCPP
}
// correct_time
arma::mat correct_time(const arma::colvec& I, const arma::colvec& V, const arma::colvec& T);
RcppExport SEXP _BatteryPrognostics_correct_time(SEXP ISEXP, SEXP VSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(correct_time(I, V, T));
    return rcpp_result_gen;
END_RCPP
}
// compare_windows_raw_cpp
Rcpp::List compare_windows_raw_cpp(const std::vector<double>& I1, const std::vector<double>& I2, const std::vector<double>& V1, const std::vector<double>& V2, const unsigned int& W, const unsigned int& R, const std::vector<double>& epsilon, const double& delta, const bool& trace, const unsigned int& trace_limit);
RcppExport SEXP _BatteryPrognostics_compare_windows_raw_cpp(SEXP I1SEXP, SEXP I2SEXP, SEXP V1SEXP, SEXP V2SEXP, SEXP WSEXP, SEXP RSEXP, SEXP epsilonSEXP, SEXP deltaSEXP, SEXP traceSEXP, SEXP trace_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type I1(I1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type I2(I2SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type trace_limit(trace_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(compare_windows_raw_cpp(I1, I2, V1, V2, W, R, epsilon, delta, trace, trace_limit));
    return rcpp_result_gen;
END_RCPP
}
// compare_windows_single_cpp
Rcpp::List compare_windows_single_cpp(const std::vector<double>& I1, const std::vector<double>& I2, const std::vector<double>& V1, const std::vector<double>& V2, const std::vector<unsigned int>& RI2, const std::vector<unsigned int>& W2, const unsigned int& W, const unsigned int& R, const std::vector<double>& epsilon, const double& delta, const bool& trace, const unsigned int& trace_limit);
RcppExport SEXP _BatteryPrognostics_compare_windows_single_cpp(SEXP I1SEXP, SEXP I2SEXP, SEXP V1SEXP, SEXP V2SEXP, SEXP RI2SEXP, SEXP W2SEXP, SEXP WSEXP, SEXP RSEXP, SEXP epsilonSEXP, SEXP deltaSEXP, SEXP traceSEXP, SEXP trace_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type I1(I1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type I2(I2SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type RI2(RI2SEXP);
    Rcpp::traits::input_parameter< const std::vector<unsigned int>& >::type W2(W2SEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type R(RSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type trace_limit(trace_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(compare_windows_single_cpp(I1, I2, V1, V2, RI2, W2, W, R, epsilon, delta, trace, trace_limit));
    return rcpp_result_gen;
END_RCPP
}
// compare_windows_ga_cpp
Rcpp::List compare_windows_ga_cpp(const std::vector<double>& I1, const std::vector<double>& I2, const std::vector<double>& V1, const std::vector<double>& V2, const std::vector<double>& Temp1, const std::vector<double>& Temp2, const unsigned int& mutation_window, const double& restrict_temperature, const double& restrict_voltage, const unsigned int& Wmin, const unsigned int& Wmax, const unsigned int& Imin, const unsigned int& N_evolution, const unsigned int& N_keep, const bool& trace, const unsigned int& trace_limit, const unsigned int& seed);
RcppExport SEXP _BatteryPrognostics_compare_windows_ga_cpp(SEXP I1SEXP, SEXP I2SEXP, SEXP V1SEXP, SEXP V2SEXP, SEXP Temp1SEXP, SEXP Temp2SEXP, SEXP mutation_windowSEXP, SEXP restrict_temperatureSEXP, SEXP restrict_voltageSEXP, SEXP WminSEXP, SEXP WmaxSEXP, SEXP IminSEXP, SEXP N_evolutionSEXP, SEXP N_keepSEXP, SEXP traceSEXP, SEXP trace_limitSEXP, SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type I1(I1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type I2(I2SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type Temp1(Temp1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type Temp2(Temp2SEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type mutation_window(mutation_windowSEXP);
    Rcpp::traits::input_parameter< const double& >::type restrict_temperature(restrict_temperatureSEXP);
    Rcpp::traits::input_parameter< const double& >::type restrict_voltage(restrict_voltageSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type Wmin(WminSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type Wmax(WmaxSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type Imin(IminSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type N_evolution(N_evolutionSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type N_keep(N_keepSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type trace_limit(trace_limitSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type seed(seedSEXP);
    rcpp_result_gen = Rcpp::wrap(compare_windows_ga_cpp(I1, I2, V1, V2, Temp1, Temp2, mutation_window, restrict_temperature, restrict_voltage, Wmin, Wmax, Imin, N_evolution, N_keep, trace, trace_limit, seed));
    return rcpp_result_gen;
END_RCPP
}
// compare_windows_variable_cpp
Rcpp::List compare_windows_variable_cpp(const std::vector<double>& I1, const std::vector<double>& I2, const std::vector<double>& V1, const std::vector<double>& V2, const unsigned int& Wmin, const unsigned int& Wmax, const double& epsilon, const double& delta, const bool& trace, const unsigned int& trace_limit);
RcppExport SEXP _BatteryPrognostics_compare_windows_variable_cpp(SEXP I1SEXP, SEXP I2SEXP, SEXP V1SEXP, SEXP V2SEXP, SEXP WminSEXP, SEXP WmaxSEXP, SEXP epsilonSEXP, SEXP deltaSEXP, SEXP traceSEXP, SEXP trace_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type I1(I1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type I2(I2SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type V1(V1SEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type V2(V2SEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type Wmin(WminSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type Wmax(WmaxSEXP);
    Rcpp::traits::input_parameter< const double& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type trace_limit(trace_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(compare_windows_variable_cpp(I1, I2, V1, V2, Wmin, Wmax, epsilon, delta, trace, trace_limit));
    return rcpp_result_gen;
END_RCPP
}
// extract_resistance_cpp
Rcpp::List extract_resistance_cpp(const std::vector<double>& I, const std::vector<double>& V, const std::vector<double>& T_s, const double& epsilon, const double& Q_max, const double& eta, const double& SOC_0);
RcppExport SEXP _BatteryPrognostics_extract_resistance_cpp(SEXP ISEXP, SEXP VSEXP, SEXP T_sSEXP, SEXP epsilonSEXP, SEXP Q_maxSEXP, SEXP etaSEXP, SEXP SOC_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type T_s(T_sSEXP);
    Rcpp::traits::input_parameter< const double& >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const double& >::type Q_max(Q_maxSEXP);
    Rcpp::traits::input_parameter< const double& >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const double& >::type SOC_0(SOC_0SEXP);
    rcpp_result_gen = Rcpp::wrap(extract_resistance_cpp(I, V, T_s, epsilon, Q_max, eta, SOC_0));
    return rcpp_result_gen;
END_RCPP
}
// particle_filtering_cpp
Rcpp::List particle_filtering_cpp(const arma::colvec& V, const arma::colvec& I, const arma::colvec& Temp, const arma::colvec& Time, const arma::colvec& OCV_c_pars, const arma::colvec& OCV_d_pars, const arma::colvec& RT_pars, const arma::colvec& Q_pars, const arma::mat& L, const double& sigma, const double& SOC_0, const unsigned int& K, const unsigned int& N, const bool& trace, const unsigned int& trace_limit);
RcppExport SEXP _BatteryPrognostics_particle_filtering_cpp(SEXP VSEXP, SEXP ISEXP, SEXP TempSEXP, SEXP TimeSEXP, SEXP OCV_c_parsSEXP, SEXP OCV_d_parsSEXP, SEXP RT_parsSEXP, SEXP Q_parsSEXP, SEXP LSEXP, SEXP sigmaSEXP, SEXP SOC_0SEXP, SEXP KSEXP, SEXP NSEXP, SEXP traceSEXP, SEXP trace_limitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Temp(TempSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Time(TimeSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type OCV_c_pars(OCV_c_parsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type OCV_d_pars(OCV_d_parsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type RT_pars(RT_parsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type Q_pars(Q_parsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type L(LSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type SOC_0(SOC_0SEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const bool& >::type trace(traceSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type trace_limit(trace_limitSEXP);
    rcpp_result_gen = Rcpp::wrap(particle_filtering_cpp(V, I, Temp, Time, OCV_c_pars, OCV_d_pars, RT_pars, Q_pars, L, sigma, SOC_0, K, N, trace, trace_limit));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BatteryPrognostics_SOCADEFKFilterCpp", (DL_FUNC) &_BatteryPrognostics_SOCADEFKFilterCpp, 14},
    {"_BatteryPrognostics_ADUKF_Cpp", (DL_FUNC) &_BatteryPrognostics_ADUKF_Cpp, 13},
    {"_BatteryPrognostics_SOCDEFKFilterCpp", (DL_FUNC) &_BatteryPrognostics_SOCDEFKFilterCpp, 13},
    {"_BatteryPrognostics_RCKCpp", (DL_FUNC) &_BatteryPrognostics_RCKCpp, 14},
    {"_BatteryPrognostics_correct_time", (DL_FUNC) &_BatteryPrognostics_correct_time, 3},
    {"_BatteryPrognostics_compare_windows_raw_cpp", (DL_FUNC) &_BatteryPrognostics_compare_windows_raw_cpp, 10},
    {"_BatteryPrognostics_compare_windows_single_cpp", (DL_FUNC) &_BatteryPrognostics_compare_windows_single_cpp, 12},
    {"_BatteryPrognostics_compare_windows_ga_cpp", (DL_FUNC) &_BatteryPrognostics_compare_windows_ga_cpp, 17},
    {"_BatteryPrognostics_compare_windows_variable_cpp", (DL_FUNC) &_BatteryPrognostics_compare_windows_variable_cpp, 10},
    {"_BatteryPrognostics_extract_resistance_cpp", (DL_FUNC) &_BatteryPrognostics_extract_resistance_cpp, 7},
    {"_BatteryPrognostics_particle_filtering_cpp", (DL_FUNC) &_BatteryPrognostics_particle_filtering_cpp, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_BatteryPrognostics(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
