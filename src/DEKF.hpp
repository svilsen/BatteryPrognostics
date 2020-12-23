#ifndef DEFK_H
#define DEFK_H

#include <RcppArmadillo.h>

template <typename T> int sgn(T val);

class DEKFModel {
private:
    //
    arma::colvec I, V, ocv_parameters;
    double dt, C_max, I_non_zero;  
    bool trace;
    unsigned int trace_limit, S, N, K, T;
    
    // 
    double V_OC(const double & SOC_);
    
    double V_OC_(const double & SOC_);
    
    double h(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_);
    
    arma::colvec f(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_);
    
    void update_F(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_);
    
    void update_H_x();
    
    void update_H_theta(const unsigned int & t);
        
public:
    //
    arma::colvec theta, X, VT, SOC;
    arma::mat F, H_x, H_theta, Q_x, Q_theta, P_x, P_theta, R_x, R_theta, theta_trace;
    double LogLikelihood;
    
    //
    DEKFModel(const arma::colvec & I_, const arma::colvec & V_, const arma::colvec & theta_,
              const arma::colvec & ocv_parameters_, const std::vector<arma::mat> & P, 
              const std::vector<arma::mat> & Q, const std::vector<arma::mat> & R, 
              const double & dt_, const unsigned int K_, const double SOC_intial_, const double C_max_, 
              const bool & trace_, const unsigned int & trace_limit_);
    
    // 
    void Filter();
};

#endif