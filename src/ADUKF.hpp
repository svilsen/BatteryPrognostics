#ifndef ADUKF_H
#define ADUKF_H

#include <RcppArmadillo.h>

template <typename T> int sgn(T val);

class ADUKFModel {
private:
    arma::colvec I, V, V_t, Temp, Time, OCV, X, W_m_X, W_c_X; // , W_m_Theta, W_c_Theta;
    arma::mat X_sigma, V_sigma, R_X, R_V; // , Theta_sigma;
    double alpha, beta, kappa, lambda, gamma, dt, sign_I;
    
    bool trace;
    unsigned int trace_limit, S, N, M, K, T;
    
    //
    double V_OC(const double & SOC_, const double & alpha_, const double & sign_);
    
    double R(const unsigned int & t, const double & SOC_, 
             const double & beta_0, const double & beta_1, 
             const double & beta_2, const double & beta_3);
    
    double Q(const unsigned int & t, 
             const double & beta_0, const double & beta_1, 
             const double & beta_2);
    
    //
    arma::colvec f(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_);
    double h(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_);
    
    //
    void calculate_sigma_points(const arma::colvec & z, const arma::mat & P_z, arma::mat & Z);
    void update_f(const unsigned int & t, arma::mat & X_, const arma::colvec & theta_);
    void update_h(const unsigned int & t, const arma::mat & X_, const arma::colvec & theta_);
    
    void calculate_sigma_mean(arma::colvec & z, const arma::mat & Z);
    void calculate_sigma_variance(const arma::colvec & X_, arma::mat & P_X_);
    void calculate_sigma_covariance(const arma::colvec & X_, const arma::colvec & V_,
                                    arma::mat & P_X_V_, arma::mat & P_V_);
    
public: 
    arma::colvec Theta, Vhat, SOC;
    arma::mat P_X, P_X_V, P_V;
    
    ADUKFModel(const arma::colvec & I_, const arma::colvec & V_, 
               const arma::colvec & Temp_, const arma::colvec & Time_, 
               const arma::colvec & OCV_, const arma::colvec & Theta_,
               const arma::colvec & sigma_,
               const arma::mat & R_X_, const arma::mat & R_V_, 
               const double & SOC_, const unsigned int & K_,
               const bool & trace_, const unsigned int & trace_limit_);
    
    void Filter();
};

#endif
