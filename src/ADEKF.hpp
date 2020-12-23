#ifndef ADEFK_H
#define ADEFK_H

#include <RcppArmadillo.h>

// Simple sign function
template <typename T> int sgn(T val);

// Sigmoid function and derivative
template <typename T> T sigmoid(const T & x);
template <typename T> T d_sigmoid(const T & x);
template <typename T> T inv_sigmoid(const T & p);
template <typename T> T d_inv_sigmoid(const T & p);


class ADEKFModel {
private:
    //
    arma::colvec I, V, Temp, Time, ocv_parameters, E;
    double dt, I_non_zero;  
    bool trace, dual;
    unsigned int trace_limit, S, N, K, T, M;
    
    arma::mat F, H_x, H_theta, C;
    
    // VOC and derivative
    double VOC(const double & SOC_);
    double d_VOC(const double & SOC_);
    
    // 
    double RT(const unsigned int & t, const double & SOC_, 
              const double & beta_0, const double & beta_1, 
              const double & beta_2, const double & beta_3);
    double Q(const unsigned int & t, const double & SOC_, 
             const double & beta_0, const double & beta_1, 
             const double & beta_2);
    
    //
    arma::colvec f(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_);
    void update_F(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_);
    
    //
    double h(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_);
    void update_H_x(const unsigned int & t, const arma::colvec & theta_);
    void update_H_theta(const unsigned int & t);
    
public:
    //
    arma::colvec theta, X, VT, SOC;
    arma::mat Q_x, Q_theta, P_x, P_theta, R_x, R_theta, S_x, S_theta, theta_trace;

    //
    ADEKFModel(const arma::colvec & I_, const arma::colvec & V_, 
               const arma::colvec & Temp_, const arma::colvec & Time_, 
               const arma::colvec & theta_, const arma::colvec & ocv_parameters_, 
               const std::vector<arma::mat> & P_, const std::vector<arma::mat> & Q_, const std::vector<arma::mat> & R_,
               const unsigned int & K_, const double & SOC_intial_, const bool & dual_,
               const bool & trace_, const unsigned int & trace_limit_);
    
    // 
    void Filter();
};

#endif
