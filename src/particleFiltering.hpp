#ifndef particleFiltering_H
#define particleFiltering_H

#include <RcppArmadillo.h>

// Sampling from multivariate normal distribution with armadillo
arma::mat mvrnorm(const int & n, const arma::vec & mu, const arma::mat & L);
    
// The pdf of a univariate normal distribution
double dnorm(const double & y, const double & mu, const double & sigma);  

arma::colvec mean(const arma::mat & X, const unsigned int & i);

class ParticleFiltering {
private:
    arma::colvec V, Temp, Time, I, OCV_c_pars, OCV_d_pars, RT_pars, Q_pars, Ri, Ci;
    arma::mat L, X_current;
    double OCV, R0, RT, Q, latest_non_zero, sigma, dt;
    unsigned int N, T, K, trace_limit;
    bool trace;
    
    // The parameter update functions
    void f_OCV(const double & SOC_t);
    void f_RT(const double & I_t, const double & SOC_t); 
    void f_Q(const double & I_t, const double & Temp_t);
    
    //
    double soc_transform(const double & SOC_t);
    double inverse_soc_transform(const double & SOC_t);
    
    void update_soc_parameters(const double & I_t, const double & SOC_t);
    
    arma::colvec f(const double & I_t, const arma::colvec & x_k_t_1);
    double g(const double & I_t, const arma::colvec & x_k_t); 
    
public:
    arma::colvec VT_mu, VT_sd;
    std::vector<arma::colvec> X;
    
    // Constructor
    ParticleFiltering(const arma::colvec & V_, const arma::colvec & I_, const arma::colvec & Temp_, const arma::colvec & Time_,
                      const arma::colvec & OCV_c_pars_, const arma::colvec & OCV_d_pars_, 
                      const arma::colvec & RT_pars_, const arma::colvec & Q_pars_, 
                      const arma::mat & L, const double & sigma, const double & SOC_0_,
                      const unsigned int & K_, const unsigned int & N_,  
                      const bool & trace_, const unsigned int & trace_limit_);
    
    // Filtering
    void Filter();
};

#endif