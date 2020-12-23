#include <RcppArmadillo.h>
#include "ADEKF.hpp"

// Simple sign function
template <typename T> int sgn(const T val) 
{
    return(T(0) < val) - (val < T(0));
} 

// Sigmoid function and derivative
template <typename T> T sigmoid(const T & x) 
{
    const T exp_x = std::exp(-x);
    const T exp_x_1 = T(1) + std::exp(-x);
    
    return T(1) / exp_x_1;
}

template <typename T> T d_sigmoid(const T & x) 
{
    const T exp_x = std::exp(-x);
    const T exp_x_1 = T(1) + std::exp(-x);
    
    return exp_x / (exp_x_1 * exp_x_1);
}

template <typename T> T inv_sigmoid(const T & p) 
{
    return std::log(p) - std::log(T(1) - p);
}

template <typename T> T d_inv_sigmoid(const T & p) 
{
    return T(1) / p - T(1) / (p - T(1));
}

// VOC and derivative w.r.t. X
double ADEKFModel::VOC(const double & SOC_) 
{
    double alpha_ = std::exp(theta[8 * K + 4]);
    double sign_I = sgn(I_non_zero);
    
    double _SOC_ = SOC_; // sigmoid(SOC_);
    double OCV = ocv_parameters[0];
    double soc_power = 1.0;
    for (unsigned int m = 1; m < M; m++) 
    {
        soc_power = soc_power * _SOC_;
        OCV += ocv_parameters[m] * soc_power;
    }
    
    OCV += alpha_ * sign_I;
    return OCV;
}

double ADEKFModel::d_VOC(const double & SOC_) //, const double & eta_, const unsigned int & t) 
{
    double _SOC_ = SOC_; // sigmoid(SOC_);
    double OCV = ocv_parameters[1];
    double soc_power = 1.0;
    for (unsigned int m = 2; m < M; m++) 
    {
        soc_power = soc_power * _SOC_;
        OCV += m * ocv_parameters[m] * soc_power;
    }
    
    // OCV = OCV * d_sigmoid(SOC_);
    return OCV; 
}

// RT
double ADEKFModel::RT(const unsigned int & t, const double & SOC_, 
                      const double & beta_0, const double & beta_1, const double & beta_2, const double & beta_3) 
{
    double _SOC_ = SOC_; // sigmoid(SOC_);
    double log_RT = beta_0 + beta_1 * std::log(_SOC_) + beta_2 - std::exp(beta_3) * std::abs(I[t]);
    return std::exp(log_RT);
}

double ADEKFModel::Q(const unsigned int & t, const double & SOC_, 
                     const double & beta_0, const double & beta_1, const double & beta_2) 
{
    const double s = beta_0 * std::exp(-std::exp(beta_1) * (std::abs(I[t]) / Temp[t])) + beta_2;
    return s;
}

// Updating the state equation
arma::colvec ADEKFModel::f(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
{
    arma::colvec X_new = arma::zeros(X_.size());
    const double & Q_ = Q(t, X_[0], theta_[8 * K + 5], theta_[8 * K + 6], theta_[8 * K + 7]);
    const double & SOC_ = X_[0]; // sigmoid(X_[0]);
    const double & SOC_new = SOC_ + Q_ * I[t] / 3600;
    
    if (SOC_new <= 1e-6) {
        X_new[0] = 1e-6;
    }
    else if (SOC_new >= (1 - 1e-6)) {
        X_new[0] = 1.0 - 1e-6;
    }
    else {
        X_new[0] = SOC_new; // inv_sigmoid(SOC_new); 
    }
    
    for (unsigned int k = 0; k < K; k++) 
    {
        const double & Rk = RT(t, X_[0], theta_[4 * k + 4], theta_[4 * k + 5], theta_[4 * k + 6], theta_[4 * k + 7]);
        const double & Ck = RT(t, X_[0], theta_[4 * K + 4 * k + 4], theta_[4 * K + 4 * k + 5], theta_[4 * K + 4 * k + 6], theta_[4 * K + 4 * k + 7]);
        const double & omega_k = 1.0 / (Rk * Ck);
        
        const double & Vk_ = std::exp(-dt * omega_k);
        X_new[k + 1] = X_[k + 1] * Vk_ + Rk * (1.0 - Vk_) * I[t];
    }
    
    return X_new;
}

void ADEKFModel::update_F(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
{
    F = arma::zeros(K + 1, K + 1);
    
    const double & Q_ = Q(t, X_[0], theta_[8 * K + 5], theta_[8 * K + 6], theta_[8 * K + 7]);
    const double & SOC_ = X_[0]; // sigmoid(X_[0]);
    const double & SOC_new = SOC_ + Q_ * I[t] / 3600;
    
    if (SOC_new <= 1e-6) {
        F(0, 0) = 0.0;
    }
    else if (SOC_new >= (1 - 1e-6)) {
        F(0, 0) = 0.0;
    }
    else {
        F(0, 0) = 1.0; // SOC_ * SOC_ * std::exp(-X_[0]) / (SOC_new - SOC_new * SOC_new);
    }
    
    for (unsigned int k = 0; k < K; k++) 
    {
        const double & Rk = RT(t, X_[0], theta_[4 * k + 4], theta_[4 * k + 5], theta_[4 * k + 6], theta_[4 * k + 7]);
        const double & Ck = RT(t, X_[0], theta_[4 * K + 4 * k + 4], theta_[4 * K + 4 * k + 5], theta_[4 * K + 4 * k + 6], theta_[4 * K + 4 * k + 7]);
        const double & omega_k = 1.0 / (Rk * Ck);
        
        F(k + 1, k + 1) = std::exp(-dt * omega_k);
    }
}

// Updating the observation equation
double ADEKFModel::h(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
{
    const double & R0 = RT(t, X_[0], theta_[0], theta_[1], theta_[2], theta_[3]);
    double V_ = VOC(X_[0]) + I[t] * R0; 
    for (unsigned int k = 0; k < K; k++)
    {
        V_ += X_[k + 1];
    }
    
    return V_;
}

void ADEKFModel::update_H_x(const unsigned int & t, const arma::colvec & theta_) 
{
    H_x = arma::ones(1, K + 1);
    
    const double SOC_ = X[0]; // sigmoid(X[0]);
    if (SOC_ <= 1e-6) {
        H_x(0, 0) = 0.0;
    }
    else if (SOC_ >= (1 - 1e-6)) {
        H_x(0, 0) = 0.0;
    }
    else {
        H_x(0, 0) = d_VOC(X[0]); // * SOC_ * SOC_ * std::exp(-X[0]);
    }
}

void ADEKFModel::update_H_theta(const unsigned int & t) 
{
    // H_theta
    arma::mat H = arma::zeros(1, S);
    
    const double & R0 = RT(t, X[0], theta[0], theta[1], theta[2], theta[3]);
    const double & SOC_ = X[0]; // sigmoid(X[0]);
    
    H(0, 0) = I[t] * R0;
    H(0, 1) = I[t] * R0 * std::log(SOC_);
    H(0, 2) = I[t] * R0 * std::log(1.0 - SOC_);
    H(0, 3) = I[t] * R0 * (-std::abs(I[t])) * exp(theta[3]);
    H(0, 8 * K + 4) = sgn(I_non_zero) * std::exp(theta[8 * K + 4]);
    
    // X_theta
    arma::mat x_theta = arma::zeros(K + 1, S);
    
    const double & Q_ = Q(t, X[0], theta[8 * K + 5], theta[8 * K + 6], theta[8 * K + 7]);   
    const double & SOC_change = dt * I[t] / 3600;
    const double & SOC_new = SOC_ + Q_ * SOC_change;
    
    // const double & SOC_denominator = SOC_new - SOC_new * SOC_new;
    const double & Q_internal = std::exp(theta[8 * K + 5]) * (std::abs(I[t]) / Temp[t]);
    
    x_theta(0, 8 * K + 5) = (SOC_change * std::exp(-Q_internal)); // / SOC_denominator;
    x_theta(0, 8 * K + 6) = (-theta[8 * K + 5] * std::exp(-Q_internal) * Q_internal * SOC_change); // / SOC_denominator;
    x_theta(0, 8 * K + 7) = SOC_change; // / SOC_denominator;
    
    for (unsigned int k = 0; k < K; k++) 
    {
        const double & Rk = RT(t, X[0], theta[4 * k + 4], theta[4 * k + 5], theta[4 * k + 6], theta[4 * k + 7]);
        const double & Ck = RT(t, X[0], theta[4 * K + 4 * k + 4], theta[4 * K + 4 * k + 5], theta[4 * K + 4 * k + 6], theta[4 * K + 4 * k + 7]);
        const double & omega_k = 1.0 / (Rk * Ck);
        
        const double & RCk = X[k + 1] * std::exp(-dt * omega_k) * dt * omega_k - Rk * std::exp(-dt * omega_k) * dt * omega_k;
        const double & RRCk = RCk + Rk * (1.0 - std::exp(-dt * omega_k)) * I[t];
        x_theta(k + 1, 4 * k + 4) = RRCk;
        x_theta(k + 1, 4 * k + 5) = RRCk * std::log(SOC_new);
        x_theta(k + 1, 4 * k + 6) = RRCk * std::log(1.0 - SOC_new);
        x_theta(k + 1, 4 * k + 7) = RRCk * (-std::abs(I[t])) * std::exp(theta[4 * k + 7]);
        
        x_theta(k + 1, 4 * K + 4 * k + 4) = RCk;
        x_theta(k + 1, 4 * K + 4 * k + 5) = RCk * std::log(SOC_new);
        x_theta(k + 1, 4 * K + 4 * k + 6) = RCk * std::log(1.0 - SOC_new);
        x_theta(k + 1, 4 * K + 4 * k + 7) = RCk * (-std::abs(I[t])) * std::exp(theta[4 * K + 4 * k + 7]);
    }
    
    arma::mat H_ = H_x * x_theta;
    H_theta = H + H_;
}

ADEKFModel::ADEKFModel(const arma::colvec & I_, const arma::colvec & V_, const arma::colvec & Temp_, const arma::colvec & Time_, 
                       const arma::colvec & theta_, const arma::colvec & ocv_parameters_, 
                       const std::vector<arma::mat> & P_, const std::vector<arma::mat> & Q_, const std::vector<arma::mat> & R_,
                       const unsigned int & K_, const double & SOC_intial_, const bool & dual_, 
                       const bool & trace_, const unsigned int & trace_limit_) :
    I(I_), V(V_), Temp(Temp_), Time(Time_), theta(theta_), ocv_parameters(ocv_parameters_), K(K_), 
    dual(dual_), trace(trace_), trace_limit(trace_limit_), N(1)
{
    T = I.size();
    S = theta.size();
    M = ocv_parameters.size();
    
    VT = arma::colvec(T);
    SOC = arma::colvec(T);
    theta_trace = arma::mat(S, T);
    
    X = arma::zeros(K + 1);
    X[0] = SOC_intial_;
    
    Q_x = Q_[0];
    P_x = P_[0];
    R_x = R_[0];
    
    Q_theta = Q_[1];
    P_theta = P_[1];
    R_theta = R_[1];
    
    C = arma::zeros(1, 1);
    E = arma::zeros(1);
}

void ADEKFModel::Filter() 
{
    double V_hat;
    arma::colvec X_, theta_;
    arma::mat P_x_, P_theta_, L_x, L_theta, HL_x, HL_theta, HPH_x, HPH_theta, C_;

    for (unsigned int t = 0; t < T; t++) 
    {
        if (trace & ((t == 0) | ((t % trace_limit) == 0) | (t == (T - 1))))
        {
            Rcpp::Rcout << "Iteration: " << t + 1 << " / " << T << "\n";
        }
        
        if (std::abs(I[t]) > 0.001) 
        {
            I_non_zero = I[t];
        }
        
        dt = 0.0;
        if (t > 0) 
        {
            dt = Time[t] - Time[t - 1];
        }
        
        update_F(t, X, theta);
        
        //// Time-update
        // X
        X_ = f(t, X, theta);
        P_x_ = F * P_x * F.t() + Q_x;
        
        //// Measurement-update
        // X
        update_H_x(t, theta);
        
        S_x = H_x * P_x_ * H_x.t() + R_x;
        L_x = P_x_ * H_x.t() * arma::inv(S_x);
        
        V_hat = h(t, X_, theta);
        X = X_ + L_x * (V[t] - V_hat);
        if (X[0] <= 1e-6) {
            X[0] = 1e-6;
        }
        else if (X[0] >= (1.0 - 1e-6)) {
            X[0] = 1.0 - 1e-6;
        }
        
        HL_x = arma::eye(K + 1, K + 1) - L_x * H_x;
        P_x = HL_x * P_x_ * HL_x.t() + L_x * R_x * L_x.t();
        
        // Theta
        if (dual)
        {
            //// Time-update
            theta_ = theta;
            P_theta_ = P_theta + Q_theta;
            
            //// Measurement-update
            update_H_theta(t);
            
            S_theta = H_theta * P_theta_ * H_theta.t() + R_theta;
            
            L_theta = P_theta_ * H_theta.t() * arma::inv(S_theta);
            
            theta = theta_ + L_theta * (V[t] - V_hat);
            
            HL_theta = arma::eye(S, S) - L_theta * H_theta;
            P_theta = HL_theta * P_theta_ * HL_theta.t() + L_theta * R_theta * L_theta.t();
        }
        
        VT[t] = h(t, X, theta);
        
        // Trace output
        SOC[t] = X[0];
        theta_trace.col(t) = theta;
        
        E[0] = VT[t] - V[t];
        C += E * E.t() / T;
    }
    
    //// Update the Q and R covariance matrices
    // X
    HPH_x = H_x * P_x * H_x.t();
    if (HPH_x(0, 0) > C(0, 0)) {
        HPH_x(0, 0) = C(0, 0) - 1e-6;
    }
    
    R_x = C - HPH_x;
    Q_x = L_x * C * L_x.t();
    
    // Theta
    if (dual) 
    {
        HPH_theta = H_theta * P_theta * H_theta.t();
        if (HPH_theta(0, 0) > C(0, 0)) {
            HPH_theta(0, 0) = C(0, 0) - 1e-6;
        }
        
        R_theta = C - HPH_theta;
        Q_theta = L_theta * C * L_theta.t();
    }
}


//[[Rcpp::export]] 
Rcpp::List SOCADEFKFilterCpp(const arma::colvec & I, const arma::colvec & V, const arma::colvec & Temp, const arma::colvec & Time,  
                             const arma::colvec & theta_0, const arma::colvec & ocv_0, const double & SOC_0,
                             const std::vector<arma::mat> & P, const std::vector<arma::mat> & Q, const std::vector<arma::mat> & R,
                             const unsigned int & K, const bool & dual, const bool & trace, const unsigned int & trace_limit) 
{
    ADEKFModel ADEKF(I, V, Temp, Time, theta_0, ocv_0, P, Q, R, K, SOC_0, dual, trace, trace_limit);
    ADEKF.Filter();
    
    std::vector<arma::mat> Q_(2), P_(2), R_(2), S_(2);
    Q_[0] = ADEKF.Q_x;
    Q_[1] = ADEKF.Q_theta;
    
    P_[0] = ADEKF.P_x;
    P_[1] = ADEKF.P_theta;
    
    R_[0] = ADEKF.R_x;
    R_[1] = ADEKF.R_theta;
    
    S_[0] = ADEKF.S_x;
    S_[0] = ADEKF.S_theta;
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = ADEKF.VT,
                              Rcpp::Named("SOC") = ADEKF.SOC,
                              Rcpp::Named("Theta") = ADEKF.theta, 
                              Rcpp::Named("ThetaTrace") = ADEKF.theta_trace, 
                              Rcpp::Named("Q") = Q_,
                              Rcpp::Named("P") = P_,
                              Rcpp::Named("R") = R_, 
                              Rcpp::Named("S") = S_);
}

