#include <RcppArmadillo.h>
#include "ADUKF.hpp"

template <typename T> int sgn(const T val) 
{
    return(T(0) < val) - (val < T(0));
} 

template <typename T> T sigmoid(const T & x) 
{
    return T(1) / (T(1) + std::exp(-x));
}

template <typename T> T inv_sigmoid(const T & p) 
{
    return std::log(p) - std::log(T(1) - p);
}

arma::mat cholesky_decomposition(arma::mat & A) 
{
    const unsigned int & n = A.n_cols;
    arma::mat L = arma::mat(n, n, arma::fill::zeros);
    
    double S;
    for (unsigned int i = 0; i < n; i++) 
    {
        S = A(i, i);
        if (i > 0) 
        {
            for (unsigned int k = 0; k < i; k++) 
            {
                S = S - L(i, k) * L(i, k);
            }
        }
        
        L(i, i) = std::pow(S, 0.5);
        if (i < (n - 1)) 
        {
            for (unsigned int j = (i + 1); j < n; j++)
            {
                S = A(j, i);
                if (i > 0) 
                {
                    for (unsigned int k = 0; k < i; k++) 
                    {
                        S = S - L(i, k) * L(j, k);
                    }
                }
                
                L(j, i) = S / L(i, i);
            }
        }
    }
    
    return L;
}

double ADUKFModel::V_OC(const double & SOC_, const double & alpha_, const double & sign_) {
    double res = OCV[0];
    double soc_power = 1.0;
    for (unsigned int m = 1; m < M; m++) 
    {
        soc_power = soc_power * SOC_;
        res += OCV[m] * soc_power;
    }
    
    res += alpha_ * sign_;
    return res;
}

double ADUKFModel::R(const unsigned int & t, const double & SOC_,
                     const double & beta_0, const double & beta_1, 
                     const double & beta_2, const double & beta_3) 
{
    double log_RT = beta_0 + beta_1 * std::log(SOC_) + beta_2 - std::exp(beta_3) * std::abs(I[t]);
    return std::exp(log_RT);
}

double ADUKFModel::Q(const unsigned int & t, const double & beta_0, const double & beta_1, const double & beta_2) 
{
    const double s = beta_0 * std::exp(-std::exp(beta_1) * (std::abs(I[t]) / Temp[t])) + beta_2;
    return s;
}

arma::colvec ADUKFModel::f(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
{
    arma::colvec X_new = arma::zeros(X_.size());
    const double & Q_ = Q(t, theta_[8 * K + 5], theta_[8 * K + 6], theta_[8 * K + 7]);
    const double & SOC_ = sigmoid(X_[0]); 
    double SOC_new = SOC_ + Q_ * I[t] / 3600;
    
    if (SOC_new > (1.0 - 2e-16)) {
        SOC_new = 1.0 - 2e-16;
    }
    else if (SOC_new < 2e-16) {
        SOC_new = 2e-16;
    }
    
    X_new[0] = inv_sigmoid(SOC_new);
    
    for (unsigned int k = 0; k < K; k++) 
    {
        const double & Rk = R(t, sigmoid(X_new[0]), theta_[4 * k + 4], theta_[4 * k + 5], theta_[4 * k + 6], theta_[4 * k + 7]);
        const double & Ck = R(t, sigmoid(X_new[0]), theta_[4 * K + 4 * k + 4], theta_[4 * K + 4 * k + 5], theta_[4 * K + 4 * k + 6], theta_[4 * K + 4 * k + 7]);
        const double & omega_k = 1.0 / (Rk * Ck);
        
        const double & Vk_ = std::exp(-dt * omega_k);
        X_new[k + 1] = X_[k + 1] * Vk_ + Rk * (1.0 - Vk_) * I[t];
    }
    
    return X_new;
}

double ADUKFModel::h(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
{
    const double & SOC_ = sigmoid(X_[0]);
    const double & alpha_ = std::exp(theta_[8 * K + 4]);
    const double & R0 = R(t, SOC_, theta_[0], theta_[1], theta_[2], theta_[3]);
    const double & V_OC_ = V_OC(SOC_, alpha_, sign_I);
    double V_ = V_OC_ + I[t] * R0; 
    
    for (unsigned int k = 0; k < K; k++)
    {
        V_ += X_[k + 1];
    }
    
    return V_;
}

void ADUKFModel::calculate_sigma_points(const arma::colvec & z, const arma::mat & P_z, arma::mat & Z) 
{
    const unsigned int & z_dim = z.size();
    Z = arma::mat(z_dim, 2 * N + 1);
    Z.col(0) = z;
    
    arma::mat P_z_ = P_z + 1e-8 * arma::eye(P_z.n_rows, P_z .n_cols);
    arma::mat P_z_sqrt = gamma * cholesky_decomposition(P_z_);
    for (unsigned int n = 0; n < N; n++)
    {
        Z.col(n + 1) = z + P_z_sqrt.col(n);
        Z.col(N + n + 1) = z - P_z_sqrt.col(n);
    }
}

void ADUKFModel::update_f(const unsigned int & t, arma::mat & X_, const arma::colvec & theta_) 
{
    for (unsigned int n = 0; n < (2 * N + 1); n++)
    {
        X_.col(n) = f(t, X_.col(n), theta_);
    }
}

void ADUKFModel::update_h(const unsigned int & t, const arma::mat & X_, const arma::colvec & theta_) 
{
    V_sigma = arma::mat(1, 2 * N + 1);
    for (unsigned int n = 0; n < (2 * N + 1); n++)
    {
        V_sigma.col(n) = h(t, X_.col(n), theta_);
    }
}

void ADUKFModel::calculate_sigma_mean(arma::colvec & z, const arma::mat & Z) 
{
    z = arma::zeros(z.size());
    for (unsigned int n = 0; n < (2 * N + 1); n++)
    {
        const double & W_m = W_m_X[n];
        z += W_m * Z.col(n);
    }
}

void ADUKFModel::calculate_sigma_variance(const arma::colvec & X_, arma::mat & P_X_) 
{
    P_X_ = arma::zeros(N, N);
    for (unsigned int n = 0; n < (2 * N + 1); n++)
    {
        const arma::colvec X_e = X_sigma.col(n) - X_;
        const double & W_c = W_c_X[n];
        
        P_X_ += W_c * (X_e * X_e.t());
    }
    
    P_X_ += R_X;
}

void ADUKFModel::calculate_sigma_covariance(const arma::colvec & X_, const arma::colvec & V_,
                                            arma::mat & P_X_V_, arma::mat & P_V_) 
{
    P_X_V_ = arma::zeros(N, 1);
    P_V_ = arma::zeros(1, 1);
    
    for (unsigned int n = 0; n < (2 * N + 1); n++)
    {
        const arma::colvec X_e = X_sigma.col(n) - X_;
        const arma::colvec V_e = V_sigma.col(n) - V_;
        
        const double & W_c = W_c_X[n];
        
        P_X_V_ += W_c * (X_e * V_e.t());
        P_V_ += W_c * (V_e * V_e.t());
    }
    
    P_V_ = 0.5 * P_V_ + 0.5 * P_V_.t();
    P_V_ = P_V_ + 1e-8 * arma::eye(P_V_.n_rows, P_V_.n_cols);
    P_V_ += R_V;
}


ADUKFModel::ADUKFModel(const arma::colvec & I_, const arma::colvec & V_, 
                       const arma::colvec & Temp_, const arma::colvec & Time_, 
                       const arma::colvec & OCV_, const arma::colvec & Theta_,
                       const arma::colvec & sigma_,
                       const arma::mat & R_X_, const arma::mat & R_V_, 
                       const double & SOC_, const unsigned int & K_,
                       const bool & trace_, const unsigned int & trace_limit_) :
    I(I_), V(V_), Temp(Temp_), Time(Time_), OCV(OCV_), Theta(Theta_), 
    R_X(R_X_), R_V(R_V_), K(K_), trace(trace_), trace_limit(trace_limit_)
{
    T = I.size();
    N = K + 1;
    S = Theta.size();
    M = OCV_.size();
    
    SOC = arma::zeros(T);
    X = arma::zeros(N);
    V_t = arma::zeros(1);
    
    X[0] = inv_sigmoid(SOC_);
    SOC[0] = SOC_;
    
    Vhat = arma::zeros(T);
    
    alpha = sigma_[0];
    beta = sigma_[1];
    kappa = sigma_[2];
    lambda = alpha * alpha * (N + kappa);
    gamma = std::pow(N + lambda, 0.5);
    
    const double & scaling_factor = 1.0 / (2.0 * (N + lambda));
    W_m_X = scaling_factor * arma::colvec(2 * N + 1, arma::fill::ones);
    W_c_X = scaling_factor * arma::colvec(2 * N + 1, arma::fill::ones);
    
    W_m_X[0] = lambda / (N + lambda);
    W_c_X[0] = W_m_X[0] + (1 - alpha * alpha + beta);
    
    P_X = R_X; 
    P_X_V = arma::ones(N, 1); 
    P_V = R_V;
}

void ADUKFModel::Filter() 
{
    Vhat[0] = h(0, X, Theta);
    for (unsigned int t = 1; t < T; t++)
    {
        if (trace & ((t == 0) | (((t + 1) % trace_limit) == 0) | (t == (T - 1))))
        {
            Rcpp::Rcout << "Iteration: " << t + 1 << " / " << T << "\n";
        }
        
        if (std::abs(I[t]) > 0.001)
        {
            sign_I = sgn(I[t]);
        }
        
        dt = 0.0;
        if (t > 0)
        {
            dt = Time[t] - Time[t - 1];
        }
        
        //// Time update...
        calculate_sigma_points(X, P_X, X_sigma);
        update_f(t, X_sigma, Theta);
        
        calculate_sigma_mean(X, X_sigma);
        calculate_sigma_variance(X, P_X);
        
        calculate_sigma_points(X, P_X, X_sigma);
        
        update_h(t, X_sigma, Theta);
        calculate_sigma_mean(V_t, V_sigma);
        
        // Measurement update...
        calculate_sigma_covariance(X, V_t, P_X_V, P_V);
        const arma::mat & inverse_P_V = arma::inv(P_V);
        
        arma::mat G = P_X_V * inverse_P_V;
        
        X = X + G * (V[t] - V_t[0]);
        P_X = P_X - G * P_V * G.t();
        
        // Trace...
        Vhat[t] = h(t, X, Theta);
        SOC[t] = sigmoid(X[0]);
        
        if (trace & ((t == 0) | (((t + 1) % trace_limit) == 0) | (t == (T - 1))))
        {
            Rcpp::Rcout << "\tV (simulated): " << Vhat[t] << "\n"
                        << "\tV (measured): " << V[t] << "\n"
                        << "\tE: " << Vhat[t] - V[t] << "\n";
        }
    }
}

//[[Rcpp::export]] 
Rcpp::List ADUKF_Cpp(const arma::colvec & I, const arma::colvec & V, const arma::colvec & Temp, const arma::colvec & Time,
                     const arma::colvec & OCV_0, const arma::colvec & Theta_0,
                     const arma::colvec & sigma_0,
                     const arma::mat & R_X, const arma::mat & R_V, const double & SOC_0,
                     const unsigned int & K, const bool & trace, const unsigned int & trace_limit) 
{
    ADUKFModel ADUKF(I, V, Temp, Time, OCV_0, Theta_0, sigma_0, R_X, R_V, SOC_0, K, trace, trace_limit);
    ADUKF.Filter();
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = ADUKF.Vhat,
                              Rcpp::Named("SOC") = ADUKF.SOC,
                              Rcpp::Named("Theta") = ADUKF.Theta);
}

