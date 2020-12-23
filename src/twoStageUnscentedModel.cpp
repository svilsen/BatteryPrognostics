#include <RcppArmadillo.h>

class TwoStageUnscentedModel {
private:
    arma::colvec I, V, Rk, Omegak;
    arma::mat R, Q, V_tilde, theta_tilde;
    
    double dt, eta, C_max, gamma, xi, kappa, lambda, epsilon, R0, V_OC;
    
    bool trace;
    unsigned int trace_limit, S, N, K, T;
    
    arma::colvec f(const unsigned int & t, const arma::colvec & theta_) 
    {
        arma::colvec theta_new = theta_;
        theta_new[0] = theta_[0] + dt * (I[t] / (3600 * C_max));
        for (unsigned int k = 0; k < K; k++) 
        {
            const double & V_k_ = theta_[k + 1];
            const double & omega_k = Omegak[k]; // = std::exp(-std::abs(theta_[K + 2 + k]));
            const double & r_k = Rk[k]; // std::exp(theta_[2 * K + 2 + k]);
            
            const double & Vk_ = std::exp(-dt * omega_k);
            theta_new[1 + k] = V_k_ * Vk_ + r_k * (1.0 - Vk_) * I[t];
        }
        
        return theta_new;
    }
    
    double h(const unsigned int & t, const arma::colvec & theta_) 
    {
        double V_ = V_OC + R0 * I[t];
        for (unsigned int k = 0; k < K; k++)
        {
            V_ += theta_[k + 1];
        }
        
        return V_; 
    }
    
    void calculate_sigma_points(const unsigned int & t, const arma::colvec & theta_, const arma::mat & P_theta_) 
    {
        unsigned int S_ = theta_.size();
        
        theta_tilde = arma::mat(S_, 2 * S_ + 1);
        theta_tilde.col(0) = theta_;
        
        arma::mat P_theta_sqroot = std::pow(S_ + lambda, 0.5) * arma::sqrtmat_sympd(P_theta_); 
        for (unsigned int s = 0; s < S_; s++) 
        {
            theta_tilde.col(s + 1) = theta_ + P_theta_sqroot.col(s);
            theta_tilde.col(S_ + s + 1) = theta_ - P_theta_sqroot.col(s);
        }
    }
    
    void update_theta_tilde(const unsigned int & t) 
    {
        const unsigned int & S_ = theta_tilde.n_cols;
        for (unsigned int s = 0; s < S_; s++) 
        {
            theta_tilde.col(s) = f(t, theta_tilde.col(s));
        }
    }
    
    void update_V_tilde(const unsigned int & t) 
    {
        const unsigned int & S_ = theta_tilde.n_cols;
        
        V_tilde = arma::mat(1, S_);
        for (unsigned int s = 0; s < S_; s++) 
        {
            V_tilde.col(s) = h(t, theta_tilde.col(s));
        }
    }
    
    void update_theta_V_filter(arma::colvec & theta_, arma::colvec & V_) 
    {
        const unsigned int & S_ = theta_.size();
        lambda = (S_ + kappa) * gamma * gamma  - S_;
        
        theta_ = arma::zeros(S_);
        V_ = arma::zeros(1);
        for (unsigned int s = 0; s < (2 * S_ + 1); s++) 
        {
            double W_m;
            if (s == 0) 
            {
                W_m = lambda / (S_ + lambda);
            }
            else 
            {
                W_m = 1.0 / (2 * (S_ + lambda));
            }
            
            theta_ += W_m * theta_tilde.col(s);
            V_ += W_m * V_tilde.col(s);
        }
    }
    
    void update_covariance_matrices_filter(arma::mat & P_theta_, arma::mat & P_theta_y, arma::mat & P_y, 
                                           const arma::colvec & theta_, const arma::colvec & VT_) 
    {
        const unsigned int & S_ = theta_.size();
        lambda = (S_ + kappa) * gamma * gamma  - S_;
        
        P_theta_ = arma::zeros(S_, S_);
        P_theta_y = arma::zeros(S_, 1);
        P_y = arma::zeros(1, 1);
        
        for (unsigned int s = 0; s < (2 * S_ + 1); s++) 
        {
            double W_c;
            if (s == 0) 
            {
                W_c = lambda / (S_ + lambda) + (1.0 - gamma * gamma + xi);
            }
            else 
            {
                W_c = 1.0 / (2 * (S_ + lambda));
            }
            
            
            const arma::colvec theta_error = theta_tilde.col(s) - theta_;
            const arma::colvec VT_error = V_tilde.col(s) - VT_;
            
            P_theta_ += W_c * (theta_error * theta_error.t());
            P_theta_y += W_c * (theta_error * VT_error.t());
            P_y += W_c * (VT_error * VT_error.t());
        }
        
        P_y = 0.5 * P_y + 0.5 * P_y;
        P_y = P_y + 1e-6 * arma::eye(P_y.n_rows, P_y.n_cols);
    }
    
    void update_internal_parameters(const arma::colvec & theta_) 
    {
        V_OC = 2.0 + std::abs(theta_[K + 1]);
        
        for (unsigned int k = 0; k < K; k++) 
        {
            Rk[k] = std::exp(theta_[K + 2 + k]);
            Omegak[k] = std::exp(-(theta_[2 * K + 2 + k] * theta_[2 * K + 2 + k]));
        }
        
        R0 = std::exp(theta_[3 * K + 2]);
    } 
    
    void rescale_theta_(arma::colvec & theta_, arma::mat & P_theta_, unsigned int S_) 
    {
        unsigned int theta_size = theta_.size();
        if (theta_size != S_) 
        {
            theta_.resize(S_);
            P_theta_.resize(S_, S_);
            
            if (S_ == S) 
            {
                theta_[K + 1] = V_OC - 2.0;
                
                for (unsigned int k = 0; k < K; k++) 
                {
                    theta_[K + 2 + k] = std::log(Rk[k]);
                    theta_[2 * K + 2 + k] = std::pow(-std::log(Omegak[k]), 0.5);
                }
                
                theta_[3 * K + 2] = std::log(R0);
            }
        }
    }
    
    void update_theta_trace(const unsigned int & t) 
    {
        theta_trace(0, t) = theta[0];
        
        theta_trace(K + 1, t) = V_OC;
        for (unsigned int k = 0; k < K; k++) 
        {
            theta_trace(k + 1, t) = theta[k + 1];
            theta_trace(K + 2 + k, t) = Rk[k];
            theta_trace(2 * K + 2 + k, t) = Omegak[k];
        }
        
        theta_trace(3 * K + 2, t) = R0;
    }
    
public :
    arma::colvec theta, SOC, Vhat;
    arma::mat P_theta, theta_trace;
    
    TwoStageUnscentedModel(const arma::colvec & I_, const arma::colvec & V_, const arma::colvec & theta_, 
                const arma::mat & P_theta_, const double & dt_, const unsigned int K_,
                const double SOC_0_, const double C_max_, const double & eta_,
                const double & gamma_, const double & xi_, const double & kappa_, const double & epsilon_,
                const arma::mat & R_, const arma::mat & Q_, const bool & trace_, const unsigned int & trace_limit_) :
        I(I_), V(V_), theta(theta_), dt(dt_), K(K_), 
        trace(trace_), trace_limit(trace_limit_), C_max(C_max_), eta(eta_), 
        gamma(gamma_), xi(xi_), kappa(kappa_), epsilon(epsilon_),
        R(R_), Q(Q_), N(1)
    {
        T = I.size();
        S = theta.size();
        
        // SOC = arma::colvec(T);
        // SOC[0] = SOC_0_;
        
        Vhat = arma::colvec(T);
        theta_trace = arma::mat(S, T);
        
        Rk = arma::zeros(K);
        Omegak = arma::zeros(K);
        
        update_internal_parameters(theta);
    }
    
    void Filter() 
    {
        arma::colvec theta_, V_;
        arma::mat P_theta_, P_theta_y, P_y;
        
        Vhat[0] = h(0, theta);
        
        double percentage_error = std::abs((Vhat[0] - V[0]) / V[0]);
        for (unsigned int t = 1; t < T; t++) 
        {
            if (trace & ((t == 0) | (((t + 1) % trace_limit) == 0) | (t == (T - 1))))
            {
                Rcpp::Rcout << "Iteration: " << t + 1 << " / " << T << "\n";
            }
            
            if (percentage_error < epsilon) 
            {
                theta_ = theta;
                P_theta_ = P_theta;
                rescale_theta_(theta_, P_theta_, K + 1);
                
                //
                calculate_sigma_points(t, theta_, P_theta_);
                
                //// Time update
                update_theta_tilde(t);
                update_V_tilde(t);
                
                update_theta_V_filter(theta_, V_);
                update_covariance_matrices_filter(P_theta_, P_theta_y, P_y, theta_, V_);
                
            }
            else 
            {
                rescale_theta_(theta_, P_theta, S);
                update_internal_parameters(theta);
            }
            
            percentage_error = std::abs((Vhat[t] - V[t]) / V[t]);
            update_theta_trace(t);
         }
    }
};

