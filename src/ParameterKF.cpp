#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

class ParameterKF {
private:
    //
    arma::colvec I, V;
    arma::mat V_k, R, Q;
    double dt, eta, C_max, gamma, xi, kappa, lambda; // 
    bool trace;
    unsigned int traceLimit, S, N, K, T, lag, V_dim, M;
    
    void update_SOC(const double & t) 
    {
        const double SOC_ = SOC[t - 1] + dt * (I[t] / (3600 * C_max));
        
        if (SOC_ < 0.0)
        {
            SOC[t] = 0.0;
        }
        else if (SOC_ > 1.0) 
        {
            SOC[t] = 1.0;
        }
        
        SOC[t] = SOC_;
    }
    
    void update_V(const unsigned int & t, arma::colvec & V_, arma::mat & R_) 
    {
        if (lag >= (t + 1)) 
        {
            V_dim = t + 1;
            R_ = arma::zeros(V_dim, V_dim);
        }
        
        V_ = arma::colvec(V_dim);
        for (unsigned int t_ = 0; t_ < V_dim; t_++) 
        {
            V_[t_] = V[t - t_];
            
            if (lag >= (t + 1)) 
                R_(t_, t_) = R(t_, t_);
        }
    }
    
    arma::colvec f(const unsigned int & t, const arma::colvec & theta_) 
    {
        return theta_;
    }
    
    arma::colvec h(const unsigned int & t, const arma::colvec & theta_) 
    {
        const double & V_OC = 2.0 + theta_[0];
        arma::colvec V_T = arma::zeros(V_dim);
        
        for (unsigned int t_ = 0; t_ < V_dim; t_++) 
        {
            const double R_0 = std::exp(theta_[S - 1]);
            
            const double & V_thevenin = V_OC + R_0 * I[t - t_];
            V_T[t_] += V_thevenin;
            
            for (unsigned int k = 0; k < K; k++)
            {
                double tau_k, R_k, V_k_ = 0.0;
                
                tau_k = std::exp(theta_[1 + k]);
                R_k = std::exp(theta_[K + 1 + k]);
                
                if ((t - t_) > 0) 
                {
                    V_k_ = V_k(k, t - t_ - 1);
                }
                
                const double & Vk_ = std::exp(-dt / tau_k);
                const double V_k_1 = V_k_ * Vk_ + R_k * (1.0 - Vk_) * I[t - t_];
                
                V_T[t_] += V_k_1;
                
                if (t_ == 0)
                {
                    V_k(k, t) = V_k_1;
                }
            }
        }
        
        // for (unsigned int t_ = 0; t_ < V_dim; t_++)
        // {
        //     double R_0;
        //     if (t_ == 0) 
        //     {
        //         R_0 = std::exp(theta_[S - 1]);
        //     }
        //     else 
        //     {
        //         R_0 = std::exp(theta_trace(S - 1, t - t_));
        //     }
        //     
        //     const double & V_thevenin = V_OC + R_0 * I[t - t_];
        //     V_T[t_] += V_thevenin;
        //     
        //     for (unsigned int k = 0; k < K; k++)
        //     {
        //         double tau_k, R_k, V_k_ = 0.0;
        //         if (t_ == 0) 
        //         {
        //             tau_k = std::exp(theta_[1 + k]);
        //             R_k = std::exp(theta_[K + 1 + k]);
        //         }
        //         else 
        //         {
        //             tau_k = std::exp(theta_trace(1 + k, t - t_));
        //             R_k = std::exp(theta_trace(K + 1 + k, t - t_));
        //         }
        //         
        //         if ((t - t_) > 0) 
        //         {
        //             V_k_ = V_k(k, t - t_ - 1);
        //         }
        //         
        //         const double & Vk_ = std::exp(-dt / tau_k);
        //         const double V_k_1 = V_k_ * Vk_ + R_k * (1.0 - Vk_) * I[t - t_];
        //         
        //         V_T[t_] += V_k_1;
        //         
        //         if (t_ == 0)
        //         {
        //             V_k(k, t) = V_k_1;
        //         }
        //     }
        // }
        
        return V_T; 
    }
    
    void calculate_sigma_points(const unsigned int & t, const arma::colvec & theta_) 
    {
        theta_tilde = arma::mat(S, 2 * S + 1);
        theta_tilde.col(0) = theta_;
        
        arma::mat P_theta_sqroot = std::pow(S + lambda, 0.5) * arma::sqrtmat_sympd(P_theta[t]); 
        for (unsigned int s = 0; s < S; s++) 
        {
            theta_tilde.col(s + 1) = theta_ + P_theta_sqroot.col(s);
            theta_tilde.col(S + s + 1) = theta_ - P_theta_sqroot.col(s);
        }
    }
    
    void update_theta_tilde(const unsigned int & t) 
    {
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            theta_tilde.col(s) = f(t, theta_tilde.col(s));
        }
    }
    
    void update_VT_tilde(const unsigned int & t) 
    {
        VT_tilde = arma::mat(V_dim, 2 * S + 1);
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            VT_tilde.col(s) = h(t, theta_tilde.col(s));
        }
    }
    
    // Filtering
    void update_theta_VT_filter(arma::colvec & theta_, arma::colvec & VT_) 
    {
        theta_ = arma::zeros(S);
        VT_ = arma::zeros(V_dim);
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            theta_ += W_m[s] * theta_tilde.col(s);
            VT_ += W_m[s] * VT_tilde.col(s);
        }
    }
    
    void update_covariance_matrices_filter(arma::mat & P_theta_, arma::mat & P_theta_y, arma::mat & P_y, 
                                     const arma::colvec & theta_, const arma::colvec & VT_, 
                                     const arma::mat & R_) 
    {
        P_theta_ = arma::zeros(S, S);
        P_theta_y = arma::zeros(S, V_dim);
        P_y = arma::zeros(V_dim, V_dim);
        
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            const double W_c_s = W_c[s];
            const arma::colvec theta_error = theta_tilde.col(s) - theta_;
            const arma::colvec VT_error = VT_tilde.col(s) - VT_;
            
            P_theta_ += W_c_s * (theta_error * theta_error.t());
            P_theta_y += W_c_s * (theta_error * VT_error.t());
            P_y += W_c_s * (VT_error * VT_error.t());
        }
        
        // P_theta_ += Q;
        // P_y += R_;
        P_y = 0.5 * P_y + 0.5 * P_y;
        P_y = P_y + 1e-6 * arma::eye(P_y.n_rows, P_y.n_cols);
    }
    
    // Smoothing
    void ma_theta(const unsigned int & t, arma::colvec & theta_) 
    {
        int t_ = static_cast<int>(t);
        int ML = (t_ - M);
        if (ML < 0) {
            ML = 0;
        }
        
        int MU = (t_ + M);
        if (MU > (T - 1)) {
            MU = T - 1;
        }
        
        ML = t_ - ML;
        MU = MU - t_;
        
        double normaliser = 1.0;
        for (unsigned int m = 1; m < std::max(ML, MU); m++) 
        {
            if (m < ML) 
            {
                const double & m_ = static_cast<double>(m);
                const arma::colvec theta_m = theta_trace.col(t - m);
                normaliser += 1.0; // (1.0 / (m_ + 1.0));
                theta_ += theta_m; // (theta_m / (m_ + 1.0));
            }
            
            if (m < MU) 
            {
                const double & m_ = static_cast<double>(m);
                const arma::colvec theta_m = theta_trace.col(t + m);
                normaliser += 1.0; // (1.0 / (m_ + 1.0));
                theta_ += theta_m; // (theta_m / (m_ + 1.0));
            }
        }
        
        theta_ = theta_ / normaliser;
    }
    
    void update_theta_smoother(arma::colvec & theta_) 
    {
        theta_ = arma::zeros(S);
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            theta_ += W_m[s] * theta_tilde.col(s);
        }
    }
    
    void update_covariance_matrices_smoother(arma::mat & P_theta_, arma::mat & D_,
                                             const unsigned int & t_, const arma::mat & theta_tilde_) 
    {
        P_theta_ = arma::zeros(S, S);
        D_ = arma::zeros(S, S);
        
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            const double W_c_s = W_c[s];
            const arma::colvec theta_t = theta_trace.col(t_);
            const arma::colvec theta_t_ = theta_trace.col(t_ + 1);
            
            const arma::colvec theta_error_t_ = theta_tilde.col(s) - theta_t_;
            const arma::colvec theta_error_t = theta_tilde_.col(s) - theta_t;
            
            P_theta_ += W_c_s * (theta_error_t_ * theta_error_t_.t());
            D_ += W_c_s * (theta_error_t * theta_error_t.t());
        }
        
        P_theta_ += Q; 
        P_theta_ = 0.5 * P_theta_ + 0.5 * P_theta_;
        P_theta_ = P_theta_ + 1e-6 * arma::eye(P_theta_.n_rows, P_theta_.n_cols);
    }
    
public:
    //
    arma::colvec theta, theta_, VT, SOC, W_c, W_m, P_y_trace;
    arma::mat theta_tilde, VT_tilde, theta_trace;
    std::vector<arma::mat> P_theta;
    
    //
    ParameterKF(const arma::colvec & I_, const arma::colvec & V_, const arma::colvec & theta_, 
                const arma::mat & P_theta_, const double & dt_, const unsigned int K_,
                const double SOC_0_, const double C_max_, const double & eta_,
                const double & gamma_, const double & xi_, const double & kappa_, 
                const unsigned int & lag_, const arma::mat & R_0_, const arma::mat & Q_0_, 
                const unsigned int & M_, const bool & trace_, const unsigned int & traceLimit_) :
        I(I_), V(V_), theta(theta_), dt(dt_), K(K_), 
        trace(trace_), traceLimit(traceLimit_), C_max(C_max_), eta(eta_), 
        gamma(gamma_), xi(xi_), kappa(kappa_), 
        R(R_0_), Q(Q_0_), lag(lag_), M(M_), N(1)
    {
        T = I.size();
        S = theta.size();
        
        SOC = arma::colvec(T);
        SOC[0] = SOC_0_;

        VT = arma::colvec(T);
        P_y_trace = arma::colvec(T);
        theta_trace = arma::mat(S, T);
        V_k = arma::mat(K, T);
        
        // Set-up for sigma-points and UT-weights
        const double & gamma_2 = gamma * gamma;
        lambda = (S + kappa) * gamma_2  - S;
        
        W_m = arma::colvec(2 * S + 1);
        W_c = arma::colvec(2 * S + 1);
        
        W_m[0] = lambda / (S + lambda);
        W_c[0] = W_m[0] + (1.0 - gamma_2 + xi);
        
        for (unsigned int s = 0; s < 2 * S; s++) 
        {
            double weight_s = 1.0 / (2 * (S + lambda));
            W_c[s + 1] = weight_s;
            W_m[s + 1] = weight_s;
        }
        
        theta_tilde = arma::mat(S, 2 * S + 1);
        
        // 
        P_theta = std::vector<arma::mat>(T + 1);
        P_theta[0] = P_theta_;
    }
    
    // 
    void Filter() 
    {
        // 
        arma::colvec theta_, V_lag, VT_;
        arma::mat P_theta_, P_theta_y, P_y, R_;
        
        //
        VT[0] = h(0, theta)(0, 0);
        
        // 
        for (unsigned int t = 1; t < T; t++) // 
        {
            if (trace & ((t == 0) | (((t + 1) % traceLimit) == 0) | (t == (T - 1))))
            {
                Rcpp::Rcout << "Iteration: " << t + 1 << " / " << T << "\n";
            }
            
            //// Updating SOC and V_lag
            update_SOC(t);
            update_V(t, V_lag, R_);
            
            //// Sigma points
                calculate_sigma_points(t, theta);
                
                //// Time update
                update_theta_tilde(t);
                update_VT_tilde(t);
                
                update_theta_VT_filter(theta_, VT_);
                update_covariance_matrices_filter(P_theta_, P_theta_y, P_y, theta_, VT_, R_);
                P_y_trace[t] = P_y(0, 0);
                
                //// Measurement update
                const arma::mat P_y_inverse = arma::inv(P_y);
                const arma::mat K = P_theta_y * P_y_inverse;
                
                theta = theta_ + K * (V_lag - VT_);
                
                P_theta[t + 1] = P_theta_ - K * P_y * K.t();
                P_theta[t + 1] = 0.5 * P_theta[t + 1] + 0.5 * P_theta[t + 1].t();
                P_theta[t + 1] = P_theta[t + 1] + 1e-6 * arma::eye(P_theta[t + 1].n_rows, P_theta[t + 1].n_cols);
            
            //// Logging
            VT[t] = h(t, theta)(0, 0);
            theta_trace.col(t) = theta;
            
            if (trace & ((t == 0) | (((t + 1) % traceLimit) == 0) | (t == (T - 1))))
            {
                Rcpp::Rcout << "\tTheta_:\n\t" << theta_.t()
                            << "\tTheta:\n\t" << theta.t()
                            << "\tV = " << V[t] << " :: V_hat = " << VT[t] << " :: I = " << I[t] << "\n"
                            << "-------------------------------\n";
            }
        }
    }
    
    void Smoother() 
    {
        double VT_;
        arma::colvec theta_ma, theta_, theta_t, V_;
        arma::mat P_theta_, D_, theta_tilde_, R_;
        
        unsigned int t_;
        for (unsigned int t = T - 1; t >= 1; t--) 
        {
            t_ = t - 1;
            VT_ = VT[t_];
            
            if (trace & ((t_ == 0) | (((t_ + 1) % traceLimit) == 0) | (t_ == (T - 2))))
            {
                Rcpp::Rcout << "Iteration: " << t_ + 1 << " / " << T << "\n";
            }
            
            update_V(t_, V_, R_);
            theta_t = theta_trace.col(t_);
            theta_ma = theta_t;
            ma_theta(t, theta_ma);
            
            ////
            calculate_sigma_points(t_, theta_ma);
            theta_tilde_ = theta_tilde;
            
            update_theta_smoother(theta_);
            update_covariance_matrices_smoother(P_theta_, D_, t_, theta_tilde_);
            
            ////
            const arma::mat P_theta_inverse = arma::inv(P_theta_);
            const arma::mat K = D_ * P_theta_inverse;
            
            const arma::colvec theta_s = theta_t + K * (theta_t - theta_);
            const arma::mat P_theta_s = P_theta[t_] - K * (P_theta[t_ + 1] - P_theta_) * K.t();
            
            theta_trace.col(t_) = theta_s;
            P_theta[t_] = P_theta_s;
            
            ////
            VT[t_] = h(t_, theta_s)(0, 0);
            
            if (trace & ((t_ == 0) | (((t_ + 1) % traceLimit) == 0) | (t_ == (T - 2))))
            {
                Rcpp::Rcout << "\tTheta_t:\n\t" << theta_t.t()
                            << "\tTheta_:\n\t" << theta_.t()
                            << "\tTheta_t^s:\n\t" << theta_s.t()
                            << "\tV = " << V[t_] << " :: V_hat = " << VT_ 
                            << " :: V_hat^s = " << VT[t_] << " :: I = " << I[t_] << "\n"
                            << "-------------------------------\n";
            }
        }
    }
};

//[[Rcpp::export]] 
Rcpp::List ParameterUKFCpp(const arma::colvec & I, const arma::colvec & V, const arma::colvec & theta_0, const arma::mat & P_0,
                           const double & SOC_0, const double & C_max, const double & eta,
                           const arma::colvec & sigma_point_pars, const unsigned int & lag,
                           const double & dt, const unsigned int K, const arma::mat & R_0, const arma::mat & Q_0,
                           const bool & trace, const unsigned int & traceLimit) 
{
    ParameterKF PKF(I, V, theta_0, P_0, dt, K, SOC_0, C_max, eta, 
                    sigma_point_pars[0], sigma_point_pars[1], sigma_point_pars[3], 
                    lag, R_0, Q_0, 0, trace, traceLimit);
    PKF.Filter();
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = PKF.VT,
                              Rcpp::Named("SOC") = PKF.SOC,
                              Rcpp::Named("Theta") = PKF.theta, 
                              Rcpp::Named("ThetaTrace") = PKF.theta_trace);
}


//[[Rcpp::export]] 
Rcpp::List ParameterUKSCpp(const arma::colvec & I, const arma::colvec & V, const arma::colvec & theta_0, const arma::mat & P_0,
                           const double & SOC_0, const double & C_max, const double & eta,
                           const arma::colvec & sigma_point_pars, const unsigned int & lag, const unsigned int & M,
                           const double & dt, const unsigned int K, const arma::mat & R_0, const arma::mat & Q_0,
                           const bool & trace, const unsigned int & traceLimit) 
{
    ParameterKF PKF(I, V, theta_0, P_0, dt, K, SOC_0, C_max, eta, 
                    sigma_point_pars[0], sigma_point_pars[1], sigma_point_pars[3],
                    lag, R_0, Q_0, M, trace, traceLimit);
    
    PKF.Filter();
    PKF.Smoother();
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = PKF.VT,
                              Rcpp::Named("SOC") = PKF.SOC,
                              Rcpp::Named("Theta") = PKF.theta, 
                              Rcpp::Named("ThetaTrace") = PKF.theta_trace);
}
