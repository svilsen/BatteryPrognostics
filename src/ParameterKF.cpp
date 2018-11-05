#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

class ParameterKF {
private:
    //
    arma::colvec I, V;
    double dt, eta, C_max, gamma, xi, kappa, lambda, alpha; // 
    bool trace;
    unsigned int traceLimit, S, N, K, T;
    
    void update_SOC(const double & t) 
    {
        SOC[t + 1] = SOC[t] + dt * (I[t] / (3600 * C_max));
    }
    
    arma::colvec f(const unsigned int & t, const arma::colvec & theta_) 
    {
        arma::colvec f_theta = theta_;
        const double & I_t = I[t];
        
        const double & R_0 = std::exp(theta_[S - 1]);
        for (unsigned int k = 0; k < K; k++)
        {
            const double & V_k = theta_[k + 1];
            
            const double & R_k = std::exp(theta_[2 * K + 1 + k]);
            const double & tau_k = std::exp(theta_[K + 1 + k]);
            
            const double & Vk_ = std::exp(-dt / tau_k); 
            f_theta[k + 1] = V_k * Vk_ + R_k * (1.0 - Vk_) * I_t;
        }

        return f_theta;
    }
    
    double h(const unsigned int & t, const arma::colvec & theta_) 
    {
        const double & V_OC = theta_[0];
        const double & R_0 = std::exp(theta_[S - 1]);
        
        double V_T = V_OC + R_0 * I[t];
        for (unsigned int k = 0; k < K; k++)
        {
            const double & V_k = theta_[k + 1];
            V_T += V_k;
        }
        
        return V_T; 
    }
    
    void calculate_sigma_points() 
    {
        theta_tilde = arma::mat(S, 2 * S + 1);
        theta_tilde.col(0) = theta;
        
        arma::mat P_theta_sqroot = std::pow(S + lambda, 0.5) * arma::sqrtmat_sympd(P_theta); 
        for (unsigned int s = 0; s < S; s++) 
        {
            theta_tilde.col(s + 1) = theta + P_theta_sqroot.col(s);
            theta_tilde.col(S + s + 1) = theta - P_theta_sqroot.col(s);
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
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            VT_tilde.col(s) = h(t, theta_tilde.col(s));
        }
    }
    
    void update_theta_VT_(arma::colvec & theta_, double & VT_) 
    {
        theta_ = arma::zeros(S, 1);
        VT_ = 0.0;
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            theta_ += W_m[s] * theta_tilde.col(s);
            VT_ += W_m[s] * VT_tilde(0, s);
        }
    }
    
    void update_covariance_matrices_(arma::mat & P_theta_, arma::mat & P_theta_y, double & P_y, 
                                     const arma::colvec & theta_, const double & VT_t) 
    {
        P_theta_ = arma::zeros(S, S);
        P_theta_y = arma::zeros(S, 1);
        P_y = 0.0;
        
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            const double W_c_s = W_c[s];
            const arma::colvec theta_error = theta_tilde.col(s) - theta_;
            const double VT_error = VT_tilde(0, s) - VT_t;
            
            P_theta_ += W_c_s * (theta_error * theta_error.t());
            P_theta_y += W_c_s * (theta_error * VT_error);
            P_y += W_c_s * (VT_error * VT_error);
        }
    }
    
public:
    //
    arma::colvec theta, theta_, VT, SOC, W_c, W_m, P_y_trace;
    arma::mat P_theta, theta_tilde, VT_tilde, theta_trace;
    
    //
    ParameterKF(const arma::colvec & I_, const arma::colvec & V_, const arma::colvec & theta_, 
                const arma::mat & P_theta_, const double & dt_, const unsigned int K_,
                const double SOC_0_, const double C_max_, const double & eta_,
                const double & gamma_, const double & xi_, const double & kappa_, const double & alpha_,
                const bool & trace_, const unsigned int & traceLimit_) :
        I(I_), V(V_), theta(theta_), P_theta(P_theta_), dt(dt_), K(K_), 
        trace(trace_), traceLimit(traceLimit_), C_max(C_max_), eta(eta_), 
        gamma(gamma_), xi(xi_), kappa(kappa_), alpha(alpha_), N(1)
    {
        T = I.size();
        S = theta.size();
        
        SOC = arma::colvec(T + 1);
        SOC[0] = SOC_0_;

        VT = arma::colvec(T);
        P_y_trace = arma::colvec(T);
        theta_trace = arma::mat(S, T);
        
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
        VT_tilde = arma::mat(1, 2 * S + 1);
    }
    
    // 
    void Filter() 
    {
        // 
        double VT_, P_y;
        arma::colvec theta_;
        arma::mat P_theta_, P_theta_y;
        
        // 
        for (unsigned int t = 0; t < T; t++) 
        {
            if (trace & ((t == 0) | (((t + 1) % traceLimit) == 0) | (t == (T - 1))))
            {
                Rcpp::Rcout << "Iteration: " << t + 1 << " / " << T << "\n";
            }
            
            //// Updating SOC
            update_SOC(t);
            
            //// Sigma points
            calculate_sigma_points();
            
            //// Time update
            update_theta_tilde(t);
            update_VT_tilde(t);
            
            update_theta_VT_(theta_, VT_);
            update_covariance_matrices_(P_theta_, P_theta_y, P_y, theta_, VT_);
            
            P_y_trace[t] = P_y;
            
            //// Measurement update
            const arma::mat K = P_theta_y / P_y;
            theta = theta_ + K * (V[t] - VT_);
            theta[0] = alpha * theta_[0] + (1 - alpha) * theta[0];
            
            P_theta = P_theta_ - P_y * K * K.t();
            P_theta = 0.5 * P_theta + 0.5 * P_theta.t();
            P_theta = P_theta + 1e-6 * arma::eye(P_theta.n_rows, P_theta.n_cols);
            
            //// Logging
            VT[t] = h(t, theta);
            theta_trace.col(t) = theta;
            
            if (trace & ((t == 0) | (((t + 1) % traceLimit) == 0) | (t == (T - 1))))
            {
                Rcpp::Rcout << "\tTheta:\n\t" << theta.t()
                            << "\tV = " << V[t] << " :: V_hat = " << VT[t] << "\n"
                            << "-------------------------------\n";
            }
        }
    }
};

//[[Rcpp::export]] 
Rcpp::List ParameterUKFCpp(const arma::colvec & I, const arma::colvec & V, const arma::colvec & theta_0, const arma::mat & P_0,
                           const double & SOC_0, const double & C_max, const double & eta,
                           const arma::colvec & sigma_point_pars, const double & alpha,
                           const double & dt, const unsigned int K, const bool & trace, const unsigned int & traceLimit) 
{
    ParameterKF PKF(I, V, theta_0, P_0, dt, K, SOC_0, C_max, eta, 
                    sigma_point_pars[0], sigma_point_pars[1], sigma_point_pars[3], 
                    alpha, trace, traceLimit);
    PKF.Filter();
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = PKF.VT,
                              Rcpp::Named("SOC") = PKF.SOC,
                              Rcpp::Named("Theta") = PKF.theta, 
                              Rcpp::Named("ThetaTrace") = PKF.theta_trace, 
                              Rcpp::Named("P_theta") = PKF.P_theta, 
                              Rcpp::Named("P_y") = PKF.P_y_trace);
}
