#include <RcppArmadillo.h>

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/math/constants/constants.hpp>


template <typename T> int sgn(T val) 
{
    return(T(0) < val) - (val < T(0));
} 

arma::mat square_root_arma_matrix(const arma::mat & M) 
{
    arma::mat U, V;
    arma::vec s;
    
    arma::svd(U, s, V, M);
    
    arma::mat S = arma::diagmat(arma::pow(s, 0.5));
    arma::mat X = U * S * V.t();
    return X;
}


class ParameterKF {
private:
    //
    arma::colvec I, V;
    double dt, eta, C_max; // 
    bool trace;
    unsigned int traceLimit, S, N, K, T;
    
    void update_SOC(const double & t) 
    {
        SOC[t + 1] = SOC[t] + dt * (I[t] / (3600 * C_max));
    }
    
    arma::colvec update_VK(const double & t, const arma::colvec & VK_) 
    {
        arma::colvec VK_1(VK.size());
        for (unsigned int k = 0; k < K; k++)
        {
            const double & R_k = theta[2 * k + 1];
            const double & omega_k = theta[2 * k + 2];
            VK_1[k + 1] = VK_[k] + dt * (omega_k * R_k * I[t] - omega_k * VK_[k]);
        }
        
        return VK_;
    }
    
    double h(const unsigned int & t, const arma::colvec & theta_, const arma::colvec & VK_) 
    {
        const double & V_OC_ = theta_[S - 1];
        const double & R0 = theta_[0];
        double V_ = V_OC_ + I[t] * R0; // V_OC(SOC[t + 1])
        
        for (unsigned int k = 0; k < K; k++)
        {
            V_ += VK_[k + 1];
        }
        
        return V_;
    }
    
    void calculate_sigma_points() 
    {
        theta_tilde = arma::mat(S, 2 * S + 1);
        theta_tilde.col(0) = theta;
        
        arma::mat P_theta_sqroot = square_root_arma_matrix((S + lambda) * P_theta);
        for (unsigned int s = 0; s < S; s++) 
        {
            theta_tilde.col(s + 1) = theta + P_theta_sqroot.col(s);
            theta_tilde.col(S + s + 1) = theta - P_theta_sqroot.col(s);
        }
    }
    
    void update_theta_tilde() 
    {
        // theta_tilde = theta_tilde;
    }
    
    void update_VT_tilde(const unsigned int & t, const arma::colvec & VK_) 
    {
        for (unsigned int s = 0; s < (2 * S + 1); s++) 
        {
            VT_tilde.col(s) = h(t, theta_tilde.col(s), VK_);
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
            
            P_theta += W_c_s * (theta_error * theta_error.t());
            P_theta_y += W_c_s * (theta_error * VT_error);
            P_y += W_c_s * (VT_error * VT_error);
        }
    }
    
public:
    //
    arma::colvec theta, theta_, VT, SOC, VK, W_c, W_m;
    arma::mat P_theta, theta_tilde, VT_tilde, theta_trace;
    double lambda;
    
    //
    ParameterKF(const arma::colvec & I_, const arma::colvec & V_, const arma::colvec & theta_, 
                const arma::mat & P_theta_, const double & dt_, const unsigned int K_,
                const double SOC_0_, const double C_max_, const double & eta_,
                const bool & trace_, const unsigned int & traceLimit_) :
        I(I_), V(V_), theta(theta_), P_theta(P_theta_), dt(dt_), K(K_), 
        trace(trace_), traceLimit(traceLimit_), C_max(C_max_), eta(eta_), N(1)
    {
        T = I.size();
        S = theta.size();
        
        SOC = arma::colvec(T + 1);
        SOC[0] = SOC_0_;

        VK = arma::zeros(K);
        
        VT = arma::colvec(T);
        theta_trace = arma::mat(S, T);
        
        // Set-up for sigma-points and UT-weights
        const double alpha_tilde = 0.01;
        const double beta_tilde = 2.0;
        const double kappa = 10.0;
        
        const double & alpha_tilde_2 = alpha_tilde * alpha_tilde;
        lambda = (S + kappa) * alpha_tilde_2  - S;
        W_m = arma::colvec(2 * S + 1);
        W_c = arma::colvec(2 * S + 1);
        W_m[0] = lambda / (S + lambda);
        W_c[0] = W_m[0] + (1.0 - alpha_tilde_2 + beta_tilde);
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
        boost::mt19937 rng;
        boost::random::uniform_01<> uniform_real;
        boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> > generate_uniform_real(rng, uniform_real);
        
        // 
        double VT_, P_y;
        arma::colvec theta_, VK_;
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
            VK_ = update_VK(t, VK);
            
            //// Sigma points
            calculate_sigma_points();
            
            //// Time update
            update_theta_tilde();
            update_VT_tilde(t, VK_);
            
            update_theta_VT_(theta_, VT_);
            update_covariance_matrices_(P_theta_, P_theta_y, P_y, theta_, VT_);
            
            //// Measurement update
            const arma::mat K = P_theta_y / P_y;
            theta = theta_ + K * (V[t] - VT_);
            P_theta = P_theta_ - P_y * K * K.t();
            
            //// Logging
            VK = update_VK(t, VK);
            VT[t] = h(t, theta, VK);
            theta_trace.col(t) = theta;
        }
    }
};

//[[Rcpp::export]] 
Rcpp::List ParameterUKFCpp(const arma::colvec & I, const arma::colvec & V, const arma::colvec & theta_0, const arma::mat & P_0,
                           const double & SOC_0, const double & C_max, const double & eta,
                           const double & dt, const unsigned int K, const bool & trace, const unsigned int & traceLimit) 
{
    ParameterKF PKF(I, V, theta_0, P_0, dt, K, SOC_0, C_max, eta, trace, traceLimit);
    PKF.Filter();
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = PKF.VT,
                              Rcpp::Named("SOC") = PKF.SOC,
                              Rcpp::Named("Theta") = PKF.theta, 
                              Rcpp::Named("ThetaTrace") = PKF.theta_trace, 
                              Rcpp::Named("P") = PKF.P_theta);
}
