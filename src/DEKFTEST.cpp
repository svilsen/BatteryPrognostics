#include <RcppArmadillo.h>

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

class DEKFModel {
private:
    //
    arma::colvec I, V, ocv_parameters;
    double dt, C_max; // eta, 
    bool trace;
    unsigned int traceLimit, S, N, K, T;
    
    double h(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
    {
        const double & R0 = std::exp(theta_[0]);
        double V_ = theta_[theta_.size() - 1] + I[t] * R0;
        for (unsigned int k = 0; k < K; k++)
        {
            V_ += X_[k + 1];
        }
        
        return V_;
    }
    
    arma::colvec f(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
    {
        arma::colvec X_new = arma::zeros(X_.size());
        X_new[0] = X_[0] + theta_[S - 2] * dt / C_max * I[t];
        
        for (unsigned int k = 0; k < K; k++) 
        {
            const double & Rk = std::exp(theta[2 * k + 1]);
            const double & tau_k = theta[2 * k + 2];
            
            const double & Vk_ = std::exp(-dt / tau_k);
            X_new[k + 1] = X_[k + 1] * Vk_ + Rk * (1.0 - Vk_) * I[t];
        }
        
        return X_new;
    }
    
    void update_F(const unsigned int & t) 
    {
        F = arma::zeros(K + 1, K + 1);
        F(0, 0) = 1.0;
        for (unsigned int k = 0; k < K; k++) 
        {
            F(k + 1, k + 1) = std::exp(-dt / theta[2 * k + 2]);
        }
    }
    
    void update_H_x() 
    {
        H_x = arma::ones(1, K + 1);
        //H_x(0, 0) = V_OC_(X[0]);
        H_x(0, 0) = 0.0;
    }
    
    void update_H_theta(const unsigned int & t) 
    {
        arma::mat H = arma::zeros(1, S);
        H(0, 0) = I[t];
        H(0, S - 1) = 1;
        
        arma::mat x_theta = arma::zeros(K + 1, S);
        for (unsigned int k = 0; k < K; k++) 
        {
            const double Rk = std::exp(theta[2 * k + 1]);
            const double tau_k = theta[2 * k + 2];
            const double tau_squared = tau_k * tau_k;
            
            x_theta(k + 1, 2 * k + 1) = I[t] * (std::exp(dt / tau_squared) - 1.0);
            x_theta(k + 1, 2 * k + 2) = (dt / tau_squared) *
                (X[k + 1] + Rk * I[t]) * std::exp(-dt / tau_k);
        }
        
        x_theta(0, S - 2) = dt * I[t] / C_max;
        
        arma::mat H_ = H_x * x_theta;
        H_theta = H + H_;
    }
    
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
              const bool & trace_, const unsigned int & traceLimit_) :
        I(I_), V(V_), theta(theta_), ocv_parameters(ocv_parameters_), dt(dt_), K(K_), trace(trace_), traceLimit(traceLimit_), N(1)
    {
        C_max = 3600 * C_max_;
        // eta = 1.0;
        
        T = I.size();
        S = theta.size();
        
        VT = arma::colvec(T);
        SOC = arma::colvec(T);
        theta_trace = arma::mat(S, T);
        
        X = arma::zeros(K + 1);
        X[0] = SOC_intial_;
        
        Q_x = Q[0];
        P_x = P[0];
        R_x = R[0];
        
        Q_theta = Q[1];
        P_theta = P[1];
        R_theta = R[1];
    }
    
    // 
    void Filter() 
    {
        boost::mt19937 rng;
        boost::random::uniform_01<> uniform_real;
        boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> > generate_uniform_real(rng, uniform_real);
        
        double V_hat;
        arma::colvec X_, theta_;
        arma::mat P_x_, P_theta_, L_x, L_theta, S_x, S_theta, HL_x, HL_theta;
        LogLikelihood = 0;
        for (unsigned int t = 0; t < T; t++) 
        {
            if (trace & ((t == 0) | ((t % traceLimit) == 0) | (t == (T - 1))))
            {
                Rcpp::Rcout << "Iteration: " << t + 1 << " / " << T << "\n";
            }
            
            update_F(t);
            
            //// Time-update
            theta_ = theta;
            X_ = f(t, X, theta_);
            
            P_x_ = F * P_x * F.t() + Q_x;
            P_theta_ = P_theta + Q_theta;
            
            //// Measurement-update
            // X
            update_H_x();
            
            S_x = H_x * P_x_ * H_x.t() + R_x;
            
            L_x = P_x_ * H_x.t() * arma::inv(S_x);
            
            V_hat = h(t, X_, theta_);
            X = X_ + L_x * (V[t] - V_hat);
            
            HL_x = arma::eye(K + 1, K + 1) - L_x * H_x;
            P_x = HL_x * P_x_ * HL_x.t() + L_x * R_x * L_x.t();
            
            // Theta
            update_H_theta(t);
            
            S_theta = H_theta * P_theta_ * H_theta.t() + R_theta;
            L_theta = P_theta_ * H_theta.t() * arma::inv(S_theta);
            
            theta = theta_ + L_theta * (V[t] - V_hat);
            
            HL_theta = arma::eye(S, S) - L_theta * H_theta;
            P_theta = HL_theta * P_theta_ * HL_theta.t() + L_theta * R_theta * L_theta.t();
            
            VT[t] = h(t, X, theta);
            SOC[t] = X[0];
            theta_trace.col(t) = theta;
            
            const double & error = V[t] - VT[t];
            LogLikelihood += (std::abs(S_x(0, 0)) + (error * error) / S_x(0, 0));
        }
    }
};

//[[Rcpp::export]] 
Rcpp::List SODEFKTESTCpp(const arma::colvec & I, const arma::colvec & V, const arma::colvec & theta_0,
                         const arma::colvec & ocv_parameters, const double & SOC_0, const double & C_max,
                         const std::vector<arma::mat> & P, const std::vector<arma::mat> & Q, const std::vector<arma::mat> & R,
                         const double & dt, const unsigned int K, const bool & trace, const unsigned int & traceLimit) 
{
    DEKFModel DEKF(I, V, theta_0, ocv_parameters, P, Q, R, dt, K, SOC_0, C_max, trace, traceLimit);
    DEKF.Filter();
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = DEKF.VT,
                              Rcpp::Named("SOC") = DEKF.SOC,
                              Rcpp::Named("Theta") = DEKF.theta, 
                              Rcpp::Named("ThetaTrace") = DEKF.theta_trace, 
                              Rcpp::Named("LogLikelihood") = DEKF.LogLikelihood);
}
