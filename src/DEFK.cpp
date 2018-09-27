#include <RcppArmadillo.h>

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

class DEKFModel {
    private:
        //
        arma::colvec I, V;
        double dt, eta, C_max;
        bool trace;
        unsigned int traceLimit, S, N, K, T;
        
        // 
        double V_OC(const double & SOC_) 
        {
            const double p0 = 3.1345894764792;
            const double p1 = 0.0090247123562;
            const double p2 = -0.0001458437195;
            const double p3 = 0.0000007995438;
            
            const double SOC_2 = SOC_ * SOC_;
            const double SOC_3 = SOC_ * SOC_2;
            return p0 + p1 * SOC_ + p2 * SOC_2 + p3 * SOC_3;
        }
        
        double V_OC_(const double & SOC_) 
        {
            const double p1 = 0.0090247123562;
            const double p2 = -0.0001458437195;
            const double p3 = 0.0000007995438;
            
            const double SOC_2 = std::pow(SOC_, 2.0);
            return p1 + 2.0 * p2 * SOC_ + 3.0 * p3 * SOC_2;
        }
        
        double h(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
        {
            double V_ = V_OC(X_[0]) + I[t] * std::exp(theta_[0]);
            for (unsigned int k = 0; k < K; k++)
            {
                V_ += X_[k + 1];
            }
            
            return V_;
        }
        
        arma::colvec f(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
        {
            arma::colvec X_new = arma::zeros(X_.size());
            X_new[0] = X_[0] + eta * dt / C_max * I[t];
            
            for (unsigned int k = 0; k < K; k++) 
            {
                const double & Vk_ = std::exp(-dt / theta_[2 * k + 2]);
                X_new[k + 1] = X_[k + 1] * Vk_ + 
                    std::exp(theta_[1 + 2 * k]) * (1.0 - Vk_) * I[t];
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
            H_x(0, 0) = V_OC_(X[0]);
        }
        
        void update_H_theta(const unsigned int & t) 
        {
            arma::mat H = arma::zeros(1, S);
            H(0, 0) = I[t];
            
            arma::mat x_theta = arma::zeros(K + 1, S);
            for (unsigned int k = 0; k < K; k++) 
            {
                const double tau_squared = std::pow(theta[2 * k + 2], 2.0);
                
                x_theta(k + 1, 2 * k + 1) = I[t] * (std::exp(dt / tau_squared) - 1.0);
                x_theta(k + 1, 2 * k + 2) = (dt / tau_squared) *
                    (X[k + 1] + std::exp(theta[2 * k + 1]) * I[t]) * std::exp(-dt / theta[2 * k + 2]);
            }
            
            arma::mat H_ = H_x * x_theta;
            H_theta = H + H_;
        }
        
    public:
        //
        arma::colvec theta, X, VT, SOC;
        arma::mat F, H_x, H_theta, Q_x, Q_theta, P_x, P_theta, R_x, R_theta, theta_trace;

        //
        DEKFModel(const arma::colvec & I_, const arma::colvec & V_, const arma::colvec theta_,
                  const double & dt_, const unsigned int K_, const double SOC_intial_,
                  const double C_max_, const bool & trace_, const unsigned int & traceLimit_) :
        I(I_), V(V_), theta(theta_), dt(dt_), K(K_), trace(trace_), traceLimit(traceLimit_), N(1)
        {
            C_max = 36 * C_max_;
            eta = 1.0;
            
            T = I.size();
            S = theta.size();
            
            VT = arma::colvec(T);
            SOC = arma::colvec(T);
            theta_trace = arma::mat(S, T);
            
            X = arma::zeros(K + 1);
            X[0] = SOC_intial_;
            
            Q_x = 1e-8 * arma::eye(K + 1, K + 1);
            P_x = 10.0 * arma::eye(K + 1, K + 1);
            R_x = 10.0 * arma::eye(N, N);
            
            Q_theta = 1e-8 * arma::eye(S, S);
            P_theta = 10.0 * arma::eye(S, S);
            R_theta = 10.0 * arma::eye(N, N);
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
            }
        }
        
        arma::vec Predict();
};


//' @title Filtering of state-of-charge and health.
//' 
//' @description Dual extended Kalman filter for determining the latent state-of-charge and health state-equation. 
//' 
//' @param I A current profile.
//' @param V The observed voltage.
//' @param theta_0 The initial values of the parameter vector.
//' @param SOC_0 The initial value of the SOC.
//' @param C_max The maximum value of the capacity.
//' @param dt Step-size.
//' @param K Number of RC branches in the EEC.
//' @param trace TRUE/FALSE: Show trace?
//' @param traceLimit Used to limit trace.
//' 
//' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
//' @export
//[[Rcpp::export]] 
Rcpp::List SODEFK(const arma::colvec & I, const arma::colvec & V, const arma::colvec & theta_0, const double & SOC_0,
                  const double & C_max, const double & dt, const unsigned int K, const bool & trace, const unsigned int & traceLimit) 
{
    DEKFModel DEKF(I, V, theta_0, dt, K, SOC_0, C_max, trace, traceLimit);
    DEKF.Filter();
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = DEKF.VT,
                              Rcpp::Named("SOC") = DEKF.SOC,
                              Rcpp::Named("Theta") = DEKF.theta, 
                              Rcpp::Named("ThetaTrace") = DEKF.theta_trace);
}
