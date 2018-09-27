#include <RcppArmadillo.h>

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

class EFKModel {
    private:
        //
        arma::colvec I, IF, V;
        double dt, eta, Cap_max, pi;
        bool trace;
        unsigned int traceLimit, N, K, T;
        
        // 
        double V_OC(const double & SOC_, const unsigned int & t) 
        {
            //
            double IF_sign = 0.0; 
            if (IF[t] > 0.0) 
            {
                IF_sign++; 
            }
            else 
            {
                IF_sign--;
            }
            
            //
            double _SOC_ = SOC_ / 100.0;;
            if (_SOC_ < 0.0) 
            {
                _SOC_ = 0.0;
            }
            else if (_SOC_ > 1.0) 
            {
                _SOC_ = 1.0;
            }
            
            const arma::colvec & parameters_OCV = parameters[2 * K + 2];
            
            const double SOC_2 = _SOC_ * _SOC_;
            const double SOC_3 = _SOC_ * SOC_2;
            return parameters_OCV[0] + parameters_OCV[1] * _SOC_ + parameters_OCV[2] * SOC_2 + parameters_OCV[3] * SOC_3 + IF_sign * parameters_OCV[4];
        }
        
        double V_OC_(const double & SOC_) 
        {
            double _SOC_ = SOC_ / 100.0;;
            if (_SOC_ < 0.0) 
            {
                _SOC_ = 0.0;
            }
            else if (_SOC_ > 1.0) 
            {
                _SOC_ = 1.0;
            }
            
            const arma::colvec & parameters_OCV = parameters[2 * K + 2];
            
            const double SOC_2 = _SOC_ * _SOC_;
            return parameters_OCV[1] + 2.0 * parameters_OCV[2] * _SOC_ + 3.0 * parameters_OCV[3] * SOC_2;
        }
        
        double h(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
        {
            double V_ = V_OC(X_[0], t) + I[t] * theta_[0];
            for (unsigned int k = 0; k < K; k++)
            {
                V_ += X_[k + 1];
            }
            
            return V_;
        }
        
        arma::colvec f(const unsigned int & t, const arma::colvec & X_, const arma::colvec & theta_) 
        {
            arma::colvec X_new = arma::zeros(X_.size());
            X_new[0] = X_[0] + eta * dt / Cap_max * I[t];
            
            for (unsigned int k = 0; k < K; k++) 
            {
                const double & Vk_ = std::exp(-dt / theta_[2 * k + 2]);
                X_new[k + 1] = X_[k + 1] * Vk_ + 
                    std::exp(theta_[1 + 2 * k]) * (1.0 + Vk_) * I[t];
            }
            
            return X_new;
        }
        
        void update_F(const unsigned int & t, const arma::colvec & theta_) 
        {
            F = arma::zeros(K + 1, K + 1);
            F(0, 0) = 1.0;
            for (unsigned int k = 0; k < K; k++) 
            {
                F(k + 1, k + 1) = std::exp(-dt / theta_[2 * k + 2]) + theta_[2 * k + 1] * I[t];
            }
        }
        
        void update_H_x() 
        {
            H_x = arma::ones(1, K + 1);
            H_x(0, 0) = V_OC_(X[0]);
        }
        
        // Updating theta_ 
        double double_exponential(const double & I_, const double & SOC_, const double & a, const double & b, const double & c, const double & d) 
        {
            return a * std::exp(-b * I_) + c + d * std::abs(SOC_ + 1.0);
        }
        
        void update_theta(const unsigned int & t, arma::mat & theta_) 
        {
            arma::mat parameters_i;
            
            // R0
            parameters_i = parameters[0];
            theta_(0, t) = double_exponential(std::abs(I[t]), X[0] / 100, parameters_i(0, 0), parameters_i(1, 0), parameters_i(2, 0), parameters_i(3, 0));
            
            // Rk and Ck
            for (unsigned int k = 0; k < K; k++)
            {
                parameters_i = parameters[2 * k + 1];
                theta_(2 * k + 1, t) = double_exponential(std::abs(I[t]), X[0] / 100, parameters_i(0, 0), parameters_i(1, 0), parameters_i(2, 0), parameters_i(3, 0));
                
                parameters_i = parameters[2 * k + 2];
                const double & normalised_Ck = double_exponential(std::abs(I[t]), X[0] / 100, parameters_i(0, 0), parameters_i(1, 0), parameters_i(2, 0), parameters_i(3, 0));
                theta_(2 * k + 2, t) = normalised_Ck; // std::pow(10.0, k) * 
            }
            
            // Cap
            parameters_i = parameters[2 * K + 1];
            Cap_max = 36 * double_exponential(std::abs(I[t]), 0.0, parameters_i(0, 0), parameters_i(1, 0), parameters_i(2, 0), 0.0);
            if (Cap_max < 36 * 1.0) 
            {
                Cap_max = 36 * 2.55;
            }
            else if (Cap_max > 36 * 4.0) 
            {
                Cap_max = 36 * 2.60;
            }
        }
        
    public:
        //
        arma::colvec X, VT, VOC, SOC, S_x;
        arma::mat F, H_x, Q_x, P_x, R_x, theta;
        
        std::vector<arma::mat> parameters;

        //, 
        EFKModel(const arma::colvec & I_, const arma::colvec & IF_, const arma::colvec & V_, 
                 const std::vector<arma::mat> & parameters_0, const double & dt_, const unsigned int & K_, const double & SOC_intial_, 
                 const bool & trace_, const unsigned int & traceLimit_) :
            I(I_), IF(IF_), V(V_), parameters(parameters_0), dt(dt_), K(K_), trace(trace_), traceLimit(traceLimit_), N(1)
        {
            pi = 2.0 * std::acos(0.0);

            Cap_max = 36 * 2.568633;
            eta = 1.0;

            T = I.size();

            VT = arma::colvec(T);
            // VOC = arma::colvec(T);
            SOC = arma::colvec(T);
            // S_x = arma::colvec(T);

            X = arma::zeros(K + 1);
            X[0] = SOC_intial_;

            Q_x = 1e-8 * arma::eye(K + 1, K + 1);
            P_x = 10.0 * arma::eye(K + 1, K + 1);
            R_x = 10.0 * arma::eye(N, N);
            
            theta = arma::mat(2 * K + 1, T);
        }
        
        // 
        void Filter() 
        {
            boost::mt19937 rng;
            boost::random::uniform_01<> uniform_real;
            boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> > generate_uniform_real(rng, uniform_real);
            
            double V_hat;
            arma::colvec X_;
            arma::mat P_x_, L_x, HL_x, S_x_;
            
            for (unsigned int t = 0; t < T; t++) 
            {
                if (trace & ((t == 0) | ((t % traceLimit) == 0) | (t == (T - 1))))
                {
                    Rcpp::Rcout << "Iteration: " << t + 1 << " / " << T << "\n";
                }
                
                update_theta(t, theta);
                update_F(t, theta.col(t));
                
                //// Time-update
                X_ = f(t, X, theta.col(t));
                P_x_ = F * P_x * F.t() + Q_x;
                
                //// Measurement-update
                // X
                update_H_x();
                S_x_ = H_x * P_x_ * H_x.t() + R_x;
                L_x = P_x_ * H_x.t() * arma::inv(S_x_);
                
                V_hat = h(t, X_, theta.col(t));
                X = X_ + L_x * (V[t] - V_hat);
                
                HL_x = arma::eye(K + 1, K + 1) - L_x * H_x;
                P_x = HL_x * P_x_ * HL_x.t() + L_x * R_x * L_x.t();
                
                // VOC[t] = V_OC(X[0]);
                VT[t] = h(t, X, theta.col(t));
                SOC[t] = X[0];
                // S_x[t] = S_x_(0, 0);
            }
        }
        
        arma::vec Predict();
};

//' @title Filtering of state-of-charge and health.
//' 
//' @description Extended Kalman filter for determining the latent state-of-charge and health state-equation. 
//' 
//' @param I A current profile.
//' @param V The observed voltage.
//' @param parameters_0 The initial parameter values.
//' @param SOC_0 The initial value of the SOC.
//' @param dt Step-size.
//' @param K Number of RC branches in the EEC.
//' @param trace TRUE/FALSE: Show trace?
//' @param traceLimit Used to limit trace.
//' 
//' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
//' @export
//[[Rcpp::export]] 
Rcpp::List SOEFK(const arma::colvec & I, const arma::colvec & IF, const arma::colvec & V,  const std::vector<arma::mat> & parameters_0, const double & SOC_0,
                 const double & dt, const unsigned int K, const bool & trace, const unsigned int & traceLimit) 
{
    EFKModel EKF(I, IF, V, parameters_0, dt, K, SOC_0, trace, traceLimit);
    EKF.Filter();
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = EKF.VT, 
                              Rcpp::Named("SOC") = EKF.SOC);
}


// Rcpp::Named("V_OC") = EKF.VOC,
//     Rcpp::Named("SOC") = EKF.SOC,
//     Rcpp::Named("Theta") = EKF.theta, 
//     Rcpp::Named("Sigma") = EKF.S_x, 
//     Rcpp::Named("P") = EKF.P_x, 
//     Rcpp::Named("Q") = EKF.Q_x, 
//     Rcpp::Named("R") = EKF.R_x