#include <RcppArmadillo.h>

class ParameterEKF {
    private:
        //
        arma::colvec I, V, e;
        arma::mat A, H, C;
        double dt, Q_max, deltaI;  
        bool trace;
        unsigned int traceLimit, S, N, K, T, W;
        
        //
        void update_SOC(const double & t) 
        {
            const double SOC_ = SOC[t] - dt * (I[t] / (3600 * Q_max));
            
            if (SOC_ < 0.0)
            {
                SOC[t + 1] = 0.0;
            }
            else if (SOC_ > 1.0)
            {
                SOC[t + 1] = 1.0;
            }
            else {
                SOC[t + 1] = SOC_;
            }
            
            // SOC[t + 1] = SOC_;
        }
        
        // 
        double h(const arma::colvec & theta_) 
        {
            const double V_ = theta_[1];
            return V_;
        }
        
        double V_dot(const unsigned int & t, const arma::colvec & theta_, const double & V_) 
        {
            const double V_OC = 2.0 + theta_[0];
            const double omega_1 = theta_[3];
            const double R_1 = std::exp(theta_[4]);
            const double R_0 = std::exp(theta_[5]);
            
            double V_dot_;
            double deltaI_ = 0;
            if (t > 0) 
            {
                deltaI_ = I[t] - I[t - 1];
            }
            
            V_dot_ = V_OC * omega_1 - V_ * omega_1 - omega_1 * R_0 * I[t] - 
                omega_1 * R_1 * I[t] - R_0 * deltaI_;
            
            return V_dot_;
        }
        
        double V_k_dot(const unsigned int & t, const arma::colvec & theta_, const unsigned int & k) 
        {
            const double V_k = theta_[2 + k];
            const double omega_k = theta_[2 + K + k];
            const double R_k = std::exp(theta_[2 + 2 * K + k]);
            
            double V_n_dot_ = -omega_k * V_k + omega_k * R_k * I[t];
            return V_n_dot_;
        }
        
        arma::colvec f_theta(const unsigned int & t, const arma::colvec & theta_) 
        {
            const double V_ = theta_[1];
            
            arma::colvec theta_new = arma::zeros(S, 1);
            theta_new[1] = V_dot(t, theta_, V_); 
            
            for (unsigned int k = 0; k < K; k++) 
            {
                theta_new[2 + k] = V_k_dot(t, theta_, k);
            }
            
            return theta_new;
        }
        
        arma::colvec f(const unsigned int & t, const arma::colvec & theta_) 
        {
            arma::colvec theta_change = f_theta(t, theta_);
            return theta_ + dt * theta_change;
        }
        
        void update_A(const unsigned int & t, const arma::colvec & theta_) 
        {
            //
            const double V_OC = 2.0 + theta_[0];
            const double V_ = theta_[1];
            const double V_1 = theta_[2];
            const double omega_1 = theta_[3];
            const double R_1 = std::exp(theta_[4]);
            const double R_0 = std::exp(theta_[5]);
            
            //
            arma::mat A_ = arma::zeros(S, S);
            
            // 
            A_(1, 0) = omega_1;
            A_(1, 1) = -omega_1;
            A_(1, 3) = V_OC - V_1 - (R_0 + R_1) * I[t];
            A_(1, 4) = -omega_1 * deltaI * R_1;
            A_(1, 5) = (-omega_1 * I[t] - deltaI) * R_0; 
            
            //
            A_(2, 2) = -omega_1;
            A_(2, 3) = -V_1 + R_1 * I[t];
            A_(2, 4) = omega_1 * I[t] * R_1;
            
            //
            A = arma::eye(S, S) + dt * A_;
        }
        
        void update_H() 
        {
            H = arma::zeros(1, S);
            H(0, 1) = 1.0;
        }
        
        //
        void update_C(const unsigned int & t) 
        {
            double W_ = t + 1;
            unsigned int t_ = 0;
            if (t >= W)
            {
                t_ = t - W;
                W_ = W;
            }

            double C_ = 0.0;
            for (unsigned int w = (t + 1); w >= (t_ + 1); w--)
            {
                const double error_w = e[w - 1];
                const double error_outer = (error_w * error_w) / W_;
                
                C_ += error_outer; 
            }
            
            C = arma::mat(1, 1);
            C(0, 0) = C_;
        }
        
    public:
        //
        arma::colvec theta, VT, VT_, SOC;
        arma::mat Q, P, R, theta_trace;
        double LogLikelihood;
        
        //
        ParameterEKF(const arma::colvec & I_, const arma::colvec & V_, const arma::colvec & theta_,
                     const arma::mat & P_, const arma::mat & Q_, const arma::mat & R_,
                     const double & dt_, const unsigned int K_, const double SOC_0_, const double Q_max_,
                     const unsigned int W_, const bool & trace_, const unsigned int & traceLimit_) :
        I(I_), V(V_), theta(theta_), P(P_), Q(Q_), R(R_), dt(dt_), K(K_), Q_max(Q_max_), W(W_),
        trace(trace_), traceLimit(traceLimit_), N(1)
        {
            T = I.size();
            S = theta.size();
            
            SOC = arma::colvec(T + 1);
            SOC[0] = SOC_0_;
            
            VT = arma::colvec(T);
            theta_trace = arma::mat(S, T);
            e = arma::colvec(T);
        }
        
        // 
        void Filter() 
        {
            double VT_;
            arma::colvec theta_;
            arma::mat P_, L, LC, S_;
            
            update_H();
            for (unsigned int t = 0; t < T; t++) 
            {
                if (trace & ((t == 0) | (((t + 1) % traceLimit) == 0) | (t == (T - 1))))
                {
                    Rcpp::Rcout << "Iteration: " << t + 1 << " / " << T << "\n";
                }
                
                update_SOC(t);
                
                //// Time-update
                theta_ = f(t, theta);
                update_A(t, theta_);
                
                P_ = A * P * A.t() + Q;
                VT_ = h(theta_);
                
                //// Measurement-update
                const arma::mat HPH = H * P_ * H.t();
                S_ = HPH + R;
                const arma::mat S_1 = arma::inv(S_);
                
                L = (P_ * H.t()) * arma::inv(S_);
                LC = arma::eye(S, S) - L * H;
                
                theta = theta_ + L * (V[t] - VT_);
                P = LC * P_ * LC.t() + L * R * L.t();
                
                //// Logging
                VT[t] = h(theta);
                theta_trace.col(t) = theta;
                e[t] = VT[t] - VT_;
                    
                //// Adaptive coveriance update
                update_C(t);
                
                if (C(0, 0) < HPH(0, 0)) // 
                {
                    C(0, 0) = C(0, 0) + HPH(0, 0) + 1e-2;
                }
                
                Q = L * C * L.t();
                Q = Q + 1e-6 * arma::eye(Q.n_rows, Q.n_cols);
                
                R = C - HPH;
                R = 0.5 * R + 0.5 * R.t();
                R = R + 1e-6 * arma::eye(R.n_rows, R.n_cols);
                
                if (trace & ((t == 0) | (((t + 1) % traceLimit) == 0) | (t == (T - 1))))
                {
                    Rcpp::Rcout << "\tTheta:\n\t" << theta.t()
                                << "\tV = " << V[t] << " :: V_hat = " << VT[t] << " :: E = " << e[t] 
                                << " :: I = " << I[t] << " :: SOC: " << SOC[t + 1] << "\n"
                                << "-------------------------------\n";
                }
            }
        }
};

//[[Rcpp::export]] 
Rcpp::List ParameterEKFCpp(const arma::colvec & I, const arma::colvec & V, const arma::colvec & theta_0,
                           const double & SOC_0, const double & Q_max,
                           const arma::mat & P, const arma::mat & Q, const arma::mat & R,
                           const double & dt, const unsigned int K, const unsigned int W, 
                           const bool & trace, const unsigned int & traceLimit) 
{
    ParameterEKF PEKF(I, V, theta_0, P, Q, R, dt, K, SOC_0, Q_max, W, trace, traceLimit);
    PEKF.Filter();
    
    return Rcpp::List::create(Rcpp::Named("V") = V,
                              Rcpp::Named("V_hat") = PEKF.VT,
                              Rcpp::Named("SOC") = PEKF.SOC,
                              Rcpp::Named("Theta") = PEKF.theta, 
                              Rcpp::Named("ThetaTrace") = PEKF.theta_trace, 
                              Rcpp::Named("Q") = PEKF.Q,
                              Rcpp::Named("P") = PEKF.P,
                              Rcpp::Named("R") = PEKF.R);
}

