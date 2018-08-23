#include <RcppArmadillo.h>

#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

class StateSpaceModel {
    private:
        // 
        arma::colvec I, V;
        arma::mat F, G, H, Q;
        
        double R;
        unsigned int N;
        
        bool trace;
        
        // 
        double OCV(const double & SOC) 
        {
            return 10.56 * std::pow(SOC, 5.0) - 32.11 * std::pow(SOC, 4.0) + 36.30 * std::pow(SOC, 3.0) - 18.54 * std::pow(SOC, 2.0) + 
                4.91 * std::pow(SOC, 5.0) + 3.03;
        }
        
        void newF(const unsigned int & n) 
        {
            F = arma::eye(3, 3);
            F(1, 1) = std::exp(-1.0 / (R1[n] * C1[n]));
            F(2, 2) = std::exp(-1.0 / (R2[n] * C2[n]));
        }
        
        void newG(const unsigned int & n) 
        {
            G = arma::ones(3, 1);
            G(0, 0) = -1.0 / C[n];
            G(1, 0) = R1[n] * (1.0 - std::exp(-1.0 / (R1[n] * C1[n])));
            G(2, 0) = R2[n] * (1.0 - std::exp(-1.0 / (R2[n] * C2[n])));
        }
        
        double newR0(const unsigned int & n) 
        {
            const double normaliser = 0.008167463;
            const double a = 1.0; // 791.41353;
            const double b = 1.00559;
            return normaliser * (b * std::exp(-std::abs(I[n]) / a));
        }
        
        double newR1(const unsigned int & n) 
        {
            const double normaliser = 0.004252113;
            const double a = 13.78734;
            const double b = 1.33075;
            return normaliser * (b * std::exp(-std::abs(I[n]) / a));
        }
        
        double newR2(const unsigned int & n) 
        {
            const double normaliser = 0.006004791;
            const double a = 10.66015;
            const double b = 1.42979;
            return normaliser * (b * std::exp(-std::abs(I[n]) / a));
        }
        
        double newC(const unsigned int & n) 
        {
            const double a = 2.575;
            const double b = 0.09273;
            const double c = 121.7;
            return a * std::exp(-std::pow((I[n] - b) / c, 2.0));
        }
        
        double newC1(const unsigned int & n) 
        {
            const double normaliser = 707185.6;
            const double a = 0.9997;
            const double b = 4.5971;
            return normaliser * (b * std::exp(-std::abs(I[n]) / a));
        }
        
        double newC2(const unsigned int & n) 
        {
            const double normaliser = 820723.9;
            const double a = 0.9044;
            const double b = 4.7572;
            return normaliser * (b * std::exp(-std::abs(I[n]) / a));
        }
        
        void newParameters(const unsigned int & n)
        {
            C[n] = newC(n);
            R0[n] = newR0(n);
            
            R1[n] = newR1(n);
            C1[n] = newC1(n);
            
            R2[n] = newR2(n);
            C2[n] = newC2(n);
        }
        
    public:
        arma::mat X, P;
        arma::colvec V_T, C, R0, C1, R1, C2, R2;
        
        //
        StateSpaceModel(const arma::vec & I_, const arma::vec & V_, const bool & trace_) : 
            I(I_), V(V_), N(I.size()), trace(trace_)
        {
            arma::mat R0 = arma::zeros(N + 1), R1 = arma::zeros(N + 1), R2 = arma::zeros(N + 1), 
                C = arma::zeros(N + 1), C1 = arma::zeros(N + 1), C2 = arma::zeros(N + 1);
            
            newParameters(0);
            newF(0);
            newG(0);
            
            H = arma::ones(1, 3);
            H(0, 1) = H(0, 2) = -1.0;
            
            Q = arma::eye(3, 3);
            R = 1.0;
                
            X = arma::zeros(N, 3);
            P = arma::zeros(3, 3);
            V_T = arma::zeros(N);
        }
        
        // 
        void Filter() 
        {
            boost::mt19937 rng;
            boost::random::uniform_01<> uniform_real;
            boost::variate_generator<boost::mt19937 &, boost::random::uniform_01<> > generate_uniform_real(rng, uniform_real);
            
            arma::colvec V_step, X_step;
            arma::mat K, P_step;
            double S = 1.0;
            
            arma::mat X_est = arma::zeros(3, 1);
            X_est(0, 0) = generate_uniform_real();
            X_est(1, 0) = generate_uniform_real();
            X_est(2, 0) = generate_uniform_real();
            
            for (unsigned int n = 0; n < N; n++) 
            {
                if (trace)
                {
                    Rcpp::Rcout << "Iteration: " << n + 1 << " / " << N << "\n";
                }
                
                //     
                newParameters(n);
                
                // 
                X_step = F * X_est + G * I[n];
                P_step = F * P * F.t() + Q;
                
                H(0, 0) = OCV(X_step(0, 0)) / X_step(0, 0);
                
                S = (H * P_step * H.t()).eval()(0, 0) + R;
                K = P_step * H.t() / S;
                
                V_step =  H * X_step - R0[n] * I[n];
                
                //
                X_est = X_step + (V[n] - V_step).eval()(0, 0) * K;
                P = P_step - K * H * P_step;
                
                newF(n);
                newG(n);
                
                V_T[n] = OCV(X_est(0, 0)) - X_est(1, 0) - R0[n] * I[n];
                X.row(n) = X_est.t();
            }
        }
        
        arma::vec Predict();
};


//' @title State-of-charge filtering
//' 
//' @description Kalman filter for determining the latent state-of-charge state-equation. 
//' 
//' @param I A current pR0file.
//' @param V The observed voltage.
//' @param trace TRUE/FALSE: Show trace?
//' 
//' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
//' @export
//[[Rcpp::export]] 
Rcpp::List SOC(const arma::colvec & I, const arma::colvec & V, const bool & trace) 
{
    StateSpaceModel SSM(I, V, trace);
    SSM.Filter();
    
    arma::mat X = SSM.X;
    return Rcpp::List::create(Rcpp::Named("V") = V, 
                              Rcpp::Named("V_hat") = SSM.V_T, 
                              Rcpp::Named("SOC") = X.col(0), 
                              Rcpp::Named("V_p") = X.col(1), 
                              Rcpp::Named("P") = SSM.P, 
                              Rcpp::Named("C") = SSM.C, 
                              Rcpp::Named("C1") = SSM.C1, 
                              Rcpp::Named("C2") = SSM.C2, 
                              Rcpp::Named("R0") = SSM.R0, 
                              Rcpp::Named("R1") = SSM.R1,
                              Rcpp::Named("R2") = SSM.R2);
}
