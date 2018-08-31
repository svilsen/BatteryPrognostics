#include <RcppArmadillo.h>

class LookupModel {
    private :
        arma::mat Cap, R0, OCV;
        arma::colvec V, I, IF, IC, SOCList, IList;
        
        std::vector<arma::mat> Rk, Ck;
        
        int unsigned N, K, N_SOC, N_I, traceLimit;
        double dt;
        bool trace;
        
        double VK(const double & t_n, const double & Vk_n, const double & dt, const double & I_n, const double & Rk_n, const double & Ck_n) 
        {
            const double tau_n = std::exp(std::log(Rk_n) + std::log(Ck_n));
            const double Vk_n_1 = Vk_n + dt * (I_n / Ck_n - Vk_n / tau_n);
            return Vk_n_1;
        }
        
        double SOC(const double & t_n, const double & SOC_n, const double & dt, const double & I_n, const double & Cap_n) 
        {
            const double & SOC_n_1 = SOC_n + dt * (I_n / (3600 * Cap_n));
            return SOC_n_1;
        }
        
        double VT(const double & OCV_n, const double & I_n, const double & R0_n, const std::vector<double> & Vk_n) 
        {
            double VT_n = OCV_n + I_n * R0_n;
            for (unsigned int k = 0; k < K; k++) 
            {
                VT_n += Vk_n[k];
            }
            
            return VT_n;
        }
        
        unsigned int updateIndex(const arma::colvec & M, const double & v) 
        {
            const unsigned int & N_M = M.n_rows;
            unsigned int i_ = 0;
            
            double M_i_ = std::floor(M[i_] * 20) / 20;
            const double v_i_ = std::floor(v * 20) / 20;
            while ((i_ < (N_M - 1)) & (M_i_ < v_i_)) 
            {
                i_++;
                M_i_ = std::floor(M[i_] * 20) / 20;
            }

            return i_;
        }
        
    public : 
        arma::colvec T_hat, V_hat, SOC_hat; 
        
        LookupModel(const arma::mat & Cap_, const arma::mat & R0_, const std::vector<arma::mat> & Rk_,
                    const std::vector<arma::mat> & Ck_, const arma::mat & OCV_, const arma::colvec & V_, 
                    const arma::colvec & I_, const arma::colvec & IF_, const arma::colvec & IC_, 
                    const double & dt_, const double & SOC_, const arma::colvec & SOCList_, const arma::colvec IList_, 
                    const bool & trace_, const unsigned int & traceLimit_) : 
        Cap(Cap_), R0(R0_), Rk(Rk_), Ck(Ck_), K(Rk.size()), OCV(OCV_), V(V_), I(I_), IF(IF_), IC(IC_), N(I.size()), dt(dt_), 
        SOCList(SOCList_), IList(IList_), N_SOC(SOCList.size()), N_I(IList.size()), trace(trace_), traceLimit(traceLimit_)
        { 
            T_hat = arma::zeros(N / dt);
            V_hat = arma::zeros(N / dt);
            SOC_hat = arma::zeros(N / dt);
            SOC_hat(0) = SOC_;
        };
        
        void SimulateModel() 
        {
            unsigned int i_c = 0, i_f = 0, j = 0;
            
            std::vector<double> Vk_n(K, 0.0);
            V_hat[0] = 0.0;
            T_hat[0] = 0.0;
            for (unsigned int n = 0; n < N / dt - 1; n++) 
            {
                if (trace & ((n == 0) | (n == (N / dt - 2)) | (((n + 1) % traceLimit) == 0))) {
                    Rcpp::Rcout << "Iteration: " << n + 1 << " / " << (N / dt - 1) << "\n";
                }
                
                T_hat[n + 1] = T_hat[n] + dt; 
                
                const unsigned int & s = std::floor(dt * n);
                const double & I_n = I[s];
                const double & IC_n = IC[s];
                const double & IF_n = IF[s];
                
                i_f = updateIndex(IList, IF_n);
                i_c = updateIndex(IList, IC_n);
                
                // 
                SOC_hat[n + 1] = SOC(T_hat[n + 1], SOC_hat[n], dt, I_n, Cap[i_f]);
                j = updateIndex(SOCList, SOC_hat[n + 1]);
                
                const double & OCV_n = OCV(j, i_f);
                const double & R0_n = R0(j, i_c);
                
                for (std::size_t k = 0; k < K; k++) 
                {
                    const double & Rk_n = Rk[k](j, i_c);
                    const double & Ck_n = Ck[k](j, i_c); // / 25;        
                    Vk_n[k] = VK(T_hat[n + 1], Vk_n[k], dt, I_n, Rk_n, Ck_n);
                }
                
                V_hat[n + 1] = VT(OCV_n, I_n, R0_n, Vk_n);
            }
        }
};

//[[Rcpp::export]] 
Rcpp::List RCKCpp(const arma::colvec & I, const arma::colvec & IC, const arma::colvec & IF, const arma::colvec & V,
                  const arma::mat & R0, const std::vector<arma::mat> & Rk, const std::vector<arma::mat> & Ck,
                  const arma::mat & Cap, const arma::mat & OCV, const arma::colvec SOCList, const arma::colvec IList, 
                  const double & dt, const double & SOCStart, const bool & trace, const unsigned int & traceLimit) 
{
    LookupModel LM(Cap, R0, Rk, Ck, OCV, V, I, IF, IC, dt, SOCStart, SOCList, IList, trace, traceLimit);
    LM.SimulateModel();
    
    const arma::colvec & V_Error = V - LM.V_hat;
    return Rcpp::List::create(Rcpp::Named("V") = V, 
                              Rcpp::Named("V_hat") = LM.V_hat, 
                              Rcpp::Named("Residuals") = V_Error, 
                              Rcpp::Named("SOC") = LM.SOC_hat, 
                              Rcpp::Named("T") = LM.T_hat);
}
