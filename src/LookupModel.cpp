#include <RcppArmadillo.h>

class LookupModel {
    private :
        arma::mat Cap, R0, R1, C1, OCV;
        arma::colvec V, I, IF, IC, SOCList, IList;
        
        int unsigned N, N_SOC, N_I;
        double dt;
        
        double V1(const double & t_n, const double & V1_n, const double & dt, const double & I_n, const double & R1_n, const double & C1_n) 
        {
            const double tau_n = std::exp(std::log(R1_n) + std::log(C1_n));
            const double V1_n_1 = V1_n + dt * (I_n / C1_n - V1_n / tau_n);
            return V1_n_1;
        }
        
        double SOC(const double & t_n, const double & SOC_n, const double & dt, const double & I_n, const double & Cap_n) 
        {
            const double & SOC_n_1 = SOC_n + dt * (I_n / (3600 * Cap_n));
            return SOC_n_1;
        }
        
        double VT(const double & OCV_n, const double & I_n, const double & R0_n, const double & V1_n) 
        {
            return (OCV_n + I_n * R0_n + V1_n);
        }
        
        unsigned int updateIndex(const unsigned int & n, const arma::colvec & M, const double & v) 
        {
            const unsigned int & N_M = M.size();
            unsigned int i_ = 0;
            while ((i_ < N_M) & (M[i_] < v)) 
            {
                i_++;
            }
            
            return i_;
        }
        
    public : 
        arma::colvec T_hat, V_hat, SOC_hat; 
        
        LookupModel(const arma::mat & Cap_, const arma::mat & R0_, const arma::mat & R1_, const arma::mat & C1_, const arma::mat & OCV_, 
                    const arma::colvec & V_, const arma::colvec & I_, const arma::colvec & IF_, const arma::colvec & IC_, 
                    const double & dt_, const double & SOC_, const arma::colvec & SOCList_, const arma::colvec IList_) : 
        Cap(Cap_), R0(R0_), R1(R1_), C1(C1_), OCV(OCV_), V(V_), I(I_), IF(IF_), IC(IC_), N(I.size()), dt(dt_), 
        SOCList(SOCList_), IList(IList_), N_SOC(SOCList.size()), N_I(IList.size())
        { 
            T_hat = arma::zeros(N / dt);
            V_hat = arma::zeros(N / dt);
            SOC_hat = arma::zeros(N / dt);
            SOC_hat(0) = SOC_;
        };
        
        void SimulateModel() 
        {
            unsigned int i = 0, i_c = 0, i_f = 0, j = 0;
            
            double V1_n = 0.0;
            V_hat[0] = 0.0;
            T_hat[0] = 0.0;
            for (unsigned int n = 0; n < N / dt - 1; n++) 
            {
                T_hat[n + 1] = T_hat[n] + dt; 
                
                const unsigned int & s = std::floor(dt * n);
                const double & I_n = I[s];
                const double & IC_n = IC[s];
                const double & IF_n = IF[s];
                
                i_f = updateIndex(n, IList, IF_n);
                i_c = updateIndex(n, IList, IC_n);
                
                // 
                SOC_hat[n + 1] = SOC(T_hat[n + 1], SOC_hat[n], dt, I_n, Cap[i_f]);
                j = updateIndex(n, SOCList, SOC_hat[n + 1]);
                
                const double & OCV_n = OCV(j, i_f);
                const double & R0_n = R0(j, i_c);
                
                const double & R1_n = R1(j, i_c);
                const double & C1_n = C1(j, i_c) / 50;
                        
                V1_n = V1(T_hat[n + 1], V1_n, dt, I_n, R1_n, C1_n);
                V_hat[n + 1] = VT(OCV_n, I_n, R0_n, V1_n);
            }
        }
};

//' @title Simulate RC-2 model
//' 
//' @description Simulate the RC-2 model (double polarisation) given look-up parameter tables.
//' 
//' @param I A current profile.
//' @param IC A current profile reflecting the change in current (not the measured current).
//' @param IF A current profile fixed in regions of zero current using the most recent non-zero current.
//' @param V The observed voltage.
//' @param R0 Look-up table of 'R0' values.
//' @param R1 Look-up table of 'R1' values.
//' @param C1 Look-up table of 'C1' values.
//' @param Cap Look-up table of 'Cap' values.
//' @param OCV Look-up table of 'OCV' values.
//' @param SOCList A vector of possible SOC values in the look-up tables.
//' @param IList A vector of possible current values in the look-up tables.
//' @param dt The simulation step-size.
//' @param SOCStart The starting SOC value.
//' @param trace TRUE/FALSE: Show trace?
//' 
//' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
//' @export
//[[Rcpp::export]] 
Rcpp::List RC2(const arma::colvec & I, const arma::colvec & IC, const arma::colvec & IF, const arma::colvec & V, 
               const arma::mat & R0, const arma::mat & R1, const arma::mat & C1, const arma::mat & Cap, const arma::mat & OCV, 
               const arma::colvec SOCList, const arma::colvec IList, const double & dt, const double & SOCStart, 
               const bool & trace) 
{
    LookupModel LM(Cap, R0, R1, C1, OCV, V, I, IF, IC, dt, SOCStart, SOCList, IList);
    LM.SimulateModel();
    
    const arma::colvec & V_Error = V - LM.V_hat;
    return Rcpp::List::create(Rcpp::Named("V") = V, 
                              Rcpp::Named("V_hat") = LM.V_hat, 
                              Rcpp::Named("Residuals") = V_Error, 
                              Rcpp::Named("SOC") = LM.SOC_hat, 
                              Rcpp::Named("T") = LM.T_hat);
}
