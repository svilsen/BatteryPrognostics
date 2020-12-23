#include <RcppArmadillo.h>

//' @title Corrects time
//' @description Correction of time in the case of inconsistant time-step sizes.
//' 
//' @param I A vector containing the current.
//' @param V A vector containing the voltage.
//' @param T A vector containing the time.
//' 
//' @return A matrix containing the corrected time, and interpolated current and voltage.
//' @export
//[[Rcpp::export()]]
arma::mat correct_time(const arma::colvec & I, const arma::colvec & V, const arma::colvec & T) 
{
    unsigned int S_total = 1;   
    unsigned int S = I.size();
    for (unsigned int s = 1; s < S; s++)
    {
        double new_time = std::ceil(T[s] - T[s - 1]);
        if (new_time < 1) 
            new_time = 1;
        
        unsigned int new_time_ = static_cast<unsigned int>(new_time);
        S_total += new_time_;
    }
    
    arma::mat M(S_total, 3);
    M(0, 0) = I[0];
    M(0, 1) = V[0];
    M(0, 2) = 1;
    
    unsigned int S_total_2 = 1;   
    for (unsigned int s = 1; s < S; s++) 
    {
        double new_time = std::ceil(T[s] - T[s - 1]);
        if (new_time < 1) 
            new_time = 1;
        
        unsigned int S_s = static_cast<unsigned int>(new_time);
        for (unsigned int s_s = 0; s_s < S_s; s_s++) 
        {
            M(S_total_2 + s_s, 0) = I[s];
            M(S_total_2 + s_s, 1) = V[s];
            M(S_total_2 + s_s, 2) = S_total_2 + s_s + 1;
        }
        
        S_total_2 += S_s;
    }
    
    return M;
}

