#include <Rcpp.h>

#include "compareWindowsVariable.hpp"

//[[Rcpp::export()]]
Rcpp::List compare_windows_variable_cpp(const std::vector<double> & I1,
                                        const std::vector<double> & I2,
                                        const std::vector<double> & V1,
                                        const std::vector<double> & V2,
                                        const unsigned int & Wmin,
                                        const unsigned int & Wmax,
                                        const double & epsilon,
                                        const double & delta,
                                        const bool & trace,
                                        const unsigned int & trace_limit) 
{
    CompareWindowsVariable CWV(I1, I2, V1, V2, Wmin, Wmax, epsilon, delta, trace, trace_limit); 
    
    CWV.Reduce();
    CWV.Compare();
    
    return Rcpp::List::create(Rcpp::Named("WindowError") = CWV.WindowError, 
                              Rcpp::Named("S1") = CWV.S1,
                              Rcpp::Named("S2") = CWV.S2,
                              Rcpp::Named("W1") = CWV.W1,
                              Rcpp::Named("W2") = CWV.W2);
}