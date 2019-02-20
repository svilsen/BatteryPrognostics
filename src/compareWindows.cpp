#include <Rcpp.h>
#include "compareWindows.hpp"
#include "compareWindowsGA.hpp"

//[[Rcpp::export()]]
Rcpp::List compare_windows_cpp(const std::vector<double> & I1, 
                               const std::vector<double> & I2, 
                               const std::vector<double> & SOC1,
                               const std::vector<double> & SOC2,
                               const unsigned int & W,
                               const unsigned int & R,
                               const std::vector<double> & epsilon,
                               const double & delta, 
                               const bool & trace, 
                               const unsigned int & trace_limit) 
{
    CompareWindows CW(I1, I2, SOC1, SOC2, W, R, epsilon, delta, trace, trace_limit); 
    
    CW.Reduce();
    CW.Compare();
    
    return Rcpp::List::create(Rcpp::Named("WindowError") = CW.WindowError, 
                              Rcpp::Named("I1StartIndices") = CW.I1Start, 
                              Rcpp::Named("I2StartIndices") = CW.RI2);
}