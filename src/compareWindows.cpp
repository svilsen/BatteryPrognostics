#include <Rcpp.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include "compareWindows.hpp"

//[[Rcpp::export()]]
Rcpp::List compare_windows_cpp(const std::vector<double> & I1, 
                               const std::vector<double> & I2, 
                               const std::vector<double> & V1,
                               const std::vector<double> & V2,
                               const unsigned int & W,
                               const unsigned int & R,
                               const std::vector<double> & epsilon,
                               const double & delta, 
                               const bool & trace, 
                               const unsigned int & trace_limit) 
{
    CompareWindows CW(I1, I2, V1, V2, W, R, epsilon, delta, trace, trace_limit); 
    
    CW.Reduce();
    CW.Compare();
    
    return Rcpp::List::create(Rcpp::Named("WindowError") = CW.WindowError, 
                              Rcpp::Named("I1StartIndices") = CW.I1Start, 
                              Rcpp::Named("I2StartIndices") = CW.RI2);
}