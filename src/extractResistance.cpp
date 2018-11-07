#include <Rcpp.h>
#include "extractResistance.hpp"

//[[Rcpp::export()]]
Rcpp::List extract_resistance_cpp(const std::vector<double> & I, const std::vector<double> & V,
                                  const std::vector<double> & T_s, const double & epsilon, 
                                  const double & Q_max, const double & eta, const double & SOC_0) 
{
    ExtractResistance ER(I, V, T_s, epsilon, Q_max, eta, SOC_0);
    ER.Extraction();
    
    return Rcpp::List::create(Rcpp::Named("R") = ER.R,
                              Rcpp::Named("S") = ER.S, 
                              Rcpp::Named("SOC") = ER.SOC, 
                              Rcpp::Named("ID") = ER.ID);
}