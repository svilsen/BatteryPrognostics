#include <Rcpp.h>

#include "compareWindows.hpp"

//// Private functions
unsigned int CompareWindows::ZeroCount(const unsigned int & t, const std::vector<double> & I) 
{
    unsigned int w = 0;
    unsigned int sum_zeroes = 0;
    while ((sum_zeroes < R) & (w < W)) 
    {
        if (std::abs(I[t + w]) < epsilon[1]) 
        {
            sum_zeroes++;
        } 
        
        w++;
    }
    
    return sum_zeroes;
}


//// Constructors
CompareWindows::CompareWindows(const std::vector<double> & I1_, const std::vector<double> & I2_,
                               const std::vector<double> & V1_, const std::vector<double> & V2_,
                               const unsigned int & W_, const unsigned int & R_,
                               const std::vector<double> & epsilon_, const double & delta_,
                               const bool & trace_, const unsigned int & trace_limit_) : 
    I1(I1_), I2(I2_), T1(I1.size()), T2(I2.size()), V1(V1_), V2(V2_), W(W_), R(R_), 
    epsilon(epsilon_), delta(delta_), trace(trace_), 
    trace_limit(trace_limit_) { };

CompareWindows::CompareWindows(const std::vector<double> & I1_, const std::vector<double> & I2_,
                               const std::vector<double> & V1_, const std::vector<double> & V2_,
                               const std::vector<unsigned int> & RI2_, const std::vector<unsigned int> & W2_,
                               const unsigned int & W_, const unsigned int & R_,
                               const std::vector<double> & epsilon_, const double & delta_,
                               const bool & trace_, const unsigned int & trace_limit_) : 
    I1(I1_), I2(I2_), T1(I1.size()), T2(I2.size()), V1(V1_), V2(V2_), RI2(RI2_), W2(W2_),
    W(W_), R(R_), epsilon(epsilon_), delta(delta_), trace(trace_), 
    trace_limit(trace_limit_) 
    { 
        RT2 = RI2.size();
    }


//// Public functions
void CompareWindows::ReduceBoth()
{
    RT1 = 0;
    RT2 = 0;
    
    RI1 = std::vector<unsigned int>(T1 - W);
    RI2 = std::vector<unsigned int>(T2 - W);
    
    const unsigned int & T = std::max(T1, T2);
    for (unsigned int t = 0; t < (T - W); t++)  
    {
        if (t < (T1 - W)) 
        {
            const unsigned int zero_count = ZeroCount(t, I1); 
            if (zero_count < R) 
            {
                RI1[RT1] = t;
                RT1++;
            }
        }
        
        if (t < (T2 - W)) 
        {
            const unsigned int zero_count = ZeroCount(t, I2); 
            if (zero_count < R) 
            {
                RI2[RT2] = t;
                RT2++;
            }
        }
    }
    
    RI1.resize(RT1);
    RI1.shrink_to_fit();
    
    RI2.resize(RT2);
    RI2.shrink_to_fit();
    
    W2 = std::vector<unsigned int>(RT2, W);
    
    if (trace)
    {
        Rcpp::Rcout << "'I1' reduced from " << T1 << " to " << RT1 << "\n"
                    << "'I2' reduced from " << T2 << " to " << RT2 << "\n";
    }
}

void CompareWindows::Reduce()
{
    RT1 = 0;
    
    RI1 = std::vector<unsigned int>(T1 - W);
    for (unsigned int t = 0; t < (T1 - W); t++)  
    {
        const unsigned int zero_count = ZeroCount(t, I1); 
        if (zero_count < R) 
        {
            RI1[RT1] = t;
            RT1++;
        }
    }
    
    RI1.resize(RT1);
    RI1.shrink_to_fit();
    
    if (trace)
    {
        Rcpp::Rcout << "'I1' reduced from " << T1 << " to " << RT1 << "\n";
    }
}

void CompareWindows::Compare() 
{
    CurrentError = std::vector<double>(RT2, HUGE_VAL);
    VoltageError = std::vector<double>(RT2, HUGE_VAL);
    S1 = std::vector<int>(RT2, -1);
    
    for (unsigned int rt2 = 0; rt2 < RT2; rt2++) 
    {
        if (trace & ((rt2 == 0) || (((rt2 + 1) % trace_limit) == 0) || (rt2 == (RT2 - 1)))) 
        {
            Rcpp::Rcout << "Iteration: " << rt2 + 1 << "\n";
        }
        
        const double & V2_t2 = V2[RI2[rt2]];
        for (unsigned int rt1 = 0; rt1 < RT1; rt1++) 
        {
            const double & V1_t1 = V1[RI1[rt1]];
            const double voltage_error = std::abs(V2_t2 - V1_t1);
            if (voltage_error < delta) {
                double max_difference_current = 0.0;
                
                unsigned int w = 0;
                while ((w < W2[rt2]) & (max_difference_current < epsilon[0])) 
                {
                    const double difference_current = std::abs((I2[RI2[rt2] + w] - I1[RI1[rt1] + w]) / I1[RI1[rt1] + w]); //
                    if (difference_current > max_difference_current)
                    {
                        max_difference_current = difference_current;
                    }
                    
                    w++;
                }
                
                if (max_difference_current < epsilon[0]) 
                {
                    if ((voltage_error < VoltageError[rt2]) & (std::abs(voltage_error - VoltageError[rt2]) > 1e-8)) 
                    {
                        CurrentError[rt2] = max_difference_current;
                        VoltageError[rt2] = voltage_error;
                        S1[rt2] = RI1[rt1];
                    }
                    else if (std::abs(voltage_error - VoltageError[rt2]) < 1e-8)
                    {
                        if (max_difference_current < CurrentError[rt2]) 
                        {
                            CurrentError[rt2] = max_difference_current;
                            VoltageError[rt2] = voltage_error;
                            S1[rt2] = RI1[rt1];
                        }
                    }
                }
            }
        }
    }
}

// R wrapper functions
//[[Rcpp::export()]]
Rcpp::List compare_windows_raw_cpp(const std::vector<double> & I1,
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
    
    CW.ReduceBoth();
    CW.Compare();
    
    return Rcpp::List::create(Rcpp::Named("CurrentError") = CW.CurrentError, 
                              Rcpp::Named("VoltageError") = CW.VoltageError,
                              Rcpp::Named("S1") = CW.S1, 
                              Rcpp::Named("S2") = CW.RI2);
}

//[[Rcpp::export()]]
Rcpp::List compare_windows_single_cpp(const std::vector<double> & I1,
                                      const std::vector<double> & I2,
                                      const std::vector<double> & V1,
                                      const std::vector<double> & V2,
                                      const std::vector<unsigned int> & RI2,
                                      const std::vector<unsigned int> & W2,
                                      const unsigned int & W,
                                      const unsigned int & R,
                                      const std::vector<double> & epsilon,
                                      const double & delta,
                                      const bool & trace,
                                      const unsigned int & trace_limit) 
{
    CompareWindows CW(I1, I2, V1, V2, RI2, W2, W, R, epsilon, delta, trace, trace_limit); 
    
    CW.Reduce();
    CW.Compare();
    
    return Rcpp::List::create(Rcpp::Named("CurrentError") = CW.CurrentError, 
                              Rcpp::Named("VoltageError") = CW.VoltageError,
                              Rcpp::Named("S1") = CW.S1, 
                              Rcpp::Named("S2") = CW.RI2);
}
