#ifndef compareWindows_H
#define compareWindows_H

class CompareWindows {
private:
    std::vector<double> I1, I2, SOC1, SOC2, epsilon;
    unsigned int W, R, T1, T2, RT1, RT2, trace_limit;
    double delta;
    bool trace;
    
    unsigned int ZeroCount(const unsigned int & t, const std::vector<double> & I) 
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
    
public:
    std::vector<std::vector<double> > WindowError;
    std::vector<std::vector<double> > I1Start;
    std::vector<unsigned int> RI1, RI2;
    
    CompareWindows(const std::vector<double> & I1_, const std::vector<double> & I2_, 
                   const std::vector<double> & SOC1_, const std::vector<double> & SOC2_, 
                   const unsigned int & W_, const unsigned int & R_, 
                   const std::vector<double> epsilon_, const double & delta_, 
                   const bool & trace_, const unsigned int & trace_limit_) : 
        I1(I1_), I2(I2_), SOC1(SOC1_), SOC2(SOC2_), W(W_), R(R_), epsilon(epsilon_), delta(delta_), 
        trace(trace_), trace_limit(trace_limit_)
    { 
        T1 = I1.size();
        T2 = I2.size();
    }
    
    
    void Reduce()
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
        
        if (trace)
        {
            Rcpp::Rcout << "'I1' reduced from " << T1 << " to " << RT1 << "\n"
                        << "'I2' reduced from " << T2 << " to " << RT2 << "\n";
        }
    }
    
    void Compare() 
    {
        WindowError = std::vector<std::vector<double>>(RT2);
        I1Start = std::vector<std::vector<double>>(RT2);
        
        for (unsigned int rt2 = 0; rt2 < RT2; rt2++) 
        {
            if ((rt2 == 0) || (((rt2 + 1) % trace_limit) == 0) || (rt2 == (RT2 - 1))) 
            {
                Rcpp::Rcout << "Iteration: " << rt2 + 1 << "\n";
            }
            
            unsigned int i = 0;
            const double & SOC2_t2 = SOC2[RI2[rt2]];
            
            std::vector<double> Windows_t2;
            std::vector<double> Indices_t2;
            
            for (unsigned int rt1 = 0; rt1 < RT1; rt1++) 
            {
                const double & SOC1_t1 = SOC1[RI1[rt1]];
                if (std::abs(SOC2_t2 - SOC1_t1) < delta) {
                    double window_difference_tolerance = 0.0;
                    
                    unsigned int w = 0;
                    bool dont_stop = true;
                    while (dont_stop) 
                    {
                        const bool difference_tolerance = std::abs((I2[RI2[rt2] + w] - I1[RI1[rt1] + w]) / I1[RI1[rt1] + w]) > epsilon[0]; //
                        if (difference_tolerance)
                        {
                            window_difference_tolerance++;
                        }
                        
                        dont_stop = (w < W) & (window_difference_tolerance < R);
                        w++;
                    }
                    
                    if (window_difference_tolerance < R) {
                        Windows_t2.push_back(window_difference_tolerance);
                        Indices_t2.push_back(static_cast<double>(RI1[rt1]));
                        i++;
                    }
                }
            }
            
            WindowError[rt2] = Windows_t2;
            I1Start[rt2] = Indices_t2;
        }
    }
};

#endif