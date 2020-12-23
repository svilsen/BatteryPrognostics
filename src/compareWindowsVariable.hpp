#ifndef compareWindows_H
#define compareWindows_H

class CompareWindowsVariable {
private:
    std::vector<double> I1, I2, V1, V2;
    unsigned int Wmin, Wmax, Wmid, T1, T2, RT1, RT2, trace_limit;
    std::vector<unsigned int> RI1, RI2;
    double epsilon, delta;
    bool trace;
    
    unsigned int ZeroCount(const unsigned int & t, const std::vector<double> & I) 
    {
        unsigned int w = 0;
        unsigned int sum_zeroes = 0;
        while ((sum_zeroes < 1) & (w < Wmax)) 
        {
            if (std::abs(I[t + w]) < delta) 
            {
                sum_zeroes++;
            } 
            
            w++;
        }
        
        return sum_zeroes;
    }
    
public:
    std::vector<double> WindowError;
    std::vector<int> S1, S2;
    std::vector<unsigned int> W1, W2;
    
    CompareWindowsVariable(const std::vector<double> & I1_, const std::vector<double> & I2_,
                           const std::vector<double> & V1_, const std::vector<double> & V2_,
                           const unsigned int & Wmin_, const unsigned int & Wmax_,
                           const double & epsilon_, const double & delta_,
                           const bool & trace_, const unsigned int & trace_limit_) : 
        I1(I1_), I2(I2_), V1(V1_), V2(V2_), T1(I1.size()), T2(I2.size()), 
        Wmin(Wmin_), Wmax(Wmax_), epsilon(epsilon_), delta(delta_),
        trace(trace_), trace_limit(trace_limit_) 
    { 
        Wmid = std::ceil((Wmax + Wmin) / 2);
    }
    
    void Reduce()
    {
        RT1 = 0;
        RT2 = 0;
        
        RI1 = std::vector<unsigned int>(T1 - Wmax);
        RI2 = std::vector<unsigned int>(T2 - Wmax);
        
        const unsigned int & T = std::max(T1, T2);
        for (unsigned int t = 0; t < (T - Wmax); t++)  
        {
            if (t < (T1 - Wmax)) 
            {
                const unsigned int zero_count = ZeroCount(t, I1); 
                if (zero_count < 1) 
                {
                    RI1[RT1] = t;
                    RT1++;
                }
            }
            
            if (t < (T2 - Wmax)) 
            {
                const unsigned int zero_count = ZeroCount(t, I2); 
                if (zero_count < 1) 
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
        WindowError = std::vector<double>(RT2, delta);
        S1 = std::vector<int>(RT2, -1); 
        S2 = std::vector<int>(RT2, -1);
        W1 = std::vector<unsigned int>(RT2, 0); 
        W2 = std::vector<unsigned int>(RT2, 0);
        
        for (unsigned int rt2 = 0; rt2 < RT2; rt2++) 
        {
            if ((rt2 == 0) || (((rt2 + 1) % trace_limit) == 0) || (rt2 == (RT2 - 1))) 
            {
                Rcpp::Rcout << "Iteration: " << rt2 + 1 << "\n";
            }
            
            const double & V2_1 = V2[RI2[rt2]];
            for (unsigned int rt1 = 0; rt1 < RT1; rt1++) 
            {
                const double & V1_1 = V1[RI1[rt1]];
                if (std::abs(V1_1 - V2_1) <= 1) 
                {
                    bool override = false;
                    for (unsigned int w = Wmin; w < Wmax; w++) 
                    {
                        const double & V2_w = V2[RI2[rt2] + w];
                        const double & V1_W = V1[RI1[rt1] + Wmid];
                        if (std::abs(V1_W - V2_w) <= 1) 
                        {
                            unsigned int u = 1;
                            double percentage_error = std::abs((I1[RI1[rt1]] - I2[RI2[rt2]]) / I1[RI1[rt1]]);
                            while ((u < w) & (percentage_error < epsilon)) 
                            {
                                double new_percentage_error = std::abs((I1[RI1[rt1] + u] - I2[RI2[rt2] + u]) / I1[RI1[rt1] + u]);
                                if (I1[RI1[rt1]] < 1e-6) 
                                {
                                    new_percentage_error = HUGE_VAL;
                                }
                                
                                if (new_percentage_error > percentage_error) 
                                {
                                    percentage_error = new_percentage_error;
                                }
                                
                                u++;
                            }
                            
                            if (percentage_error < epsilon)  
                            {
                                if (W2[rt2] < w) 
                                {
                                    WindowError[rt2] = percentage_error;
                                    S1[rt2] = RI1[rt1];
                                    S2[rt2] = RI2[rt2];
                                    W1[rt2] = Wmid;
                                    W1[rt2] = w;
                                }
                                else if (percentage_error < (WindowError[rt2] / 2)) 
                                {
                                    WindowError[rt2] = percentage_error;
                                    S1[rt2] = RI1[rt1];
                                    S2[rt2] = RI2[rt2];
                                    W1[rt2] = Wmid;
                                    W1[rt2] = w;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
};

#endif