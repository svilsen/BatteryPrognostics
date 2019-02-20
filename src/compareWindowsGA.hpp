#ifndef compareWindowsGA_H
#define compareWindowsGA_H

class CompareWindowsGA {
private:
    std::vector<double> I1, I2, SOC1, SOC2, epsilon;
    unsigned int W, R, T1, T2, RT1, RT2, trace_limit;
    double delta;
    bool trace;
    
public:
    std::vector<std::vector<double> > WindowError;
    std::vector<std::vector<double> > I1Start;
    std::vector<unsigned int> RI1, RI2;
    
    CompareWindowsGA(const std::vector<double> & I1_, const std::vector<double> & I2_,
                     const std::vector<double> & SOC1_, const std::vector<double> & SOC2_,
                     const unsigned int & W_, const unsigned int & R_,
                     const std::vector<double> epsilon_, const double & delta_,
                     const bool & trace_, const unsigned int & trace_limit_) : 
        I1(I1_), I2(I2_), SOC1(SOC1_), SOC2(SOC2_), W(W_), R(R_), 
        epsilon(epsilon_), delta(delta_), 
        trace(trace_), trace_limit(trace_limit_)
    { 
        T1 = I1.size();
        T2 = I2.size();
    }
    
    
};

#endif