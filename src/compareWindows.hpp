#ifndef compareWindows_H
#define compareWindows_H

class CompareWindows {
private:
    // Objects
    std::vector<double> I1, I2, V1, V2, epsilon;
    unsigned int W, R, T1, T2, RT1, RT2, trace_limit;
    double delta;
    bool trace;
    
    // Functions
    unsigned int ZeroCount(const unsigned int & t, const std::vector<double> & I);
    
public:
    // Objects
    std::vector<double> CurrentError, VoltageError;
    std::vector<int> S1;
    std::vector<unsigned int> RI1, RI2, W1, W2;
    
    // Constructors
    CompareWindows(const std::vector<double> & I1_, const std::vector<double> & I2_, 
                   const std::vector<double> & V1_, const std::vector<double> & V2_, 
                   const unsigned int & W_, const unsigned int & R_, 
                   const std::vector<double> & epsilon_, const double & delta_, 
                   const bool & trace_, const unsigned int & trace_limit_);
    
    CompareWindows(const std::vector<double> & I1_, const std::vector<double> & I2_, 
                   const std::vector<double> & V1_, const std::vector<double> & V2_, 
                   const std::vector<unsigned int> & RI2_, const std::vector<unsigned int> & W2_,
                   const unsigned int & W_, const unsigned int & R_, 
                   const std::vector<double> & epsilon_, const double & delta_, 
                   const bool & trace_, const unsigned int & trace_limit_);
    
    // Functions
    void ReduceBoth();
    void Reduce();
    void Compare();
};

#endif