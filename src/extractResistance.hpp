#ifndef extractResistance_H
#define extractResistance_H

class ExtractResistance {
private:
    std::vector<double> I, V, T;
    unsigned int N;
    double epsilon, Q_max, eta;
    
    double sign(double val) 
    {
        return ((double(0) < val) - (val < double(0)));
    } 
    
    void update_SOC(const double & n, const double & dt) 
    {
        double SOC_new = SOC[n - 1] + dt * eta * (I[n] / (3600 * Q_max));
        if (SOC_new < 0.0) {
            SOC_new = 0.0;
        }
        else if (SOC_new > 1.0) {
            SOC_new = 1.0;
        }
        
        SOC[n] = SOC_new;
    }
    
public:
    std::vector<double> R, S, SOC, ID;
    
    ExtractResistance(const std::vector<double> & I_, const std::vector<double> & V_, 
                      const std::vector<double> & T_, const double & epsilon_, 
                      const double & Q_max_, const double & eta_, const double & SOC_0) :
    I(I_), V(V_), T(T_), N(I.size()), epsilon(epsilon_), Q_max(Q_max_), eta(eta_)
    {
        if (T_.size() == 0) 
        {
            T = std::vector<double>(N);
            for (unsigned int n = 0; n < N; n++) 
            {
                T[n] = n + 1;
            }
        }
        
        //
        R = std::vector<double>(N);
        S = std::vector<double>(N);
        
        ID = std::vector<double>(N);
        
        SOC = std::vector<double>(N);
        SOC[0] = SOC_0;
    };
    
    //
    void Extraction() 
    {
        R[0] = 0.0; 
        S[0] = -1.0;
        ID[0] = 0.0;
        
        double I_scale = 0.0, V_scale = 0.0, T_scale = 0.0, ID_scale = 0.0; 
        for (unsigned int n = 1; n < N; n++)
        {
            update_SOC(n, T[n] - T[n - 1]);
            
            //
            const double abs_I_n_1 = std::abs(I[n - 1]);
            const double abs_I_n = std::abs(I[n]);
            const double abs_delta_I_n = std::abs(I[n] - I[n - 1]);
            
            //
            if ((abs_delta_I_n < epsilon) & (abs_I_n > epsilon)) 
            {
                T_scale += T[n] - T[n - 1];
            }
            else if ((abs_delta_I_n > epsilon) & (abs_I_n_1 < epsilon) & (abs_I_n > epsilon)) 
            {
                T_scale = 0.0;
                ID_scale += 1.0;
                
                V_scale = V[n - 1];
                I_scale = I[n - 1];
            }
            else if ((abs_delta_I_n > epsilon) & (abs_I_n_1 > epsilon) & (abs_I_n > epsilon)) 
            {
                T_scale = 0.0;
                ID_scale += 1.0;
            }
            
            // 
            const double delta_V = V[n] - V_scale;
            const double delta_I = I[n] - I_scale;
            
            //
            if (std::abs(delta_I) < epsilon) 
            {
                R[n] = 0.0;
                S[n] = -1.0;
            }
            else 
            {
                R[n] = std::abs(delta_V / delta_I);
                S[n] = T_scale;
            }
            
            ID[n] = ID_scale;
        }
    }
};

#endif