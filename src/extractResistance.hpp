#ifndef extractResistance_H
#define extractResistance_H

class ExtractResistance {
private:
    std::vector<double> I, V, T; //, ocv_pars;
    unsigned int N;
    double epsilon, Q_max, eta;
    
    double sign(double val) 
    {
        return ((double(0) < val) - (val < double(0)));
    } 
    
    void update_SOC(const unsigned int & n, const double & dt) 
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
    std::vector<double> R, S, SOC, ID, NonZero, Reset, V_scale_, I_scale_;
    std::vector<bool> ZeroUpdate;
    
    ExtractResistance(const std::vector<double> & I_, const std::vector<double> & V_, 
                      const std::vector<double> & T_, const double & epsilon_, 
                      const double & Q_max_, const double & eta_, const double & SOC_0) :
    I(I_), V(V_), T(T_), N(I.size()), epsilon(epsilon_), Q_max(Q_max_), eta(eta_)//, ocv_pars(ocv_pars_)
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
        V_scale_ = std::vector<double>(N);
        I_scale_ = std::vector<double>(N);
        
        
        ID = std::vector<double>(N);
        NonZero = std::vector<double>(N);
        Reset = std::vector<double>(N);
        ZeroUpdate = std::vector<bool>(N);
        
        SOC = std::vector<double>(N);
        SOC[0] = SOC_0;
    };
    
    //
    void Extraction() 
    {
        double I_scale = 0.0, V_scale = 0.0, T_scale = 0.0, ID_scale = 0.0, NonZero_scale = 0.0, 
            lambda = 1.0, SOC_scale = 0.0, S_scale = 0, S_scale_diff = 0; 
        bool zero_update = false;
        for (unsigned int n = 0; n < 2; n++) 
        {
            R[n] = 0.0; 
            S[n] = -1.0;
            ID[n] = 0.0;
            NonZero[n] = 0.0;
            Reset[n] = 0.0;
            ZeroUpdate[n] = zero_update;
            
            V_scale_[n] = V_scale;
            I_scale_[n] = I_scale;
        }
        
        for (unsigned int n = 2; n < N; n++)
        {
            update_SOC(n, T[n] - T[n - 1]);
            
            //
            const double abs_I_n_1 = std::abs(I[n - 1]);
            const double abs_I_n = std::abs(I[n]);
            const double abs_delta_I_n_1 = std::abs(I[n - 1] - I[n - 2]);
            const double abs_delta_I_n = std::abs(I[n] - I[n - 1]);
    
            //
            if ((abs_delta_I_n < epsilon) & (abs_I_n > epsilon)) 
            {
                S_scale_diff = n - S_scale;
                T_scale += T[n] - T[n - 1];
            }
            else if ((abs_delta_I_n > epsilon) & (abs_delta_I_n_1 < epsilon)) //(abs_I_n_1 < epsilon) & (abs_I_n > epsilon)) 
            {
                lambda = 1.0;
                SOC_scale = SOC[n];
                S_scale = static_cast<double>(n);
                T_scale = 0.0;
                S_scale_diff = 0;
                
                ID_scale += 1.0;
                NonZero_scale = 0.0;
                
                V_scale = V[n - 1];
                I_scale = I[n - 1];
                
                if (abs_I_n_1 < epsilon) {
                    zero_update = true;
                }
                else {
                    zero_update = false;
                }
            }
            else if ((abs_delta_I_n > epsilon) & (abs_delta_I_n_1 > epsilon)) // (abs_I_n_1 > epsilon) & (abs_I_n > epsilon)) 
            {
                T_scale = 0.0;
                
                ID_scale += 1.0;
                NonZero_scale = 1.0;
            }
            
            // 
            const double delta_V = V[n] - V_scale;
            const double delta_I = I[n] - I_scale;
            
            //
            if (std::abs(delta_I) < epsilon) 
            {
                R[n] = 0.0;
                S[n] = -1.0;
                
                NonZero[n] = -1.0;
            }
            else 
            {
                R[n] = std::abs(delta_V / delta_I);
                S[n] = T_scale;
                
                NonZero[n] = NonZero_scale;
            }
            
            ID[n] = ID_scale;
            Reset[n] = S_scale_diff;
            ZeroUpdate[n] = zero_update;
            
            V_scale_[n] = V_scale;
            I_scale_[n] = I_scale;
        }
    }
};

#endif