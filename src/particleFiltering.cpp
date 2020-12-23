#include <RcppArmadillo.h>
#include "particleFiltering.hpp"

// Sampling from multivariate normal distribution with armadillo
arma::mat mvrnorm(const int & n, const arma::vec & mu, const arma::mat & L) {
    const int & n_cols = mu.n_elem;
    
    arma::mat Z = arma::randn(n_cols, n);
    return arma::repmat(mu, 1, n) + L.t() * Z;
}

// The pdf of a univariate normal distribution
double dnorm(const double & y, const double & mu, const double & sigma) {
    double log_density = -0.5 * (std::log(2.0) + std::log(PI)) - std::log(sigma) - (y - mu) * (y - mu) / (2 * sigma * sigma);
    return std::exp(log_density);
}

arma::colvec mean(const arma::mat & X, const unsigned int & i) 
{
    arma::colvec x;
    unsigned int N, M;
    if (i == 1) 
    {
        N = X.n_rows;
        M = X.n_cols;
        x = arma::colvec(M, arma::fill::zeros);
        for (unsigned int n = 0; n < N; n++) 
        {
            x += X.row(n);
        }
    }
    else if (i == 2) 
    {
        M = X.n_rows;
        N = X.n_cols;
        x = arma::colvec(M, arma::fill::zeros);
        for (unsigned int n = 0; n < N; n++) 
        {
            x += X.col(n);
        }
    }
    
    return x / N;
}

// The parameter update functions
void ParticleFiltering::f_OCV(const double & SOC_t) 
{
    bool charging = latest_non_zero >= 0.0;
    
    double OCV_ = 0.0;
    if (charging) 
    {
        OCV_ = OCV_c_pars[0] + OCV_c_pars[1] * std::log(SOC_t) + OCV_c_pars[2] * std::log(1.0 - SOC_t);
    }
    else 
    {
        OCV_ = OCV_d_pars[0] + OCV_d_pars[1] * std::log(SOC_t) + OCV_d_pars[2] * std::log(1.0 - SOC_t);
    }
    
    OCV = std::exp(OCV_);
}

// void ParticleFiltering::f_R0(const double & I_t, const double & SOC_t) 
// {
//     double R0_ = R0_pars[0] + R0_pars[1] * (std::log(SOC_t)) + R0_pars[2] * (std::log(1.0 - SOC_t)) + R0_pars[3] * std::abs(I_t);
//     R0 = std::exp(R0_);
// }

void ParticleFiltering::f_RT(const double & I_t, const double & SOC_t) 
{
    double RT_ = RT_pars[0] + RT_pars[1] * (std::log(SOC_t)) + RT_pars[2] * (std::log(1.0 - SOC_t)) + RT_pars[3] * std::abs(I_t);
    RT = std::exp(RT_);
    
    R0 = RT_pars[4] * RT;
    
    Ri = arma::colvec(K, arma::fill::zeros);
    for (unsigned int k = 0; k < K; k++) 
    {
        const double & rho_n = RT_pars[5 + k];
        Ri[k] = RT * rho_n;
    }
}

void ParticleFiltering::f_Q(const double & I_t, const double & Temp_t) 
{
    const double Q_internal = (Q_pars[1] - Q_pars[2] * I_t) / (Q_pars[3] * Temp_t);
    const double log_Q = Q_pars[0] - Q_internal;
    Q = std::exp(log_Q) + Q_pars[4];
}

// Update function
double ParticleFiltering::soc_transform(const double & SOC_t) 
{
    return atan(SOC_t) / PI + 0.5;
}

double ParticleFiltering::inverse_soc_transform(const double & SOC_t) 
{
    return tan(PI * (SOC_t - 0.5));
}

void ParticleFiltering::update_soc_parameters(const double & I_t, const double & SOC_t) 
{
    f_OCV(SOC_t);
    f_RT(I_t, SOC_t);
}

arma::colvec ParticleFiltering::f(const double & I_t, const arma::colvec & x_k_t_1) 
{
    arma::colvec mu_k_t(K + 1, arma::fill::zeros);
    mu_k_t[0] = soc_transform(x_k_t_1[0]) + dt * (I_t / (3600 * Q));
    mu_k_t[0] = inverse_soc_transform(mu_k_t[0]);
    
    for (unsigned int k = 0; k < K; k++) 
    {
        const double & x_k_t_1_ = x_k_t_1[k + 1];
        const double & Ri_n = Ri[k];
        const double & Ci_n = Ci[k];
        const double & tau_n = std::exp(std::log(Ri_n) + std::log(Ci_n));
        
        mu_k_t[k + 1] = x_k_t_1_ + dt * (I_t / Ci_n - x_k_t_1_ / tau_n);
    }
    
    return mu_k_t;
}

double ParticleFiltering::g(const double & I_t, const arma::colvec & x_k_t) 
{
    double VT_ = OCV + I_t * R0;
    for (unsigned int k = 0; k < K; k++) 
    {
        VT_ += x_k_t[1 + k];
    }
    
    return VT_;
}

// Constructor
ParticleFiltering::ParticleFiltering(const arma::colvec & V_, const arma::colvec & I_, 
                                     const arma::colvec & Temp_, const arma::colvec & Time_,
                                     const arma::colvec & OCV_c_pars_, const arma::colvec & OCV_d_pars_, 
                                     const arma::colvec & RT_pars_,
                                     const arma::colvec & Q_pars_, const arma::mat & L_, 
                                     const double & sigma_, const double & SOC_0_,
                                     const unsigned int & K_, const unsigned int & N_,
                                     const bool & trace_, const unsigned int & trace_limit_) 
{
    //
    V = V_;
    I = I_;
    Temp = Temp_;
    Time = Time_;
    
    OCV_c_pars = OCV_c_pars_;
    OCV_d_pars = OCV_d_pars_;
    
    RT_pars = RT_pars_;
    Q_pars = Q_pars_;
    
    L = L_;
    sigma = sigma_;
    
    K = K_;
    N = N_;
    
    trace = trace_;
    trace_limit = trace_limit_;
    
    T = Time.size(); 
    
    //
    const arma::colvec mu_0(K + 1, arma::fill::zeros);
    X_current = arma::mat(K + 1, N, arma::fill::zeros);
    for (unsigned int n = 0; n < N; n++) 
    {
        arma::colvec nn = mvrnorm(1, mu_0, L);
        X_current.col(n) = nn;
        X_current(0, n) = SOC_0_;
    }
    
    X = std::vector<arma::colvec>(T);
    X[0] = mean(X_current, 2);
    
    Ci = arma::colvec(K, arma::fill::zeros);
    for (unsigned int k = 0; k < K; k++) 
    {
        Ci[k] = RT_pars[5 + K + k];
    }
    
    VT_mu = arma::colvec(T, arma::fill::zeros);
    VT_sd = arma::colvec(T, arma::fill::zeros);
}

// 
void ParticleFiltering::Filter() 
{
    latest_non_zero = 0; 
    for (unsigned int t = 1; t < T; t++) 
    {
        if (trace & ((t + 1 == 2) || (t + 1 == T) || (((t + 1) % trace_limit) == 0))) 
        {
            Rcpp::Rcout << "Iteration: " << t + 1 << " / " << T << "\n";
        }
        
        // 
        const double & I_t = I[t];
        const double & Temp_t = Temp[t];
        const double & V_t = V[t];
        
        dt = Time[t] - Time[t - 1];
        
        // Update capacity
        f_Q(I_t, Temp_t);
        
        // Is the current non-zero?
        if (std::abs(I_t) > 1e-6) 
        {
            latest_non_zero = I_t;
        } 
        
        // Particle filtering stage 1 -- sampling states and calcualting weights.
        arma::mat X_new(K + 1, N, arma::fill::zeros);
        arma::colvec W_new(N, arma::fill::zeros);
        for (unsigned int n = 0; n < N; n++) 
        {
            const arma::colvec & x_k_t_1 = X_current.col(n);
            
            double SOC_t = soc_transform(x_k_t_1[0]);
            update_soc_parameters(I_t, SOC_t);
            
            const arma::colvec & f_x_k_t_1 = f(I_t, x_k_t_1);
            const arma::colvec & x_k_new = mvrnorm(1, f_x_k_t_1, L);
            
            SOC_t = soc_transform(x_k_new[0]);
            update_soc_parameters(I_t, SOC_t);
            const double & g_x_k_t = g(I_t, x_k_new);
            const double & w_k_t = dnorm(V_t, g_x_k_t, sigma);
            
            X_new.col(n) = x_k_new;
            W_new[n] = w_k_t;
            if (n > 0) 
            {
                W_new[n] += W_new[n - 1];
            }
        }
        
        W_new = W_new / W_new[N - 1];
        
        // Particle filtering stage 2 -- re-sampling states based on weights.
        arma::colvec VT_current(N, arma::fill::zeros);
        for (unsigned int n = 0; n < N; n++) 
        {
            arma::colvec u = arma::randu(1);
            unsigned int i = 0;
            while (W_new[i + 1] < u[0]) 
            {
                i++;
            }
            
            const arma::colvec & x_i_t = X_new.col(i);
            X_current.col(n) = x_i_t;
            
            // Storing the expected predicted voltage 
            double SOC_t = soc_transform(x_i_t[0]);
            update_soc_parameters(I_t, SOC_t);
            const double g_x_i_t = g(I_t, x_i_t);
            VT_current[n] = g_x_i_t;
        }
        
        // Calculate averages and standard deviations
        arma::colvec X_average(K + 1, arma::fill::zeros);
        double VT_sum = 0.0;
        double VT_sum_squared = 0.0;
        for (unsigned int n = 0; n < N; n++) 
        {
            X_average += X_current.col(n) / N;
            VT_sum += VT_current[n];
            VT_sum_squared += VT_current[n] * VT_current[n];
        }
        
        double VT_variance = (VT_sum_squared - std::exp(2.0 * std::log(VT_sum) - std::log(N))) / (N - 1);
        if (VT_variance < 2e-16) 
        {
            VT_variance = 2e-16;
        }
        
        X[t] = X_average;
        VT_mu[t] = VT_sum / N;
        VT_sd[t] = std::pow(VT_variance, 0.5);
        
        // Rcpp::Rcout << "Xbar = " << X_average.t() << "\n" 
        //             << "VTbar = " << VT_sum / N << "\n"
        //             << "VT_variance = " << VT_variance << "\n";
        // Rcpp::stop("...");
    }
}

// R wrapper functions
//[[Rcpp::export()]]
Rcpp::List particle_filtering_cpp(const arma::colvec & V, const arma::colvec & I, 
                                  const arma::colvec & Temp, const arma::colvec & Time, 
                                  const arma::colvec & OCV_c_pars, const arma::colvec & OCV_d_pars, 
                                  const arma::colvec & RT_pars,
                                  const arma::colvec & Q_pars, const arma::mat & L, const double & sigma,
                                  const double & SOC_0, const unsigned int & K, const unsigned int & N, 
                                  const bool & trace, const unsigned int & trace_limit) 
{
    ParticleFiltering PF(V, I, Temp, Time, OCV_c_pars, OCV_d_pars, 
                         RT_pars, Q_pars, L, sigma, 
                         SOC_0, K, N, trace, trace_limit);
    PF.Filter();
    
    return Rcpp::List::create(Rcpp::Named("Vmean") = PF.VT_mu,
                              Rcpp::Named("Vsd") = PF.VT_sd,
                              Rcpp::Named("Xmean") = PF.X);
}

