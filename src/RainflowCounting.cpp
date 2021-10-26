#include <RcppArmadillo.h>

//[[Rcpp::export()]]
arma::mat rainflow_cpp(const arma::colvec & x, const double & flm, const double & ul, const double & pls) {
    const double flmargin = ul - std::fabs(flm);
    const double x_rows = x.n_rows;
    
    arma::mat res = arma::zeros<arma::mat>(x_rows, 5);
    int pr = 0;
    int po = 0;
    int j = -1;
    
    arma::colvec a = arma::zeros<arma::colvec>(x_rows);
    for (int i = 0; i < x_rows; i++) {
        j += 1;
        a[j] = x[pr];
        
        pr += 1;
        
        while ((j >= 2) & (std::fabs(a[j - 1] - a[j - 2]) <= std::fabs(a[j] - a[j - 1]))) {
            double lrange = std::fabs(a[j - 1] - a[j - 2]);
            
            if (j == 2) {
                double avg = 0.5 * (a[0] + a[1]);
                
                double adj_range = lrange * flmargin / (ul - std::fabs(avg));
                double adj_zero_avg_range = lrange * ul / (ul - std::fabs(avg));
                
                a[0] = a[1];
                a[1] = a[2];
                
                j = 1;
                
                if (lrange > 0) {
                    res(po, 0) = lrange;
                    res(po, 1) = avg;
                    res(po, 2) = adj_range;
                    res(po, 3) = pls;
                    res(po, 4) = adj_zero_avg_range;
                    
                    po += 1;
                }
            }
            else {
                double avg = 0.5 * (a[j - 1] + a[j - 2]);
                
                double adj_range = lrange * flmargin / (ul - std::fabs(avg));
                double adj_zero_avg_range = lrange * ul / (ul - std::fabs(avg));
                
                a[j - 2] = a[j];
                
                j -= 2;
                
                if (lrange > 0) {
                    res(po, 0) = lrange;
                    res(po, 1) = avg;
                    res(po, 2) = adj_range;
                    res(po, 3) = 1.0;
                    res(po, 4) = adj_zero_avg_range;
                    
                    po += 1;
                }
            }
        }
    }
    
    for (int i = 0; i < j; i++) {
        double lrange = std::fabs(a[i] - a[i + 1]);
        
        double avg = 0.5 * (a[i] + a[i + 1]);
        double adj_range = lrange * flmargin / (ul - std::fabs(avg));
        double adj_zero_avg_range = lrange * ul / (ul - std::fabs(avg));
        
        if (lrange > 0) {
            res(po, 0) = lrange;
            res(po, 1) = avg;
            res(po, 2) = adj_range;
            res(po, 3) = pls;
            res(po, 4) = adj_zero_avg_range;
            
            po += 1;
        }
    }
    
    return res;
}