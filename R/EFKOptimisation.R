initialise_parameters <- function(K) {
    pars <- vector("list", 2 * K + 3)
    lower <- vector("list", 2 * K + 3)
    upper <- vector("list", 2 * K + 3)
    pars_index <- c()
    
    pars[[1]] <- c(0.0, 0.0, 0.0053, 0.0)
    lower[[1]] <- rep(0, length(pars[[1]]))
    upper[[1]] <- rep(Inf, length(pars[[1]]))
    pars_index <- c(pars_index, rep(1, length(pars[[1]])))
    
    for (k in 1:K) {
        pars[[2 * k]] <- c(0.0, 0.0, 0.0036, 0.0082)
        lower[[2 * k]] <- rep(0, length(pars[[2 * k]]))
        upper[[2 * k]] <- rep(Inf, length(pars[[2 * k]]))
        pars_index <- c(pars_index, rep(2*k, length(pars[[2 * k]])))
        
        pars[[2 * k + 1]] <- c(2.62, 1.86, 2.89, 0.61)
        lower[[2 * k + 1]] <- rep(0, length(pars[[2 * k + 1]]))
        upper[[2 * k + 1]] <- rep(Inf, length(pars[[2 * k + 1]]))
        pars_index <- c(pars_index, rep(2*k + 1, length(pars[[2 * k + 1]])))
    }
    
    pars[[2 * K + 2]] <- c(0.0, 0.0, 2.5, 0.0)
    lower[[2 * K + 2]] <- c(rep(0, length(pars[[2 * K + 2]]) - 1), -Inf)
    upper[[2 * K + 2]] <- rep(Inf, length(pars[[2 * K + 2]]))
    pars_index <- c(pars_index, rep(2 * K + 2, length(pars[[2 * K + 2]])))
    
    pars[[2 * K + 3]] <- c(3.11602927, 0.20411943, -0.12431090, 0.02695484, 0.01) 
    lower[[2 * K + 3]] <- c(rep(-Inf, length(pars[[2 * K + 3]]) - 1), 0.0)
    upper[[2 * K + 3]] <- rep(Inf, length(pars[[2 * K + 3]]))
    pars_index <- c(pars_index, rep(2 * K + 3, length(pars[[2 * K + 3]])))
    
    res <- list(unlist(pars), pars_index, unlist(lower), unlist(upper))
    return(res)
} 

#' Optimise EFK
#' 
#' @description Parameter optimisation of the SO-EFK function.
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param SOC_0 The initial value of the SOC.
#' @param dt Step-size.
#' @param K Number of RC branches in the EEC.
#' 
#' @return A list of the optimal parameter values.
#' @export
SOEFK_optimisation <- function(I, V, SOC_0 = 0.0, dt = 1, K = 2) {
    #
    f <- function(par, pars_index, I, IF, V, SOC_0, dt, K) {
        pars_ <- lapply(unname(split(par, pars_index)), function(x) matrix(x, ncol = 1))
        V_hat <- SOEFK(I, IF, V, pars_, SOC_0, dt, K, FALSE, 1)
        
        res <- sum((V - V_hat$V_hat)^2)
        if (any(is.nan(res)) | any(is.infinite(res))) 
            return(1000000)
        
        return(res)
    } 
    
    # Creating IF
    IDifference = diff(I)
    IF = I
    
    last_non_zero_I = I[1]
    for (k in seq(2, length(I))) {
        if (abs(IDifference[k - 1]) > 1e-6) {
            if (abs(I[k - 1]) > 1e-6) {
                last_non_zero_I = I[k - 1]
            }
        }
        
        if (abs(I[k]) < 1e-6) {
            IF[k] = last_non_zero_I
        }
    }

    # Creating initial par and bounds
    initial_parameters_bound <- initialise_parameters(K)
    initial_parameters <- initial_parameters_bound[[1]]
    pars_index <- initial_parameters_bound[[2]]
    lower_bounds <- initial_parameters_bound[[3]]
    upper_bounds <- initial_parameters_bound[[4]]
    
    # 
    opt <- optim(par = initial_parameters, fn = f, #gr = g,
                 pars_index = pars_index, I = I, IF = IF, V = V, SOC_0 = SOC_0, dt = dt, K = K,
                 method = "L-BFGS-B", lower = lower_bounds, upper = upper_bounds, 
                 control = list(trace = 1, REPORT = 1, maxit = 1000))
    
    pars_ <- lapply(unname(split(opt$par, pars_index)), function(x) matrix(x, ncol = 1))
    V_hat <- SOEFK(I, IF, V, pars_, SOC_0, dt, K, FALSE, 1)
    
    res_list <- list(Time = 1:length(V_hat$V_hat), Simulated = V_hat$V_hat, Measured = V, EstimatedParameters = pars_, 
                     LogLikelihood = f(unlist(pars_), pars_index, I, IF, V, SOC_0, dt, K)) 
    return(res_list)
}