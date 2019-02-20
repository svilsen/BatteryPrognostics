#' @title Filtering of state-of-charge and health.
#' 
#' @description Dual extended Kalman filter for determining the latent state-of-charge and health state-equation. 
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param theta_0 The initial values of the parameter vector.
#' @param ocv_0 The parameters used for the third order polynomial used to model the OCV.
#' @param SOC_0 The initial value of the SOC.
#' @param C_max The maximum value of the capacity.
#' @param eta The capacity efficiency.
#' @param P The prior conditional variance of the state equation.
#' @param Q The variance of the state equation.
#' @param R The variance of the observation error. 
#' @param dt Step-size.
#' @param K Number of RC branches in the EEC.
#' @param trace TRUE/FALSE: Show trace?
#' @param traceLimit Used to limit trace.
#' 
#' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
#' @export
SOC_DEFK_filter <- function(I, V, theta_0, ocv_0, SOC_0, C_max, eta,
                         P = NULL, Q = NULL, R = NULL,
                         dt = 1, K = 2, trace = TRUE, traceLimit = 1e6) {
    if (is.null(P)) {
        P <- list(diag(10, K + 1), diag(10, length(theta_0)))
    }
    else if ((!is.list(P)) | (length(P) != 2)) {
        stop(paste("P is a list:", is.list(P), ":: P contains two elements:", length(P) != 2))
    }
    
    if (is.null(Q)) {
        Q <- list(diag(1e-8, K + 1), diag(1e-8, length(theta_0)))
    }
    else if ((!is.list(Q)) | (length(Q) != 2)) {
        stop(paste("Q is a list:", is.list(Q), ":: Q contains two elements:", length(Q) != 2))
    }
    
    if (is.null(R)) {
        R <- list(diag(10, 1), diag(10, 1))
    }
    else if ((!is.list(R)) | (length(R) != 2)) {
        stop(paste("R is a list:", is.list(R), ":: R contains two elements:", length(R) != 2))
    }
    
    res <- SODEFKFilterCpp(I, V, theta_0, ocv_0, SOC_0, C_max, P, Q, R, dt, K, trace, traceLimit)
    return(res)
}

#' @title Filtering of parameters.
#' 
#' @description Unscented Kalman filter for determining the latent parameters. 
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param theta_0 The initial values of the parameter vector.
#' @param P_0 The prior conditional variance of the state equation.
#' @param SOC_0 The initial value of the SOC.
#' @param C_max The maximum value of the capacity.
#' @param eta The capacity efficiency.
#' @param sigma_0 A vector containing the three sigma point parameters.
#' @param alpha An OCV forgetting factor.
#' @param lag The time 'lag' used for parameter smoothing. 
#' @param dt Step-size.
#' @param K Number of RC branches in the EEC.
#' @param trace TRUE/FALSE: Show trace?
#' @param traceLimit Used to limit trace.
#' 
#' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
#' @export
UKF_parameters <- function(I, V, theta_0, P_0, SOC_0, C_max, eta = 1.0,
                         sigma_0 = c(0.1, 2, 10), lag = 4,
                         dt = 1, K = 2, R_0 = NULL, Q_0 = NULL, 
                         trace = TRUE, traceLimit = 1e6) {
    if (is.null(P_0)) {
        P_0 <- diag(c(0.1, rep(10, K), rep(1, K), 1), length(theta_0))
    }
    
    if (is.null(R_0)) {
        R_0 <- diag(10, lag)
    }
     
    if (is.null(Q_0)) {
        Q_0 <- diag(10, length(theta_0))
    }
    
    res <- ParameterUKFCpp(I, V, theta_0, P_0, SOC_0, C_max, eta, 
                           sigma_0, lag, dt, K, R_0, 
                           Q_0, trace, traceLimit)
    return(res)
}


#' @title Smoothing of parameters.
#' 
#' @description Unscented Kalman filter and smoother for determining the latent parameters. 
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param theta_0 The initial values of the parameter vector.
#' @param P_0 The prior conditional variance of the state equation.
#' @param SOC_0 The initial value of the SOC.
#' @param C_max The maximum value of the capacity.
#' @param eta The capacity efficiency.
#' @param sigma_0 A vector containing the three sigma point parameters.
#' @param lag The time 'lag' used for parameter smoothing. 
#' @param M Half the window size used for weighted moving average of parameters.
#' @param dt Step-size.
#' @param K Number of RC branches in the EEC.
#' @param trace TRUE/FALSE: Show trace?
#' @param traceLimit Used to limit trace.
#' 
#' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
#' @export
UKS_parameters <- function(I, V, theta_0, P_0, SOC_0, C_max, eta = 1.0,
                           sigma_0 = c(0.1, 2, 10), lag = 5, M = 10,
                           dt = 1, K = 2, R_0 = NULL, Q_0 = NULL,
                           trace = TRUE, traceLimit = 1e6) {
    if (is.null(P_0)) {
        P_0 <- diag(c(0.1, rep(10, K), rep(1, K), 1), length(theta_0))
    }
    
    if (is.null(R_0)) {
        R_0 <- diag(1e-03, lag)
    }
    
    if (is.null(Q_0)) {
        Q_0 <- diag(1e-02, length(theta_0))
    }
    
    res <- ParameterUKSCpp(I, V, theta_0, P_0, SOC_0, C_max, eta, 
                           sigma_0, lag, M, dt, K, R_0, 
                           Q_0, trace, traceLimit)
    return(res)
}


#' @title Filtering of state-of-charge and health.
#' 
#' @description Dual extended Kalman filter for determining the latent state-of-charge and health state-equation. 
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param theta_0 The initial values of the parameter vector.
#' @param SOC_0 The initial value of the SOC.
#' @param Q_max The maximum value of the capacity.
#' @param eta The capacity efficiency.
#' @param P The prior conditional variance of the state equation.
#' @param dt Step-size.
#' @param K Number of RC branches in the EEC.
#' @param W Size of the windows used to estimate the R and Q.
#' @param trace TRUE/FALSE: Show trace?
#' @param traceLimit Used to limit trace.
#' 
#' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
#' @export
EKF_parameters <- function(I, V, theta_0, SOC_0, Q_max, eta,
                           P = NULL, dt = 1, K = 2, W = 100,
                           trace = TRUE, traceLimit = 1e6) {
    if (is.null(P)) {
        P <- diag(c(0.01, 0.01, 6e-5, 0.001, 0.001), length(theta_0))
    }
    
    Q <- diag(0, length(theta_0))
    R <- diag(0, 1)
    
    res <- ParameterEKFCpp(-I, V, theta_0, SOC_0, Q_max, 
                           P, Q, R, dt, K, W, trace, traceLimit)
    return(res)
}
