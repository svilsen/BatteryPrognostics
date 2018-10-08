#' @title Filtering of state-of-charge and health.
#' 
#' @description Dual extended Kalman filter for determining the latent state-of-charge and health state-equation. 
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param theta_0 The initial values of the parameter vector.
#' @param ocv_parameters The parameters used for the third order polynomial used to model the OCV.
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
SODEFKFilter <- function(I, V, theta_0, ocv_parameters, SOC_0, C_max, eta,
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
    
    res <- SODEFKFilterCpp(I, V, theta_0, ocv_parameters, SOC_0, C_max, eta, 
                           P, Q, R, dt, K, trace, traceLimit)
    return(res)
}


#' @title Smoothing of state-of-charge and health.
#' 
#' @description Dual extended Kalman filter and smoother for determining the latent state-of-charge and -health. 
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param theta_0 The initial values of the parameter vector.
#' @param ocv_parameters The parameters used for the third order polynomial used to model the OCV.
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
#' @return A list containing the voltage, the smoothed voltage, the SOC, the polarisation voltages, and the simulated parameters.
#' @export
SODEFKSmooth <- function(I, V, theta_0, ocv_parameters, SOC_0, C_max, eta,
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
    
    res <- SODEFKSmoothCpp(I, V, theta_0, ocv_parameters, SOC_0, C_max, eta,
                           P, Q, R, dt, K, trace, traceLimit)
    return(res)
}


#### TEST

#' @title Filtering of state-of-charge and health.
#' 
#' @description Dual extended Kalman filter for determining the latent state-of-charge and health state-equation. 
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param theta_0 The initial values of the parameter vector.
#' @param ocv_parameters The parameters used for the third order polynomial used to model the OCV.
#' @param SOC_0 The initial value of the SOC.
#' @param C_max The maximum value of the capacity.
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
SODEFKTEST <- function(I, V, theta_0, ocv_parameters, SOC_0, C_max,
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
    
    res <- SODEFKFilterCpp(I, V, theta_0, ocv_parameters, SOC_0, C_max, 
                           P, Q, R, dt, K, trace, traceLimit)
    return(res)
}