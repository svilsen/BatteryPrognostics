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
#' @param trace_limit Used to limit trace.
#' 
#' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
#' @export
DEKF <- function(I, V, theta_0, ocv_0, SOC_0, C_max, eta,
                 P = NULL, Q = NULL, R = NULL,
                 dt = 1, K = 2, trace = TRUE, trace_limit = 1e6) {
    if (is.null(P)) {
        P <- list(diag(10, K + 1), diag(10, length(theta_0)))
    }
    else if ((!is.list(P)) | (length(P) != 2)) {
        stop(paste("P is a list:", is.list(P), ":: P contains two elements:", length(P) != 2))
    }
    
    if (is.null(Q)) {
        Q <- list(diag(1e-6, K + 1), diag(1e-6, length(theta_0)))
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
    
    res <- SOCDEFKFilterCpp(I, V, theta_0, ocv_0, SOC_0, C_max, P, Q, R, dt, K, trace, trace_limit)
    return(res)
}

#' @title Filtering of state-of-charge and health.
#' 
#' @description Dual extended Kalman filter for determining the latent state-of-charge and health state-equation. 
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param Temp The observed temperature.
#' @param Time The observed time.
#' @param theta_0 The initial values of the parameter vector.
#' @param ocv_0 The parameters used for the third order polynomial used to model the OCV.
#' @param SOC_0 The initial value of the SOC.
#' @param dual TRUE/FALSE: Should it be a dual, or a single, EKF?
#' @param P The prior conditional variance of the state equation.
#' @param Q The variance of the state equation.
#' @param R The variance of the observation error. 
#' @param K Number of RC branches in the EEC.
#' @param trace TRUE/FALSE: Show trace?
#' @param trace_limit Used to limit trace.
#' 
#' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
#' @export
ADEKF <- function(I, V, Temp, Time, theta_0, ocv_0, SOC_0, dual = TRUE,
                  P = NULL, Q = NULL, R = NULL,
                  K = 1, trace = TRUE, trace_limit = 1e6) {
    if (is.null(P)) {
        P <- list(diag(10, K + 1), diag(10, length(theta_0)))
    }
    else if ((!is.list(P)) | (length(P) != 2)) {
        stop(paste("P is a list:", is.list(P), ":: P contains two elements:", length(P) != 2))
    }
    
    if (is.null(Q)) {
        Q <- list(diag(1e-6, K + 1), diag(1e-6, length(theta_0)))
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
    
    res <- SOCADEFKFilterCpp(I, V, Temp, Time, theta_0, ocv_0, SOC_0, P, Q, R, K, dual, trace, trace_limit)
    return(res)
}


#' @title Filtering of state-of-charge and health.
#' 
#' @description Dual extended Kalman filter for determining the latent state-of-charge and health state-equation. 
#' 
#' @param I A current profile.
#' @param V The observed voltage.
#' @param Temp The observed temperature.
#' @param Time The observed time.
#' @param ocv_0 The parameters used for the third order polynomial used to model the OCV.
#' @param theta_0 The initial values of the parameter vector.
#' @param SOC_0 The initial value of the SOC.
#' @param K Number of RC branches in the EEC.
#' @param sigma_0 A vector of parameters used for the unscented transformation: (alpha, beta, kappa).
#' @param R_X The variance of the state-equation. 
#' @param R_V The variance of the observation-equation. 
#' @param trace TRUE/FALSE: Show trace?
#' @param trace_limit Used to limit trace.
#' 
#' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
#' @export
ADUKF <- function(I, V, Temp, Time, ocv_0, theta_0,
                  SOC_0 = 0.0, K = 1, sigma_0 = NULL,
                  R_X = NULL, R_V = NULL,
                  trace = TRUE, trace_limit = 1e6) {
   if (is.null(R_X)) {
        R_X <- diag(c(0.01, rep(1, K)))
    }
    else if ((!is.matrix(R_X)) & (dim(R_X)[1] != dim(R_X)[2]) & (dim(R_X)[1] != (K + 1))) {
        stop("'R_X' needs to be a '(K + 1) x (K + 1)' matrix.")
    }
    
    if (is.null(R_V)) {
        R_V <- diag(1, 1)
    }
    else if ((!is.matrix(R_V)) & (dim(R_V)[1] != dim(R_V)[2]) & (dim(R_V)[1] != 1)) {
        stop("'R_V' needs to be a '1 x 1' matrix.")
    }
    
    if (is.null(sigma_0)) {
        sigma_0 <- c(1, 2, 0)
    }
    else if (length(sigma_0) > 3)  {
        sigma_0 <- sigma_0[1:3]
    }
    else if (length(sigma_0) < 3) {
        sigma_ <- sigma_0
        sigma_0 <- c(1, 2, 0)
        sigma_0[seq_along(sigma_)] <- sigma_
    }
    
    res <- ADUKF_Cpp(I, V, Temp, Time, ocv_0, theta_0, sigma_0, R_X, R_V, SOC_0, K, trace, trace_limit)
    return(res)
}

