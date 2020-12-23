#' @title Filtering of state-of-charge and health.
#' 
#' @description Particle filter for determining the latent state-of-charge, and non-linear optimisation for determining state-of-health. 
#' 
#' @param Current The observed current profile.
#' @param Voltage The observed voltage profile.
#' @param Temperature The observed temperature profile.
#' @param Time The time increments.
#' @param SOC_0 The initial state-of-charge.
#' @param K The number of RC branches used in the state-space model.
#' @param N The number of particles.
#' @param trace TRUE/FALSE: Show trace? If 'numeric' provided, then the trace is shown every 'trace' elements.
#' @param ... Additional arguments including initialisation of parameters.
#' 
#' @return List of estimated parameters, predicted voltage, and estimated SOC.
#' @export
non_linear_particle_filter <- function(I, V, Temperature, Time, SOC_0, K = 2, N = 10000, trace = 100, ...) {
    dots <- list(...)
    if (is.null(dots$OCV)) {
        OCV_0_ <- c(0.3, -0.4, -0.4)
        OCV_0 <- matrix(rep(OCV_0_, times = 2), ncol = 2)
    }
    else {
        OCV_0 <- dots$OCV
    }
    
    if (is.null(dots$R0)) {
        R0_0 <- c(-6, -0.4, -0.4, 10)
    }
    else {
        R0_0 <- dots$R0
    }
    
    if (is.null(dots$RT)) {
        RT_0 <- c(c(-6, -0.4, -0.4, 10), rep(1 / N, N), rep(1000, N))
    }
    else {
        RT_0 <- dots$RT
    }
    
    if (is.null(dots$Q)) {
        Q_0 <- rep(1, 5)
    }
    else {
        Q_0 <- dots$Q
    }
    
    if (is.null(dots$L)) {
        L_0 <- diag(0.1, 1 + K)
    }
    else {
        L_0 <- dots$L_0
    }
    
    if (is.null(dots$sigma)) {
        sigma_0 <- 0.01
    }
    else {
        sigma_0 <- dots$sigma
    }
    
    if (is.numeric(trace)) {
        trace_limit = trace
        trace = TRUE
    } 
    else if (is.logical(trace)) {
        trace_limit = ceiling(length(V) / 10)
    }
    else {
        trace_limit = 1 # not value doesn't matter, but is needed
        trace = FALSE
    }
    
    
    res <- BatteryPrognostics:::particle_filtering_cpp(I, V, Temperature, Time,
                                                       OCV_0[, 2], OCV_0[, 1],
                                                       ## R0_0, 
                                                       RT_0, Q_0,
                                                       L_0, sigma_0, tan(pi * (SOC_0 - 0.5)),
                                                       K, N, trace, trace_limit)
    return(res)
}
