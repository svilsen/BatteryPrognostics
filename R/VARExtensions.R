VAR_resistance_parameters_control <- function(p = 1, type = c("const", "trend", "both", "none"),
                                              season = NULL, exogen = NULL, lag.max = NULL,
                                              ic = c("AIC", "HQ", "SC", "FPE")) {
    list(p = p, type = type, season = season, exogen = exogen, lag.max = lag.max, ic = ic)
}

#' @title Wrapper for VAR
#' 
#' @description A wrapper for the VAR function for the class 'resistance_parameters'. The function also interpolates weeks not found in the provided 'resistance_parameters'-object.
#' 
#' @param resistance_parameters A list containing a vector of weeks and a matrix of estimated parameters constructed using \link{estimate_parameters_resistance_weeks}.
#' @param ... Arguments passed to the \link{VAR}-function.
#' 
#' @return A \link{varest}-object.
#' @export
VAR_resistance_parameters <- function(resistance_parameters, ...) {
    dots <- list(...)
    
    if (is.null(dots$equal_variance)) {
        equal_variance <- FALSE
    }
    else {
        equal_variance <- dots$equal_variance
        dots <- dots[-which(names(dots) == "equal_variance")]
    }
    
    input <- do.call("VAR_resistance_parameters_control", dots)
    
    parameter_names <- colnames(resistance_parameters$Parameters)
    if ((!as.logical(resistance_parameters$BarlettsTest["Rejected"])) || equal_variance) {
        parameter_names <- parameter_names[-which(parameter_names == "sd")]
    }
    
    input$y <- resistance_parameters$Parameters[, parameter_names]
    VAR_fit <- do.call("VAR", input)
    return(VAR_fit)
}

#' @title Predict posterior of VAR
#' 
#' @description Simulates observations from the posterior predictive distribution of a VAR. In particular, for an object of the class '\link{varest}' from the \link{vars} package.
#' 
#' @param object A '\link{varest}' object.
#' @param n_ahead The number of forward predictions.
#' @param n_sims The number of simulations. 
#' 
#' @export
predict_posterior_VAR <- function(object, n_ahead = 1, n_sims = 10000) {
    K <- object$K
    p <- object$p
    
    B <- Bcoef(object)
    
    col_names <- colnames(object$datamat[, 1:K])
    Z <- as.matrix(object$datamat[, -c(1:K)])
    ZZ <- t(Z) %*% Z
    iZZ <- solve(ZZ)
    
    T_end <- Z[dim(Z)[1], "trend"]
    residual <- resid(object) 
    
    S_ <- t(residual) %*% residual
    
    res <- vector("list", n_sims)
    for (r in 1:n_sims) {
        Psi <- BatteryPrognostics:::riwish_(v = T_end - p, S = S_)
        Gamma <- BatteryPrognostics:::rmatrixnorm_(B, Psi, iZZ)
        
        Z_new <- c(Z[dim(Z)[1], ] %*% t(B), Z[dim(Z)[1], ])[1:(K * p)]
        u <- mvtnorm::rmvnorm(n_ahead, rep(0, K), Psi)
        beta_ahead <- matrix(NA, ncol = K, nrow = n_ahead)
        for (h in 1:n_ahead) {
            Z_new <- c(Z_new, 1, T_end + h)
            beta_ahead[h, ] <- Z_new %*% t(B) + u[h, ]
            
            Z_new <- c(beta_ahead[h, ], Z_new)
            Z_new <- Z_new[1:(K * p)]
        }
        
        colnames(beta_ahead) <- col_names
        res[[r]] <- beta_ahead
    }
    
    res_list <- list(Simulations = res, Model = object)
    class(res_list) <- "post_pred_varest"
    return(res_list)
}

#' @export
as_resistance_parameters <- function(object, ...) {
    UseMethod("as_resistance_parameters")
}

#' @export
as_resistance_parameters.default <- function(object, ...) {
    cat(paste("Function not implemented for class '", class(object), "'. \n", sep = ""))
}

#' @export
as_resistance_parameters.post_pred_varest <- function(object, ...) {
    dots <- list(...)
    
    initial <- 1
    if(!is.null(dots$initial)) {
        initial <- dots$initial
    }
    
    ##    
    original_data <- object$Model$y
    original_W <- 1:dim(original_data)[1] + initial - 1
    
    ##
    col_names <- colnames(object$Simulations[[1]])
    pred_data <- do.call("cbind", lapply(col_names, function(k) {
        apply(do.call("cbind", lapply(object$Simulations, function(xx) xx[, k])), 1, mean)
    })) 
    
    pred_W <- 1:dim(pred_data)[1] + original_W[length(original_W)]
    
    res <- list("W" = c(original_W, pred_W), "Parameters" = rbind(original_data, pred_data), 
                "BarlettsTest" = NA)
    class(res) <- "resistance_parameters"
    return(res)
}


#' @export
ggplot.post_pred_varest <- function(post_pred_var, ...) {
    dots <- list(...)
    
    K <- post_pred_var$Model$K
    W <- post_pred_var$Model$datamat[, dim(post_pred_var$Model$datamat)[2]]
    
    obs_tibble <- post_pred_var$Model$y %>% 
        as_tibble() %>% 
        mutate(Week = 1:n()) %>% 
        gather(ParameterName, Parameter, -Week)
    
    fitted_tibble <- fitted(post_pred_var$Model) %>% 
        as_tibble() %>% 
        mutate(Week = W) %>% 
        gather(ParameterName, Parameter, -Week)
    
    n_ahead <- dim(post_pred_var$Simulations[[1]])[1]
    pred_weeks <- (W[length(W)] + 1):(W[length(W)] + n_ahead)
    pred_tibble <- lapply(post_pred_var$Simulations, function(xx) xx %>% as_tibble() %>% mutate(Week = pred_weeks)) %>% 
        enframe(name = "Simulation", value = "Parameter") %>% 
        unnest() %>% 
        gather(ParameterName, Parameter, -Week, -Simulation) %>% 
        group_by(Week, ParameterName) %>% 
        summarise(LowerBound = quantile(Parameter, probs = c(0.025)), 
                  UpperBound = quantile(Parameter, probs = c(0.975)), 
                  Parameter = mean(Parameter))
    
    if (!is.null(dots$weeks) & ((length(dots$week) == 2) || (length(dots$week) == 1))) {
        weeks <- dots$weeks
        if (length(weeks) == 1) {
            weeks <- c(1, weeks)
        }
        
        obs_tibble <- obs_tibble %>% filter(Week >= weeks[1], Week <= weeks[2])
        fitted_tibble <- fitted_tibble %>% filter(Week >= weeks[1], Week <= weeks[2])
        pred_tibble <- pred_tibble %>% filter(Week >= weeks[1], Week <= weeks[2])
    }
    
    new_labels <- as_labeller(c(beta_0 = "beta[0]", beta_1 = "beta[1]", beta_2 = "beta[2]", sd = "Standard~deviation"), 
                              label_parsed)
    p <- ggplot(data = obs_tibble, aes(x = Week, y = Parameter)) + 
        geom_line(size = 0.7) + 
        geom_line(data = fitted_tibble, colour = "dodgerblue2", size = 0.8) + 
        geom_line(data = pred_tibble, colour = "red", size = 1) + 
        geom_line(data = pred_tibble, aes(y = LowerBound), colour = "red", size = 0.7, linetype = "dashed") + 
        geom_line(data = pred_tibble, aes(y = UpperBound), colour = "red", size = 0.7, linetype = "dashed") + 
        facet_wrap(~ ParameterName, scales = "free_y", labeller = new_labels) + 
        ylab("") + 
        theme_bw()
    
    return(p)
}

#' @export 
summary.post_pred_varest <- function(post_pred_var) {
    post_pars <- colnames(post_pred_var$Model$y)
    W <- post_pred_var$Model$datamat[, dim(post_pred_var$Model$datamat)[2]]
    n_ahead <- dim(post_pred_var$Simulations[[1]])[1]
    pred_weeks <- (W[length(W)] + 1):(W[length(W)] + n_ahead)
    
    pred_tibble <- lapply(post_pred_var$Simulations, function(xx) xx %>% as_tibble() %>% mutate(Week = pred_weeks)) %>% 
        enframe(name = "Simulation", value = "Parameter") %>% 
        unnest() %>% 
        gather(ParameterName, Parameter, -Week, -Simulation) %>% 
        group_by(Week, ParameterName) %>% 
        summarise(`mean` = mean(Parameter), 
                  `MAP` = BatteryPrognostics:::mode_density(Parameter),
                  `2.5%` = quantile(Parameter, probs = c(0.025)),
                  `50%` = median(Parameter), 
                  `97.5%` = quantile(Parameter, probs = c(0.975)))
    
    res_list <- structure(lapply(post_pars, function(ppv) {
        res <- pred_tibble %>% 
            ungroup() %>% 
            filter(ParameterName == ppv) %>% 
            dplyr::select(-ParameterName)
        
        if (ppv != "sd") {
            res <- res %>% mutate_at(vars(-one_of("Week")), exp)
        }
        
        return(res)
    }), .Names = post_pars)
    
    return(res_list)
}
