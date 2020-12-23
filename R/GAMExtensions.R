#' @title GAM parameter wrap 
#' 
#' @description A wrapping function creating a list of GAM fits of a 'resistance_parameters'-object.
#' 
#' @param resistance_parameters A list containing a vector of weeks and a matrix of estimated parameters constructed using \link{estimate_parameters_resistance_weeks}.
#' 
#' @return A list of GAM fits.
#' @export
gam_parameter_fit <- function(resistance_parameters) {
    res_tibble <- resistance_parameters %>% 
        as_tibble() %>% 
        spread(key = ParameterName, value = Parameter)
    
    res <- vector("list", dim(res_tibble)[2] - 1)
    res_names <- colnames(resistance_parameters$Parameters)
    for (i in 1:(dim(res_tibble)[2] - 1)) {
        formula <- as.formula(paste("`", res_names[i], "` ~ s(Time)", sep = ""))
        res[[i]] <- suppressWarnings(gam(formula, data = res_tibble))
    }
    
    names(res) <- res_names
    class(res) <- "gam_list"
    return(res)
}

#' @export
predict.gam_list <- function(object, ...) {
    res <- lapply(object, predict, ...)
    return(res)
}

#' @export
ggplot.gam_list <- function(gam_list, ...) {
    arg_list <- list(...)
    if (is.null(arg_list$type)) {
        type = "link"
    }
    else {
        type = arg_list$type
    }
    
    if (is.null(arg_list$se.fit)) {
        se.fit = TRUE
    }
    else {
        se.fit = arg_list$se.fit
    }
    
    pred_list <- predict(gam_list, type = type, se.fit = se.fit)
    
    data_tibble <- vector("list", length(gam_list))    
    parameter_tibble <- vector("list", length(gam_list))
    
    if (!is.null(arg_list$newdata)) {
        newdata = arg_list$newdata
        predict_tibble <- vector("list", length(gam_list))
    }
    else {
        newdata = NULL
    }
    
    for (i in seq_along(gam_list)) {
        Week_i = gam_list[[i]]$data$Time 
        ParameterName_i = names(gam_list[i])
        data_tibble[[i]] <- tibble(Time = Week_i, Parameter = gam_list[[i]]$y, ParameterName = ParameterName_i)
        parameter_tibble[[i]] <- tibble(Time = Week_i, Parameter = pred_list[[i]]$fit, ParameterName = ParameterName_i)
        
        if (se.fit) {
            if (is.null(arg_list$cv)) {
                cv = 1.96
            } 
            else {
                cv = arg_list$cv
            }
            
            parameter_tibble[[i]] <- parameter_tibble[[i]] %>% mutate(LowerCI = Parameter - cv * pred_list[[i]]$se.fit[, 1], 
                                                                      UpperCI = Parameter + cv * pred_list[[i]]$se.fit[, 1]) 
        }
        
        if (!is.null(newdata)) {
            predict_tibble[[i]] <- newdata %>% mutate(Parameter = predict(gam_list[[i]], newdata = newdata), ParameterName = ParameterName_i)
        }
    }
    
    data_tibble <- data_tibble %>% bind_rows() %>% 
        mutate(ParameterName = factor(ParameterName, levels = c("sd", "beta_0", "beta_1", "beta_2")))
    parameter_tibble <- parameter_tibble %>% bind_rows() %>% 
        mutate(ParameterName = factor(ParameterName, levels = c("sd", "beta_0", "beta_1", "beta_2")))
    
    new_labels <- as_labeller(c(beta_0 = "beta[0]", beta_1 = "beta[1]", beta_2 = "beta[2]", sd = "Standard~deviation"), label_parsed)
    p <- ggplot(data_tibble, aes(x = Time, y = Parameter)) + 
        facet_wrap(~ParameterName, ncol = 2, nrow = 2, scales = "free_y", labeller = new_labels) + 
        ylab("") + theme_bw()
    
    if (se.fit) {
        p <- p + 
            geom_ribbon(data = parameter_tibble, aes(ymin = LowerCI, ymax = UpperCI), colour = "transparent", fill = "grey50", alpha = 0.4) + 
            geom_line(data = parameter_tibble, colour = "dodgerblue2", size = 1.2) +
            geom_line(data = parameter_tibble, aes(y = LowerCI), colour = "dodgerblue2", size = 0.7, linetype = "dashed") + 
            geom_line(data = parameter_tibble, aes(y = UpperCI), colour = "dodgerblue2", size = 0.7, linetype = "dashed") 
    }
    else {
        p <- p + geom_line(data = parameter_tibble, colour = "dodgerblue2", size = 1.2)
    }
    
    if (!is.null(newdata)) {
        predict_tibble <- predict_tibble %>% bind_rows() %>% 
            mutate(ParameterName = factor(ParameterName, levels = c("sd", "beta_0", "beta_1", "beta_2")))
        
        p <- p + geom_line(data = predict_tibble, colour = "orangered2", size = 1.2)        
    }
    
    p <- p + geom_point()
    return(p)
}

#' @export
simulate.gam_list <- function(object, nsim = 1, seed = NULL, ...) {
    set.seed(seed)
    W <- object[[1]]$data$Time
    new_data <- data.frame(Time = seq(min(W), max(W)))
    
    res <- vector("list", length(object))
    for (j in seq_along(object)) {
        x <- object[[j]]
        fit <- predict(x, newdata = new_data)
        se <- sqrt(deviance(x) / df.residual(x))
        
        x_sim <- matrix(0, ncol = nsim, nrow = length(fit))
        for (i in 1:nsim) {
            if (names(object[j]) %in% c("beta_1", "beta_2")) {
                x_sim_i <- fit + rnorm(length(fit), mean = 0, sd = se)
                re_roll <- which(x_sim_i > 0)
                while (length(re_roll) > 0) {
                    x_sim_i[re_roll] <- fit[re_roll] + rnorm(length(re_roll), mean = 0, sd = se)
                    re_roll <- which(x_sim_i > 0)
                }
            }
            else {
                x_sim_i <- fit + rnorm(length(fit), mean = 0, sd = se)
            }
            
            x_sim[, i] <- x_sim_i
        }
        
        res[[j]] <- x_sim
    }
    
    names(res) <- names(object)
    res_list <- list(W = new_data$Time, 
                     Simulations = res, 
                     GamList = object)
    class(res_list) <- "simulate_list"
    return(res_list)
}

#' @export 
ggplot.simulate_list <- function(simulate_list) {
    W <- simulate_list$W
    
    data_tibble <- lapply(seq_along(simulate_list$Simulations), function(i) {
        sim <- simulate_list$Simulations[[i]]
        res <- sim %>% as_tibble() %>% 
            mutate(Time = W) %>% 
            gather(key = ParameterName, value = Parameter, -Time) %>% 
            mutate(Simulation = as.numeric(str_sub(ParameterName, 2)), 
                   ParameterName = names(simulate_list$Simulations[i]))
        return(res)
    }) %>% bind_rows() %>% 
        mutate(ParameterName = factor(ParameterName, levels = c("sd", "beta_0", "beta_1", "beta_2")))
    
    pred_gam_list <- predict(simulate_list$GamList, newdata = data.frame(Time = simulate_list$W)) %>% 
        enframe(name = "ParameterName", value = "Parameter") %>% 
        mutate(Time = list(W)) %>% 
        unnest() %>% 
        mutate(ParameterName = factor(ParameterName, levels = c("sd", "beta_0", "beta_1", "beta_2")))
    
    new_labels <- as_labeller(c(beta_0 = "beta[0]", beta_1 = "beta[1]", beta_2 = "beta[2]", sd = "Standard~deviation"), label_parsed)
    p <- ggplot(data_tibble, aes(x = Time, y = Parameter)) + 
        geom_line(aes(group = Simulation), alpha = 0.2) + 
        facet_wrap(~ParameterName, ncol = 2, nrow = 2, scales = "free_y", labeller = new_labels) + 
        geom_line(data = pred_gam_list, alpha = 1, colour = "dodgerblue2", size = 1.2) +
        ylab("") + theme_bw()
    
    return(p)
}


#' @title Interpolate 'resistance_parameters'
#' 
#' @description Interpolates parameters of missing weeks using GAM.
#' 
#' @param resistance_parameters A list containing a vector of weeks and a matrix of estimated parameters constructed using \link{estimate_parameters_resistance_weeks}.
#' @param smooth TRUE/FALSE: Should the fitted 'gam_list' be used to smooth the estimated parameters.
#' @param limits A vector of one or two elements indicating the limits of the 'W'.
#' 
#' @return An interpolated object of class 'resistance_parameters'.
#' @export
resistance_parameters_interpolate <- function(resistance_parameters, smooth = FALSE, 
                                              limits = NULL) {
    gam_list <- gam_parameter_fit(resistance_parameters)
    
    if (length(limits) == 1) {
        Wmin = 1
        Wmax = limits
    }
    else if (length(limits) == 2) {
        Wmin = limits[1]
        Wmax = limits[2]
    }
    else {
        Wmin = 1
        Wmax = resistance_parameters$W[length(resistance_parameters$W)]
    }
    
    W <- seq(Wmin, Wmax, 1)
    pred_all <- predict(gam_list, newdata = tibble(Time = W)) %>% 
        as_tibble() %>% 
        mutate(Time = 1:n())
    
    if (smooth) {
        est_pars_ <- pred_all %>% dplyr::select(-Time) %>% as.matrix()
    }
    else {
        not_in <- pred_all %>% filter(!(Time %in% resistance_parameters$W))
        est_pars_ <- matrix(0, ncol = ncol(resistance_parameters$Parameters), nrow = length(W)) 
        est_pars_[-not_in$Time, ] <- resistance_parameters$Parameters
        est_pars_[not_in$Time, ] <- not_in %>% dplyr::select(-Time) %>% as.matrix()
    }
    
    colnames(est_pars_) <- names(gam_list)
    res <- list(W = W, Parameters = est_pars_)
    
    res$BarlettsTest <- resistance_parameters$BarlettsTest
    
    class(res) <- 'resistance_parameters'
    return(res)
}
