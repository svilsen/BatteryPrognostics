#' @title The expected log-resistance.
#' 
#' @param beta_w The mean parameters of the resistance.
#' @param SOC A vector of the SOC.
#' 
#' @return Numeric
#' @export 
resistance_expected <- function(beta_w, SOC) {
    beta_0 <- beta_w[1]
    if (length(beta_w) == 1) {
        beta_1 <- beta_2 <- 1
    }
    else if (length(beta_w) == 2) {
        beta_1 <- beta_2 <- beta_w[2]
    }
    else {
        beta_1 <- beta_w[2]
        beta_2 <- beta_w[3]
    }
    
    mu = beta_0 + beta_1 * log(SOC) + beta_2 * log(1 - SOC)
    return(mu)
}

resistance_expected_vectorised <- function(beta_0, beta_1, beta_2, SOC) {
    if (is.null(beta_2)) {
        beta_2 <- beta_1
    }
    
    mu = beta_0 + beta_1 * log(SOC) + beta_1 * log(1 - SOC)
    return(mu)
}

#' @title Standard deviation of log-resistance
#' 
#' @param beta_w The mean parameters of the resistance.
#' @param R The extracted battery resistance.
#' @param SOC The corresponding SOC of the extracted battery resistance.
#' 
#' @return Numeric
#' @export 
resistance_sd <- function(beta_w, R, SOC) {
    mu <- resistance_expected(beta_w, SOC)
    sd_r <- sqrt(sum((log(R) - mu)^2) / (length(R) - 1))
    return(sd_r)
}

log_likelihood_resistance <- function(pars, R, SOC) {
    log_pars <- log(pars)
    
    mu <- resistance_expected(log_pars, SOC)
    sd_r <- resistance_sd(log_pars, R, SOC)
    return(-sum(dnorm(log(R), mean = mu, sd = sd_r, log = TRUE)))
}

#' @title Estimate resistance parameters
#' 
#' @description Estimates the parameters of the resistance model.
#' 
#' @param R The extracted battery resistance.
#' @param SOC The corresponding SOC of the extracted battery resistance.
#' @param symmetric TRUE/FALSE: Is the model symmetric? That is, can we assume 'b == c'?
#'  
#' @details The resistance model assumes that 'log(R) = a + b log(SOC) + c log(1 - SOC) + epsilon', 
#' where epsilon follows a normal distribution with mean zero and variance sigma^2.
#' 
#' @return The mean parameters and standard deviations of the resistance model.
#' @export 
estimate_parameters_resistance <- function(R, SOC, symmetric = FALSE) {
    symmetric_0_1 <- ifelse(symmetric, 0, 1)
    initial_pars <- rep(0.5, 2 + symmetric_0_1)
    optimal_pars <- solnp(pars = initial_pars,
                          fun = log_likelihood_resistance, R = R, SOC = SOC,
                          LB = rep(0, 2 + symmetric_0_1), UB = rep(1, 2 + symmetric_0_1),
                          control = list(trace = FALSE))
    
    log_pars <- log(optimal_pars$par)
    sd_r = resistance_sd(log_pars, R = R, SOC = SOC)
    
    res <- c(log_pars, sd_r)
    names(res) <- c(paste("beta", 0:(1 + symmetric_0_1), sep = "_"), "sd")
    return(res)
} 

#' @title Symmetry's test
#' 
#' @description Performs likelihood ratio test for symmertry.
#' 
#' @param beta_w A matrix of the estimated parameters in the non-symmetric model. 
#' 
#' @return A p-value for each week in the provided matrix.
#' @export 
symmetry_test <- function(beta_w, R, SOC, W) {
    weeks <- sort(unique(W))
    
    ## Symmetric
    cl <- makeCluster(number_of_clusters)
    registerDoParallel(cl)
    
    symmetric_pars <- foreach (i = seq_along(weeks), .combine = rbind, 
                     .packages = c('BatteryPrognostics', 'Rsolnp')) %dopar%
        estimate_parameters_resistance(R[which(W == weeks[i])], SOC[which(W == weeks[i])], TRUE)
    
    stopCluster(cl)
    rownames(symmetric_pars) <- NULL
    
    ## Non-symmetric
    cl <- makeCluster(number_of_clusters)
    registerDoParallel(cl)
    
    non_symmetric_pars <- foreach (i = seq_along(weeks), .combine = rbind, 
                     .packages = c('BatteryPrognostics', 'Rsolnp')) %dopar%
        estimate_parameters_resistance(R[which(W == weeks[i])], SOC[which(W == weeks[i])], FALSE)
    
    stopCluster(cl)
    rownames(non_symmetric_pars) <- NULL
    
    ## LR-test
    res <- rep(0, length(weeks))
    for (w in seq_along(weeks)) {
        which_w <- which(W == weeks[w])
        log_like_non_symmetric <- log_likelihood_resistance(exp(non_symmetric_pars[w, c("beta_0", "beta_1", "beta_2")]), R[which_w], SOC[which_w])
        log_like_symmetric <- log_likelihood_resistance(exp(symmetric_pars[w, c("beta_0", "beta_1")]), R[which_w], SOC[which_w])
        
        LR_stat <- -2 * (log_like_non_symmetric - log_like_symmetric)
        df <- 1
        res[w] <- 1 - pchisq(LR_stat, df)
    }
    
    names(res) <- weeks
    return(res)
}

#' @title Barlett's test
#' 
#' @description Performs Barlett's test for equal variance in 'k' groups.
#' 
#' @param n A vector containing the counts of the 'k' groups.
#' @param sd_est The estimated standard deviation of the 'k' groups.
#' 
#' @return A list of the elements need to calculate the Barlett's test and the associated p-value.
#' @export 
barletts_test <- function(n, sd_est, alpha_level) {
    N = sum(n)
    k = length(n)
    Sp2 <- sum((n - 1) * sd_est^2) / (N - k)
    logS2 <- sum((n - 1) * log(sd_est^2))
    
    Barletts <- ((N - k) * log(Sp2) - logS2) / (1 + (sum(1 / (n - 1)) - 1 / (N - k)) / (3 * (k - 1)))
    p_value <- pchisq(Barletts, k, lower.tail = F) + (1 - pchisq(Barletts, k, lower.tail = T))

    res <- c("Barletts" = Barletts, "p_value" = p_value, "Rejected" = p_value < alpha_level)
    return(res)
}

#' @title Estimate resistance parameters for all weeks
#' 
#' @description Estimates the parameters of the resistance model for every week on record.
#' 
#' @param R The extracted battery resistance.
#' @param SOC The corresponding SOC of the extracted battery resistance.
#' @param W The corresponding week of the extracted battery resistance.
#' @param number_of_clusters The number of clusters passed to 'doParallel' and 'foreach'.
#' @param symmetric TRUE/FALSE: Is the model symmetric? That is, can we assume 'b == c'?
#' @param alpha_level The significance level used for Barlett's and LR tests.
#' 
#' @details The resistance model assumes that 'log(R) = a + b log(SOC) + c log(1 - SOC) + epsilon', 
#' where epsilon follows a normal distribution with mean zero and variance sigma^2.
#' 
#' @return A matrix with a row for every week containing the mean parameters and standard deviations of the resistance model.
#' @export 
estimate_parameters_resistance_weeks <- function(R, SOC, W, number_of_clusters = 4, symmetric = FALSE, alpha_level = 0.05) {
    # Estimate parameters
    cl <- makeCluster(number_of_clusters)
    registerDoParallel(cl)
    
    weeks <- sort(unique(W))
    pars <- foreach (i = seq_along(weeks), .combine = rbind, 
                    .packages = c('BatteryPrognostics', 'Rsolnp')) %dopar%
        estimate_parameters_resistance(R[which(W == weeks[i])], SOC[which(W == weeks[i])], symmetric)
    
    stopCluster(cl)
    rownames(pars) <- NULL
    
    # Barlett's test
    table_W <- table(W)
    names(dimnames(table_W)) <- "n"
    barletts_test_ <- barletts_test(table_W, pars[, "sd"], alpha_level)
    
    res <- list(W = weeks, Parameters = pars, 
                BarlettsTest = barletts_test_)
    
    class(res) <- "resistance_parameters"
    return(res)
} 

#' @title Estimate resistance parameters for all weeks
#' 
#' @description Estimates the parameters of the resistance model for every week on record. The estimated parameters are stored in a tibble for ease of use.
#' 
#' @param R The extracted battery resistance.
#' @param SOC The corresponding SOC of the extracted battery resistance.
#' @param W The corresponding week of the extracted battery resistance.
#' @param estimated_parameters A matrix of estimated parameters containing one row for every week.
#' @param number_of_clusters The number of clusters passed to 'doParallel' and 'foreach'.
#' @param symmetric TRUE/FALSE: Is the model symmetric? That is, can we assume 'b == c'?
#' 
#' @details The resistance model assumes that 'log(R) = a + b log(SOC) + c log(1 - SOC) + epsilon', 
#' where epsilon follows a normal distribution with mean zero and variance sigma^2.
#' 
#' @return A tibble containing the mean parameters, standard deviations of the resistance model, expected log-resistance, and the minimum possible expected log-resistance.
#' @export 
estimate_parameters_resistance_weeks_tibble <- function(R, SOC, W, estimated_parameters = NULL, number_of_clusters = 4, symmetric = FALSE) {
    if (is.null(estimated_parameters)) 
        estimated_parameters <- estimate_parameters_resistance_weeks(R, SOC, W, number_of_clusters, symmetric)
    
    parameter_tibble <- tibble(Week = estimated_parameters$W, beta_0 = estimated_parameters$Parameters[, 1], beta_1 = estimated_parameters$Parameters[, 2], 
           beta_2 = estimated_parameters$Parameters[, 3], sd = estimated_parameters$Parameters[, 4])
    
    res <- tibble(Week = W, R = R, SOC = SOC) %>% 
        left_join(parameter_tibble, by = "Week") %>% 
        mutate(LogR = log(R), LogRhat = resistance_expected_vectorised(beta_0, beta_1, beta_2, SOC), 
               LogRhatMinimum = beta_0 + (beta_1 + beta_2) * log(0.5))
    
    return(res)
}

#' @title Resistance parameters list
#' 
#' @description A list wrap function for \link{resistance_parameters} class objects.
#'
#' @param ... A series of \link{resistance_parameters} objects.
#' 
#' @return A list of class 'resistance_parameters_list'.
#' @export
resistance_parameters_list <- function(...) {
    res_pars <- list(...)
    class(res_pars) <- "resistance_parameters_list"
    return(res_pars)
}

#' @export
ggplot.resistance_parameters <- function(resistance_parameters) {
    data_tibble <- resistance_parameters %>% as_tibble() %>% 
        mutate(ParameterName = factor(ParameterName, levels = c("sd", "beta_0", "beta_1", "beta_2")))
    
    new_labels <- as_labeller(c(beta_0 = "beta[0]", beta_1 = "beta[1]", beta_2 = "beta[2]", sd = "Standard~deviation"), label_parsed)
    p <- ggplot(data_tibble, aes(x = Week, y = Parameter)) + 
        geom_point() + 
        geom_line() + 
        facet_wrap(~ParameterName, ncol = 2, nrow = 2, scales = "free_y", labeller = new_labels) + 
        ylab("") + 
        theme_bw()
    
    return(p)    
}

#' @export
ggplot.resistance_parameters_list <- function(resistance_parameters_list, ...) {
    dots <- list(...)
    list_names <- dots$list_names
    if (is.null(list_names)) {
        list_names <- paste("Pars:", seq_along(resistance_parameters_list), sep = " ")
    }
    else if (length(list_names) == 1) {
        list_names <- paste(list_names, seq_along(resistance_parameters_list), sep = " ")
    }
    else if (length(list_names) != length(resistance_parameters_list)) {
        stop("The 'list_names' argument should be 'NULL', have length 1, or have length equal to the 'resistance_parameters_list'.")
    }
    
    data_tibble <- lapply(seq_along(resistance_parameters_list), function(xx) {
        resistance_parameters_list[[xx]] %>% as_tibble() %>%
            mutate(ParameterName = factor(ParameterName, levels = c("sd", "beta_0", "beta_1", "beta_2")), 
                   ListName = paste(xx))
    }) %>% bind_rows()
    
    new_labels <- as_labeller(c(beta_0 = "beta[0]", beta_1 = "beta[1]", beta_2 = "beta[2]", sd = "Standard~deviation"), label_parsed)
    p <- ggplot(data_tibble, aes(x = Week, y = Parameter, colour = ListName)) + 
        geom_point() + 
        geom_line() + 
        facet_wrap(~ParameterName, ncol = 2, nrow = 2, scales = "free_y", labeller = new_labels) + 
        scale_colour_discrete(name = "", labels = list_names) + 
        ylab("") + 
        theme_bw()
    
    return(p)    
}

## 
unnormalised_joint_probability <- function(R, SOC, beta_w, sigma_w) {
    if (length(beta_w) == 1) {
        beta_w <- c(beta_w, 1, 1)
    }
    else if (length(beta_w) == 2) {
        beta_w <- c(beta_w, beta_w[2])
    }
    
    like <- dnorm(R, mean = beta_w[1] + beta_w[2] * log(SOC) + beta_w[3] * log(1 - SOC), sd = sigma_w)
    return(like)
}

normalised_probability_w <- function(R, SOC, w, parameter_matrix) {
    n_cols <- dim(parameter_matrix)[2]
    
    beta_w = parameter_matrix[w, 1:(n_cols - 1)]
    sigma_w = parameter_matrix[w, n_cols]
    post_prob <- unnormalised_joint_probability(R, SOC, beta_w, sigma_w) / sum(apply(parameter_matrix, 1, function(xx) unnormalised_joint_probability(R, SOC, xx[1:(n_cols - 1)], xx[n_cols])))
    return(post_prob)
}

resistance_likelihood <- function(log_R, beta_w, sd_w) {
    f <- function(x) {
        mu_w = resistance_expected(beta_w, SOC = x)
        dnorm(log_R, mean = mu_w, sd = sd_w) * dunif(x, 0, 1)
    }
    
    N = 1001
    x <- seq(0 + 1/N, 1 - 1/N, length.out = N)
    y <- f(x)
    
    h  <- x[2] - x[1]
    ca <- (y[2] - y[1]) / h
    cb <- (y[N] - y[N - 1]) / h
    res <- pracma::trapz(x, y) - h^2 / 12 * (cb - ca)
    return(res)
}

resistance_of_given_week <- function(beta_w, sd_w, n_sims, burn_in, trace, traceLimit, w) {
    current_R = rnorm(1, mean = -4, sd = 0.1)
    current_likelihood <- resistance_likelihood(current_R, beta_w, sd_w)
    sampled_R <- rep(NA, n_sims)
    sampled_R[1] <- current_R
    for (i in 2:n_sims) {
        if (trace & ((i %% traceLimit) == 0))
            cat("Week:", w, ":: Iteration:", i, "/", sprintf("%i", n_sims), "\n")
        
        new_R <- rnorm(1, mean = current_R, sd = 0.01)
        new_likelihood <- resistance_likelihood(new_R, beta_w, sd_w)
        
        ratio <- min(ifelse(current_likelihood == 0, 1, new_likelihood / current_likelihood), 1)
        u = runif(1, 0, 1)
        if (u < ratio) {
            current_R <- new_R
        }
        
        sampled_R[i] <- current_R
    }
    
    res <- sampled_R[(burn_in + 1):n_sims]
    return(res)
}

#' @title Sample resistance
#' 
#' @description MH-sampling from the posterior resistance distribution dependent only on the week. 
#' 
#' @param resistance_parameters A list containing the weeks and a matrix of estimated parameters constructed using \link{estimate_parameters_resistance_weeks}.
#' @param n_sims The number of samples generated by the MH.
#' @param burn_in The number of initial samples removed as burn-in. 
#' @param return_tibble TRUE/FALSE: Should a tibble be returned.
#' @param trace TRUE/FALSE: Should a trace be shown?
#' @param traceLimit Numeric. How often should the trace be shown?
#' @export
resistance_given_week <- function(resistance_parameters, n_sims = 10000, burn_in = 1000, 
                                  return_tibble = TRUE, trace = TRUE, traceLimit = 1000) {
    if ((burn_in + 1) > n_sims) {
        n_sims = burn_in + 1
        warning("'n_sims' argument should be '1' larger than 'burn_in' argument, as all burn-in is discarded.")
    }
    
    gam_fit <- gam_parameter_fit(resistance_parameters)
    W <- resistance_parameters$W[1]:resistance_parameters$W[length(resistance_parameters$W)]
    
    parameters <- predict(gam_fit, newdata = tibble(Week = W))
    parameter_matrix <- do.call("cbind", parameters) %>% as.matrix()
    
    res <- vector("list", length(W))
    for (w in seq_along(W)) {
        beta_w <- parameter_matrix[w, -4]
        sd_w <- parameter_matrix[w, 4]
        
        res[[w]] <- resistance_of_given_week(beta_w, sd_w, n_sims, burn_in, trace, traceLimit, (resistance_parameters$W[1] - 1) + w)
    }
    
    if (return_tibble) {
        res <- res %>% 
            enframe(name = "Week", value = "R") %>% 
            mutate(Week = Week + (resistance_parameters$W[1] - 1)) %>% 
            unnest()
    }
    
    return(res)
}

#' @title Week probability 
#' 
#' @description The probability of each week given a resistance and SOC. 
#' 
#' @param R A new resistance.
#' @param SOC A new state-of-charge.
#' @param resistance_parameters A matrix of estimated parameters constructed using \link{estimate_parameters_resistance_weeks}.
#' 
#' @return A numeric vector of probabilities.
#' @export
exact_week_probablity <- function(R, SOC, resistance_parameters) {
    W <- resistance_parameters$W
    parameter_matrix <- resistance_parameters$Parameters
    
    post_prob <- sapply(1:dim(parameter_matrix)[1], function(w) {
        BatteryPrognostics:::normalised_probability_w(log(R), SOC, w, parameter_matrix)
    })
    
    res <- list(W = W, Posterior = post_prob)
    class(res) <- "post_w"
    return(res)
}

#' @title Week probability 
#' 
#' @description The probability of each week given a resistance and SOC. 
#' 
#' @param R A new resistance.
#' @param resistance_parameters A matrix of estimated parameters constructed using \link{estimate_parameters_resistance_weeks}.
#' @param SOC_dist Distribution of the SOC.
#' @param ... Additional arguments passed to...
#' 
#' @return A numeric vector of probabilities.
#' @export
week_probability <- function(R, resistance_parameters, SOC_dist = dunif, ...) {
    W <- resistance_parameters$W
    parameter_matrix <- resistance_parameters$Parameters
    
    n_cols <- dim(parameter_matrix)[2]
    g <- function(y, v) {
        which.min(abs(y - v))
    }
    
    f <- function(x, w) {
        if (class(SOC_dist) == "function") {
            probability_soc <- SOC_dist(x, ...)
        }
        else if (class(SOC_dist) %in% "histogram") {
            probability_soc <- SOC_dist$density[sapply(x, g, v = SOC_dist$mids)]
        }
        
        beta_ = parameter_matrix[w, -n_cols]
        sd_ = parameter_matrix[w, n_cols]
            
        res <- BatteryPrognostics:::unnormalised_joint_probability(R = log(R), SOC = x, beta_w = beta_, sigma_w = sd_) * 
            probability_soc
        return(res)
    }
    
    #
    N_ = 1001
    x_ <- seq(0 + 1 / N_, 1 - 1 / N_, length.out = N_)
    
    h_  <- x_[2] - x_[1]
    unnormalised <- sapply(seq_along(W), function(w) {
        y_ <- f(x_, w)
        ca <- (y_[2] - y_[1]) / h_
        cb <- (y_[N_] - y_[N_ - 1]) / h_
        
        pracma::trapz(x_, y_) - h_^2 / 12 * (cb - ca)
    })
    
    #
    normalised <- unnormalised / sum(unnormalised)
    res <- list(W = W, Posterior = normalised)
    class(res) <- "post_w"
    return(res)
}

#' Weigthed quantiles of remaining useful life.
#' 
#' @description A wrapper for the \link{weighted.quantile} function of the \link{spatstat} package designed to work with 'post_w' objects.
#' 
#' @param probability_w An object of class 'post_w', created with either the \link{exact_week_probability} or \link{week_probability}.
#' @param ... Additional arguments passed to \link{weighted.quantile}.
#' 
#' @return Numeric vector of quantiles.
#' @export
weighted_quantile_rul <- function(probabilities_w, ...) {
    spatstat::weighted.quantile(x = probabilities_w$W[length(probabilities_w$W)] - probabilities_w$W, 
                                w = probabilities_w$Posterior, 
                                ...)
}

#' @export
hdr <- function(object, ...) {
    UseMethod("hdr")
}

#' @export
hdr.default <- function(object, ...) {
    cat(paste("Function not implemented for class '", class(object), "'. \n", sep = ""))
}

#' @export
hdr.post_w <- function(object, ...) {
    W <- object$W
    normalised <- object$Posterior
    normalised_order <- order(normalised, decreasing = TRUE)
    
    dots_list <- list(...)
    if (is.null(dots_list$p)) {
        p <- 0.95
    } 
    else {
        p <- dots_list$p
    }
    
    res <- vector("list", length(p)) 
    for (i in seq_along(p)) {
        critical_region <- sort(normalised_order[1:which(cumsum(normalised[normalised_order]) >= p[i])[1]])
        region_position <- which(diff(critical_region) != 1)
        
        start_week <- c(critical_region[c(1, region_position + 1)])
        end_week <- c(critical_region[region_position], critical_region[length(critical_region)])
        
        cmf <- sapply(1:length(start_week), function(xx) sum(normalised[start_week[xx]:end_week[xx]]))
        ints <- paste("[", W[start_week], "; ", W[end_week], "]", sep = "")
        
        res[[i]] <- tibble("alpha" = p[i], "HDR" = ints, "Lower" = W[start_week], "Upper" = W[end_week], "Density" = cmf)
    }
    
    post_mode <- normalised[which.max(normalised)]
    names(post_mode) <- paste(W[[which.max(normalised)]])
    res_list <- list(mode = post_mode, HDR = res)
    return(res_list)
}

#' @title Mean residual life
#' 
#' @description Calculates the mean residual life, given the current age of the battery, and a 'post_w' class objcet.
#' 
#' @param w The current age of the battery.
#' @param probability_w An object of class 'post_w', created with either the \link{exact_week_probability} or \link{week_probability}.
#' 
#' @return Numeric.
#' @export
mean_residual_life <- function(w, probability_w) {
    which_week <- which(probability_w$W == w)
    S <- 1 - cumsum(probability_w)
    mrl <- sum(S[which_week:length(S)]) / S[which_week]
    return(mrl)
}

#' @title Expected mean residual life
#' 
#' @description Calculates the expected mean residual life of a 'post_w' class objcet.
#' 
#' @param probability_w An object of class 'post_w', created with either the \link{exact_week_probability} or \link{week_probability}.
#' 
#' @return Numeric.
#' @export
expected_residual_life <- function(post_w) {
    erl <- sum(1 - cumsum(post_w))
    return(erl)
}

#' @title Standard deviation of mean residual life
#' 
#' @description Calculates the standard deviation of the mean residual life of a 'post_w' class objcet.
#' 
#' @param probability_w An object of class 'post_w', created with either the \link{exact_week_probability} or \link{week_probability}.
#' 
#' @return Numeric.
#' @export
sd_residual_life <- function(post_w) {
    sdrl <- (2 * sum((post_w$W) * (1 - cumsum(post_w))) - sum(1 - cumsum(post_w))^2)^(1 / 2)
    return(sdrl)
}

#' @export
ggplot.post_w <- function(object, ...) {
    arg_list <- list(...)
    
    if (is.null(arg_list$include_hdr)) {
        include_hdr = FALSE
    }
    else {
        include_hdr = arg_list$include_hdr
    }
    
    if (is.null(arg_list$rul)) {
        rul = FALSE
    }
    else {
        rul = arg_list$rul
    }
    
    obj <- object %>% as_tibble()
    obj_plot <- ggplot(obj, aes_(x = ~Week, y = ~Posterior)) + 
        geom_line() + geom_point() + 
        ylab("Posterior probability") + 
        theme_bw()
    
    if (!is.null(arg_list$base_size)) {
        obj_plot <- obj_plot + theme_bw(base_size = arg_list$base_size)
    }
    
    if (include_hdr) {
        if (is.null(arg_list$p)) {
            p <- 0.95
        }
        else {
            p <- arg_list$p
            
            if (length(p) > 1) {
                p <- p[1]
                
                warning("Only the first value of 'p' is used.")
            }
        }
        
        hdr_w <- hdr(object, p = p)
        
        hdr_1 <- hdr_w$HDR[[1]]
        for (i in 1:dim(hdr_1)[1]) {
            obj_w <- obj %>% filter(Week >= hdr_1$Lower[i], Week <= hdr_1$Upper[i])
            
            obj_plot <- obj_plot + geom_ribbon(data = obj_w, 
                                               aes(ymin = 0, ymax = Posterior), fill = "blue", alpha = "0.5", colour = "transparent")
        }
        
        obj_plot <- obj_plot + ggtitle(paste(100 * p, "% High posterior density region", sep = ""))
    }
    
    if (rul) {
        week_max <- max(obj$Week)
        step_size <- ceil(week_max / 10)
        obj_plot <- obj_plot + 
            xlab("Remaining useful life") + 
            scale_x_continuous(breaks = seq(1, week_max, step_size), labels = seq(week_max - 1, 0, -step_size))
    }
    
    return(obj_plot)
}

#' @title SOC probability 
#' 
#' @description The marginal probability of SOC given a resistance and week. 
#' 
#' @param R The extracted battery resistance.
#' @param W The corresponding week of the extracted battery resistance.
#' @param resistance_parameters A list containing the weeks and a matrix of estimated parameters constructed using \link{estimate_parameters_resistance_weeks}.
#' @param n_sims The number of samples generated by the MH.
#' @param burn_in The number of initial samples removed as burn-in. 
#' @param trace TRUE/FALSE: Should a trace be shown?
#' @param traceLimit Numeric: How often should the trace be shown?
#' 
#' @return A numeric vector of probabilities.
#' @export
posterior_soc_probability <- function(R, W, resistance_parameters, burn_in, n_sims, trace, traceLimit) {
    parameter_matrix <- resistance_parameters$Parameters
    sd <- mean(resistance_parameters$Parameters[, "sd"], na.rm = TRUE)
    table_W <- table(W)
    which_non_zero <- which(table_W > 0)
    
    beta_0 <- rep(parameter_matrix[which_non_zero, "beta_0"], times = table_W[which_non_zero])
    beta_1 <- rep(parameter_matrix[which_non_zero, "beta_1"], times = table_W[which_non_zero])
    
    beta_2 <- NULL
    if ("beta_2" %in% colnames(parameter_matrix)) {
        beta_2 <- rep(parameter_matrix[which_non_zero, "beta_2"], times = table_W[which_non_zero])
    }
    
    current_soc <- runif(1, 0, 1)
    current_likelihood <- sum(dnorm(log(R), 
                                mean = BatteryPrognostics:::resistance_expected_vectorised(beta_0, beta_1, beta_2, current_soc), 
                                sd = sd, log = T))
    
    sampled_soc <- rep(NA, n_sims)
    sampled_soc[1] <- current_soc
    for (i in 2:n_sims) {
        if (trace & ((i %% traceLimit) == 0))
            cat("Iteration:", i, "/", sprintf("%i", n_sims), "\n")
        
        new_soc <- runif(1, 0, 1) # runif(1, max(current_soc - 0.1, 0), min(current_soc + 0.1, 1))
        new_likelihood <- sum(dnorm(log(R), 
                                    mean = BatteryPrognostics:::resistance_expected_vectorised(beta_0, beta_1, beta_2, new_soc), 
                                    sd = sd, log = T))
        
        ratio <- min(ifelse(current_likelihood > 0, new_likelihood / current_likelihood, 1), 1)
        u = runif(1, 0, 1)
        if (u < ratio) {
            current_soc <- new_soc
        }
        
        sampled_soc[i] <- current_soc
    }
    
    post_prob <- sampled_soc[burn_in:n_sims]
    class(post_prob) <- c("post_soc", class(post_prob))
    return(post_prob)
}

#' @title Predict resistance
#' 
#' @description The posterior of the battery resistance
#' 
#' @param post_pred_var An object of class ''.
#' @param resistance_parameters An object of class 'resistance_parameters'.
#' @param SOC A single value, an object of class 'histogram', an object of class '', or a 'function'.
#' @param R The extracted battery resistance.
#' @param W The corresponding week of the extracted battery resistance.
#' @param trace TRUE/FALSE: Show trace?
#' @param traceLimit Limits the amount the trace is printed.
#' @param equal_variance TRUE/FALSE: Is the variance assumed to be equal.
#' @param return_tibble TRUE/FALSE: Should a tibble be returned? If 'FALSE', a matrix is returned.
#' @param ... Additional information passed to 'SOC'.
#'  
#' @return A matrix, or tibble, containing predictions of every simulation for every week. 
#' @export
posterior_resistance <- function(post_pred_var, resistance_parameters,
                                 SOC = NULL, R = NULL, W = NULL,
                                 trace = FALSE, traceLimit = 1000,
                                 equal_variance = FALSE, return_tibble = FALSE,
                                 ...) {
    dots <- list(...)
    weeks <- resistance_parameters$W
    
    K = dim(post_pred_var$Simulations[[1]])[2]
    n_ahead <- dim(post_pred_var$Simulations[[1]])[1]
    n_sims <- length(post_pred_var$Simulations)
    
    if (is.null(SOC)) {
        if (is.null(R) || is.null(W)) {
            stop("'SOC' was NULL and 'R' and 'W' were not provided.")
        }
        
        if (length(R) != length(W)) {
            stop("'SOC' was NULL and 'R' and 'W' were of different length.")
        } 
        
        if (is.null(burn_in)) {
            burn_in = 1000
        }
        else {
            burn_in <- dots$burn_in
            dots <- dots[-which(names(dots) == "burn_in")]
        }
        
        cat("Simulating posterior of 'SOC'.\n")
        post_SOC <- posterior_soc_probability(R = R, W = W, resistance_parameters = resistance_parameters,
                                              burn_in = burn_in, n_sims = n_sims, trace = trace, 
                                              traceLimit = traceLimit)
        
        hist_SOC <- hist(post_SOC, breaks = "fd", plot = FALSE)
        SOC_ <- sample_histogram(hist_SOC, n_sims * (n_ahead + length(weeks)))
    }
    else if ("post_soc" %in% class(SOC)) {
        hist_SOC <- hist(SOC, breaks = "fd", plot = FALSE)
        SOC_ <- sample_histogram(hist_SOC, n_sims * (n_ahead + length(weeks)))
    }
    else if('histogram' %in% class(SOC)) {
        hist_SOC <- SOC
        SOC_ <- sample_histogram(hist_SOC, n_sims * (n_ahead + length(weeks)))
    }
    else if (is.numeric(SOC) & ((length(SOC) == 1) || (length(SOC) == n_sims * (n_ahead + length(weeks))))) {
        SOC_ <- SOC
    }
    else if (is.function(SOC)) {
        dots$n <- n_sims * n_ahead
        SOC_ <- do.call(SOC, dots)
    }
    
    sd_ <- NULL
    if (equal_variance || ((equal_variance == "bartlett") & (!as.logical(resistance_parameters$BarlettsTest["Rejected"])))) {
        sd_ <- rep(mean(resistance_parameters$Parameters[, "sd"]), n_ahead + length(weeks))
    }
    
    estimated_parameters <- resistance_parameters$Parameters
    parameter_simulations <- post_pred_var$Simulations
    par_names <- colnames(parameter_simulations[[1]])
    res <- matrix(NA, nrow = n_sims, ncol = length(weeks) + n_ahead)
    for (i in 1:n_sims) {
        parameter_simulations_i <- rbind(estimated_parameters[, par_names], parameter_simulations[[i]])
        sd_i <- sd_
        if (is.null(sd_i) & ("sd" %in% par_names)) {
            sd_i <- parameter_simulations_i[, "sd"]
        }
        
        beta_0_i <- parameter_simulations_i[, "beta_0"]
        beta_1_i <- parameter_simulations_i[, "beta_1"]
        
        beta_2_i <- NULL
        if ("beta_2" %in% par_names) {
            beta_2_i <- parameter_simulations_i[, "beta_2"]
        }
        
        SOC_i <- SOC_[((i - 1) * (n_ahead + length(weeks)) + 1):(i * (n_ahead + length(weeks)))]
        expected_log_r <- BatteryPrognostics:::resistance_expected_vectorised(beta_0 = beta_0_i, beta_1 = beta_1_i, 
                                                                              beta_2 = beta_2_i, SOC = SOC_i)
        res[i, ] <- rnorm(n_ahead + length(weeks), mean = expected_log_r, sd = sd_i)
    }
    
    if (return_tibble) {
        res_ <- t(res) %>% 
            as_tibble() %>% 
            mutate(Week = (weeks[1]):n()) %>% 
            gather(Simulation, R, -Week) %>% 
            mutate(Simulation = as.numeric(stringr::str_sub(Simulation, start = 2)))
    }
    else {
        res_ <- res
    }
    
    class(res_) <- c("post_r", class(res_))
    return(res_)
}

