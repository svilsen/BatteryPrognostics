#### (V)AR Feature modelling ----

#' @title Estimate feature AR models
#' 
#' @description Wrapper for esay estimation of feature AR models in RUL estimation.
#' 
#' @param data A data-set with each column containing the features to be modelled.
#' @param l The lag of the AR feature models (NOTE: predict and forward only works with '\code{l = 1}'). 
#' 
#' @return A list of AR models.
#' @export
feature_vars <- function(data, l = 1) {
    N <- dim(data)[1]
    M <- dim(data)[2]
    
    res <- vector("list", M)
    for (m in seq_along(res)) {
        Y_m <- data[, m]
        data_m <- data.frame(
            Y = unname(unlist(Y_m[seq(1 + l, N), ]))
        )
        
        for (ll in seq_len(l)) {
            data_m$Z <- unname(unlist(Y_m[seq(ll, N - l + ll - 1), ]))
            colnames(data_m)[colnames(data_m) == "Z"] <- paste0("Y", ll)
        }
        
        data_m <- data_m %>% 
            as_tibble() %>% 
            mutate(Index = seq_len(N - ll) - 1)
        
        res[[m]] <- lm(Y ~ 1 + ., data = data_m)
    }
    
    names(res) <- colnames(data)
    class(res) <- "featurevars"
    return(res)
}

simulate_ar <- function(object, s, nr_simulations) {
    coefs <- coef(object)
    residual <- resid(object)
    sigma <- sqrt(sum((residual)^2) / (length(residual) - length(coefs))) 
    T_end <- length(residual)
    
    res <- vector("list", nr_simulations)
    for (r in seq_len(nr_simulations)) {
        x_new <- matrix(s, ncol = 1, nrow = T_end + 1)
        for (k in seq_len(T_end)) {
            x_new[k + 1, ] <- coefs[1] + coefs[2] * x_new[k, ] + coefs[3] * (k - 1) + rnorm(1, 0, sigma)
        }
        
        colnames(x_new) <- names(s)
        res[[r]] <- x_new
    }
    
    return(res)
}

#' @export 
predict.featurevars <- function(object, ...) {
    dots <- list(...)
    if (is.null(dots$nr_simulations) || !is.numeric(dots$nr_simulations)) {
        nr_simulations <- 100
    }
    else {
        nr_simulations <- dots$nr_simulations
    }
    
    res <- lapply(seq_along(object), function(i) {
        object_i <- object[[i]]
        
        if (is.null(dots$s) || (!is.numeric(dots$s)) || ((length(dots$s) != length(object)) && (length(dots$s) != 1))) {
            s_i <- object_i$model[1, max(which(str_detect(names(object_i$model), "Y")))]
        } 
        else if (length(dots$s) == 1) {
            s_i <- dots$s[1] 
        } 
        else {
            s_i <- dots$s[i]
        }
        
        rr <- simulate_ar(object = object_i, s = s_i, nr_simulations = nr_simulations) %>% 
            enframe(name = "Simulation", value = "VALUE") %>% 
            unnest(cols = "VALUE") %>% 
            mutate(Index = rep(seq_len(n() / nr_simulations), nr_simulations), Feature = names(object[i])) %>% 
            select(Simulation, Index, Feature, VALUE) %>% 
            mutate(VALUE = VALUE[, 1])
        
        return(rr)
    })
    
    return(res)
}

#### SOH forward simulation ----
simulate_forward_ar <- function(object, s, nr_simulations, nr_forward) {
    coefs <- coef(object)
    residual <- resid(object)
    sigma <- sqrt(sum((residual)^2) / (length(residual) - length(coefs))) 
    T_end <- length(residual)
    
    res <- vector("list", nr_simulations)
    for (r in seq_len(nr_simulations)) {
        x_new <- matrix(s, ncol = 1, nrow = nr_forward + 1)
        for (k in seq_len(nr_forward)) {
            x_new[k + 1, ] <- coefs[1] + coefs[2] * x_new[k, ] + coefs[3] * (T_end + k - 1) + sqrt(k) * rnorm(1, 0, sigma)
        }
        
        colnames(x_new) <- names(s)
        res[[r]] <- x_new
    }
    
    return(res)
}

#' @title Forward simulation
#' 
#' @description \code{forward} is a generic function for the look-ahead simulation of various stochastic process models. This function will invoke a particular method which will depend on the \link{class} of the first argument.  
#' 
#' @param object An object to simulate forward in time.
#' @param ... Additional arguments needed in forward simulation.
#' 
#' @export
forward <- function(object, ...) {
    UseMethod("forward")
}

#' @export
forward.default <- function(object, ...) {
    cat(paste("Function not implemented for class '", class(object), "'. \n", sep = ""))
}

#' @export
forward.featurevars <- function(object, ...) {
    dots <- list(...)
    if (is.null(dots$nr_simulations) || !is.numeric(dots$nr_simulations)) {
        nr_simulations <- 100
    }
    else {
        nr_simulations <- dots$nr_simulations
    }
    
    if (is.null(dots$nr_ahead) || !is.numeric(dots$nr_ahead)) {
        nr_ahead <- 10
    }
    else {
        nr_ahead <- dots$nr_ahead
    }
    
    if (is.null(dots$last_index) || !is.numeric(dots$last_index)) {
        last_index <- 1
    }
    else {
        last_index <- dots$last_index
    }
    
    res <- lapply(seq_along(object), function(i) {
        object_i <- object[[i]]
        
        if (is.null(dots$s) || (!is.numeric(dots$s)) || ((length(dots$s) != length(object)) && (length(dots$s) != 1))) {
            s_i <- object_i$model[1, max(which(str_detect(names(object_i$model), "Y")))]
        } else if (length(dots$s) == 1) {
            s_i <- dots$s[1] 
        } else {
            s_i <- dots$s[i]
        }
        
        rr <- simulate_forward_ar(object = object_i, s = s_i, nr_simulations = nr_simulations, nr_forward = nr_ahead) %>% 
            enframe(name = "Simulation", value = "VALUE") %>% 
            unnest(cols = "VALUE") %>% 
            mutate(
                Index = rep(seq_len(n() / nr_simulations), nr_simulations) + last_index - 1,
                Feature = names(object[i])
            ) %>% 
            select(Simulation, Index, Feature, VALUE) %>% 
            mutate(VALUE = VALUE[, 1])
        
        return(rr)
    })
    
    return(res)
}

#### RUL estimation ----
#' @title Simulating Remaining Useful Life
#' 
#' @description The set-up, simulation, and prediction of remaining useful life (RUL) based on a feature abstraction approach.
#' 
#' @param X A matrix of features.
#' @param y A vector of health metrics.
#' @param health_model A function modelling the relationship between features and health metrics. 
#' @param n_ahead Amount of time simulated forward.
#' @param n_simulations The number of simulations used in the forward simulation.
#' @param l The lag used in feature modelling.
#' @param eol A numeric representing the end-of-life criterion.
#' @param y0 An initial value of the health metric used to calculate the end-of-life. If '\code{NULL}' (default) it uses the first value of '\code{y}'.
#' @param eol_method The method used to determine end-of-life. Takes a string \code{"first"}, \code{"last"}, or \code{"quantile"}.
#' @param return_full \code{TRUE/FALSE}: Should all estimated and simulated objects be returned?
#' @param ... A list of parameters passed to the \code{health_model} function. 
#' 
#' @return A list containing the following four elements are returned: 
#' 
#' @export
rul <- function(X, y, health_model, n_ahead = 50, n_simulations = 100, l = 1, eol = 0.8, y0 = NULL, eol_method = "first", return_full = FALSE, ...) {
    ## Set-up
    health_model_parameters <- list(...)
    if ((!is.numeric(n_ahead)) || (abs(n_ahead - round(n_ahead)) > 1e-8) || (n_ahead < 1)) {
        stop("'n_ahead' has to be a numeric, or an interger, larger than one.")
    }
    
    if ((!is.numeric(n_simulations)) || (abs(n_simulations - round(n_simulations)) > 1e-8) || (n_simulations < 1)) {
        stop("'n_simulations' has to be a numeric, or an interger, larger than one.")
    }
    
    if ((!is.numeric(l)) || (abs(l - round(l)) > 1e-8) || (l < 1)) {
        stop("'l' has to be a numeric, or an interger, larger than one.")
    }
    
    if (is.null(y0)) {
        y0 <- unname(unlist(y)[1])
    }
    
    if (is.null(return_full)) {
        return_full <- FALSE
    }
    
    if (!is.logical(return_full)) {
        stop("'return_full' has to be TRUE or FALSE.")
    }
    
    ## Estimate feature model
    estimated_feature_models <- feature_vars(as_tibble(X), l)
    expected_features <- predict(estimated_feature_models, s = NULL, nr_simulations = n_simulations) %>% 
        bind_rows() %>% 
        group_by(Index, Feature) %>% 
        summarise(VALUE = mean(VALUE), .groups = "drop") %>% 
        mutate(Feature = factor(Feature, levels = colnames(X))) %>% 
        arrange(Feature, Index) %>% 
        pivot_wider(names_from = "Feature", values_from = "VALUE") %>% 
        select(-Index)
    
    ## Estimate SOH model
    if (any(str_detect(names(formals(health_model)), "formula"))) {
        health_model_parameters$formula <- y ~ .
        health_model_parameters$data <- cbind(X, y = unlist(y))
        
        estimated_health_model <- suppressWarnings(do.call(health_model, health_model_parameters))
    } else {
        health_model_parameters$X <- X
        health_model_parameters$y <- y
        
        if (any(str_detect(names(formals(health_model)), "x"))) {
            names(health_model_parameters)[which(names(health_model_parameters) == "X")] <- "x"
        }
        
        if (any(str_detect(names(formals(health_model)), "Y"))) {
            names(health_model_parameters)[which(names(health_model_parameters) == "y")] <- "Y"
        }
        
        estimated_health_model <- suppressWarnings(do.call(health_model, health_model_parameters))
    }
    
    ## Forward simulation of feature model
    X_endpoint <- unlist(X[dim(X)[1], ])
    last_index <- dim(y)[1]
    forward_features <- forward(estimated_feature_models, s = X_endpoint, nr_simulations = n_simulations, nr_ahead = n_ahead, last_index = last_index) %>% 
        bind_rows() %>% 
        pivot_wider(names_from = "Feature", values_from = "VALUE") 
    
    X_new <- forward_features %>%
        split(., f = .$Simulation) %>% 
        lapply(., function(xx) xx %>% select(-Simulation, -Index))
    
    ## Forward simulation of SOH model
    forward_soh <- lapply(X_new, function(xn) suppressWarnings(predict(estimated_health_model, newdata = xn)))
    
    ## RUL estimation
    eol_value <- eol * y0
    if (tolower(eol_method) == "first") {
        rul <- sapply(forward_soh, function(xx) {
            min(which(xx < eol_value))
        }) %>% enframe(name = "Simulation", value = "RUL")
    } else if (tolower(eol_method) == "last") {
        rul <- sapply(forward_soh, function(xx) {
            max(which(xx > eol_value))
        }) %>% enframe(name = "Simulation", value = "RUL")
    } else if (tolower(eol_method) == "quantile") {
        rul <- forward_soh %>% 
            enframe(name = "SIM", value = "VAL") %>% 
            unnest(cols = "VAL") %>% 
            group_by(SIM) %>% 
            mutate(IND = seq_len(n()), VAL = VAL[, 1]) %>% 
            ungroup() %>% 
            group_by(IND) %>% 
            summarise(QLower = quantile(VAL, probs = 0.025), 
                      Q = quantile(VAL, probs = 0.5), 
                      QUpper = quantile(VAL, probs = 0.975), 
                      .groups = "drop") %>% 
            mutate(RUL = n_ahead - IND + 1) %>% 
            select(RUL, QLower, Q, QUpper)
    } else {
        stop("'eol_method' not implemented.")
    }
    
    res <- list(
        RUL = rul,
        EOL = list(Percentage = eol, Value = eol_value), 
        Objects = list()
    )
    
    if (return_full) {
        res$Objects <- list(
            FeatureModel = estimated_feature_models, 
            HealthModel = estimated_health_model, 
            FeatureForward = forward_features, 
            HealthForward = forward_soh
        )
    }
    
    class(res) <- "rul"
    return(res)
}

#' @export
mean.rul <- function(x, ...) {
    mean(x$RUL$RUL, ...)
}

#' @export
sd.rul <- function(x, na.rm = FALSE) {
    sd(x$RUL$RUL, na.rm = na.rm)
}

#' @export
median.rul <- function(x, na.rm = FALSE, ...) {
    median(x$RUL$RUL, na.rm = na.rm, ...)
}

#' @export
quantile.rul <- function(x, ...) {
    quantile(x$RUL$RUL, ...)
}

#' @export 
hist.rul <- function(x, ...) {
    hist(x = x$RUL$RUL, xlab = "RUL", main = "Histogram of RUL", ...)
}

#' @export
ggplot.rul <- function(rul_object, ...) {
    dots <- list(...)
    
    if (is.null(dots$bins)) {
        dots$bins <- 10
    }
    
    if (is.null(dots$errorbar)) {
        dots$errorbar <- TRUE
    }
    
    if (is.null(dots$errortext)) {
        dots$errortext <- TRUE
    }
    
    if (is.null(dots$level)) {
        dots$level <- 0.95
    }
    
    if (is.null(dots$period)) {
        dots$period <- 1
    }
    
    rul_object$RUL <- rul_object$RUL %>% mutate(RUL = RUL / dots$period)
    p <- ggplot(rul_object$RUL, aes(x = RUL))
    if (is.null(dots$density) || (!dots$density)) {
        p <- p + 
            geom_histogram(bins = dots$bins, colour = "dodgerblue2") + 
            ylab("Count")
    } else {
        p <- p + 
            geom_histogram(aes(y = ..density..), bins = dots$bins, colour = "dodgerblue2") + 
            ylab("Density")
    }
    
    if (dots$errorbar) {
        plot_max <- ifelse(is.null(dots$density) || (!dots$density), max(ggplot_build(p)$data[[1]]$count), max(ggplot_build(p)$data[[1]]$density))
        width_error <- 0.05 * plot_max
        rul_object_summary <- rul_object$RUL %>% 
            summarise(Lower = quantile(RUL, (1 - dots$level) / 2), 
                      Median = quantile(RUL, 0.5), 
                      Upper = quantile(RUL, 1 - (1 - dots$level) / 2))
        
        p <- p + 
            geom_errorbar(
                data = rul_object_summary, 
                aes(x = NULL, y = 0, xmin = Lower, xmax = Upper), 
                width = width_error,
                position = position_dodge(0.05), 
                colour = "red", 
                size = 0.8
            ) + 
            geom_point(
                data = rul_object_summary, 
                aes(x = Median, y = 0), 
                colour = "red", 
                size = 3
            )
        
        if (dots$errortext) {
            text_height <- -0.05 * plot_max
            text_lower <- grid::textGrob(paste(round(rul_object_summary$Lower, 2)), gp = grid::gpar(fontsize = 13, col = "red"))
            text_median <- grid::textGrob(paste(round(rul_object_summary$Median, 2)), gp = grid::gpar(fontsize = 13, col = "red"))
            text_upper <- grid::textGrob(paste(round(rul_object_summary$Upper, 2)), gp = grid::gpar(fontsize = 13, col = "red"))
            p <- p + 
                annotation_custom(text_lower, xmin = rul_object_summary$Lower, xmax = rul_object_summary$Lower, ymin = text_height, ymax = text_height) + 
                annotation_custom(text_median, xmin = rul_object_summary$Median, xmax = rul_object_summary$Median, ymin = text_height, ymax = text_height) + 
                annotation_custom(text_upper, xmin = rul_object_summary$Upper, xmax = rul_object_summary$Upper, ymin = text_height, ymax = text_height) + 
                coord_cartesian(clip = "off")
        }
    }
    
    p <- p +
        theme_minimal() + 
        ggtitle("Histogram of RUL")
    
    p
}




