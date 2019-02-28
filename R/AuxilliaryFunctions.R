mode_density <- function(x, na.rm = FALSE, ...) {
    if(na.rm){
        x = x[!is.na(x)]
    }
    
    y <- density(x, ...)
    res <- y$x[which.max(y$y)]
    return(res)
}

#' @title Extract resistance
#'
#' @description Extracts the resistance from a dynamic current/voltage profile.
#'
#' @param I Vector containing the current.
#' @param V Vector containing the voltage.
#' @param T_s Vector containing the time step sizes (in seconds). If 'NULL' it is set as 'seq_along(I)'.
#' @param epsilon The lower limit used for identifying change in current. Default is '1e-3'.
#' @param Q_max The maximal capacity of the battery cell.
#' @param eta The Coulombic effeciency.
#' @param SOC_0 The initial value of the state-of-charge.
#' @param convert_to_tibble Organise in \link{tibble}.
#'
#' @return A list containing the extracted resistance and their time scales extracted from the profile.
#' @export
extract_resistance <- function(I, V, T_s = NULL, epsilon = 1e-3, Q_max = 2.56, eta = 1.0, SOC_0 = 0.0, 
                               convert_to_tibble = FALSE) { 
    if (is.null(I) | is.null(V)) {
        stop("Both 'I' and 'V' has to be provided.")
    }
    
    if (is.null(T_s)) {
        T_s = seq_along(I)
    }
    
    if (is.null(epsilon)) {
        epsilon = 0.5
    }
    
    res_matrix <- BatteryPrognostics:::extract_resistance_cpp(I, V, T_s, epsilon, Q_max, eta, SOC_0) 
    class(res_matrix) <- "ext_res"
    
    if (convert_to_tibble) {
        res_tibble <- as_tibble(res_matrix)
        data_tibble <- tibble("I" = I, "V" = V, "T_s" = T_s) %>% 
            bind_cols(res_tibble) %>% 
            group_by(ID) %>% 
            mutate("SL" = max(S)) %>% 
            ungroup() %>% 
            mutate("Day" = floor(T_s / 3600 / 24)) %>% 
            group_by(ID, NonZero) %>% 
            mutate(Relax = 1:n()) %>% 
            ungroup() %>% 
            mutate(Relax = ifelse(NonZero > -1, -1, Relax))
        
        data_tibble$Relax[which(data_tibble$S == 0)] <- 
            data_tibble$Relax[which(data_tibble$S == 0) - 1] 
        
        data_tibble$SLP <- data_tibble$SL
        data_tibble$SLP[which(data_tibble$S == 0)] <- 
            data_tibble$SLP[which(data_tibble$S == 0) - 1] 
        
        data_tibble <- data_tibble %>% 
            group_by(ID, NonZero) %>% 
            mutate(Relax = ifelse(NonZero > -1, Relax[1], Relax), 
                   SLP = ifelse(NonZero > -1, SLP[1], SLP)) %>% 
            ungroup()
        
        return(data_tibble)
    }
    else {
        return(res_matrix)
    }
}

#' @export
as_tibble.ext_res <- function(x, ...) {
    res_tibble <-x[1:2] %>% enframe() %>%
        unnest() %>%
        group_by(name) %>%
        mutate(TT = 1:n()) %>%
        spread(key = "name", value = "value") %>%
        select(-TT) %>%
        mutate(SOC = x$SOC,
               ID = x$ID, 
               NonZero = x$NonZero, 
               Reset = x$Reset) 
    
    return(res_tibble)
}

#' @title 2d-density plot of the extracted resistance
#'
#' @description A wrapping function calling ggplot with 'stat_density_2d' and 'ggMarginal' (from the 'ggExtra'-package).
#'
#' @param ext_res_tibble \link{tibble} made with the \link{extract_resistance}-function.
#' @param restrict_dim_y Quantiles used to restrict the y-axis. See details.
#' @param include_marginals TRUE/FALSE: Should marginal plots be included? Note: requires 'ggExtra'.
#' @param include_current TRUE/FALSE: Include a histogram of the current? Note requires 'gridExtra'.
#' @param dims The number of dimensions to be plotted.
#' @param dim_labels Labels for the 'x' and 'y' dimensions.
#' @param facet_variable String passed to 'facet_wrap'.
#' @param plot_type Character of the type of plot. Takes the values 'density', 'hex', and 'bin'.
#' @param n_bins The number of bins passed to the 'hex' and 'bin' plot types.
#'
#' @details If 'restrict_dim_y' has length 1 it is taken as a restriction on the upper quantile. If it is 'NULL' the y-axis is not restricted.
#'
#' @return ggplot-object.
#' @export
plot_extracted_resistance <- function(ext_res_tibble, restrict_dim_y = NULL, include_marginals = TRUE, include_current = TRUE,
                                      dims = c("SOC", "R"), dim_labels = c(expression("SOC [%]"), expression(paste("R", " [", Omega, "]"))),
                                      facet_variable = NULL, plot_type = "hex", n_bins = NULL) {
    if (length(dims) == 0) {
        dims <- c("SOC", "R")
    }
    
    if ((!is.character(dims) & !is.numeric(dims)) | (length(dims) != 2)) {
        stop("'dims' must be a character or numeric vector of length '2'.")
    }
    
    if(is.numeric(dims)) {
        dims <- names(ext_res_tibble)[dims]
    }
    
    ext_res_tibble_ <- ext_res_tibble
    if (!is.null(restrict_dim_y)) {
        if (length(restrict_dim_y) == 1) {
            ext_res_tibble_ <- ext_res_tibble_ %>% 
                mutate(Upper = quantile(!!sym(dims[2]), probs = restrict_dim_y)) %>% 
                filter((!!sym(dims[2])) <= Upper) %>% 
                dplyr::select(-Upper)
        }
        else if (length(restrict_dim_y) == 2) {
            ext_res_tibble_ <- ext_res_tibble_ %>% 
                mutate(Lower = quantile(!!sym(dims[2]), probs = restrict_dim_y[1]), 
                       Upper = quantile(!!sym(dims[2]), probs = restrict_dim_y[2])) %>% 
                filter((!!sym(dims[2])) >= Lower, (!!sym(dims[2])) <= Upper) %>% 
                dplyr::select(-Lower, -Upper)
        }
        else {
            stop("'restrict_dim_x' should have length '1' or '2', or be 'NULL'.")
        }
    }
    
    legend_position = "right"
    if (include_marginals) 
        legend_position = "bottom"
    
    if (is.null(n_bins)) {
        upper_limit = 25
        if (plot_type == "density") {
            upper_limit = 75
        }
        
        n_bins <- min(round(dim(ext_res_tibble_)[1] / 10, -1), upper_limit)
    }
    
    if (plot_type == "density") {
        p <- ggplot(ext_res_tibble_, aes(x = !!sym(dims[1]), y = !!sym(dims[2]))) + 
            stat_density2d(aes(fill = ..level..), geom = "polygon", n = n_bins) +
            scale_fill_continuous(name = "Count") + 
            theme_bw() + 
            theme(legend.position = legend_position, 
                  panel.background = element_rect(fill = "#0f2134"), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
    }
    else if (plot_type == "hex") {
        p <- ggplot(ext_res_tibble_, aes(x = !!sym(dims[1]), y = !!sym(dims[2]))) + 
            stat_bin_hex(bins = n_bins) +
            scale_fill_continuous(name = "Count") +
            theme_bw() + 
            theme(legend.position = legend_position, 
                  panel.background = element_rect(fill = "#0f2134"), #132B43 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
    }
    else if (plot_type == "bin") {
        p <- ggplot(ext_res_tibble_, aes(x = !!sym(dims[1]), y = !!sym(dims[2]))) + 
            stat_bin_2d(bins = n_bins) +
            scale_fill_continuous(name = "Count") +
            theme_bw() + 
            theme(legend.position = legend_position, 
                  panel.background = element_rect(fill = "#0f2134"), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())
    }
    else {
        stop("'plot_type' should take the value 'density', 'hex', or 'bin'.")
    }
    
    if (length(dim_labels) == 1) {
        p <- p + ylab(dim_labels)
    } 
    else if (length(dim_labels) == 2) {
        p <- p + xlab(dim_labels[1]) + ylab(dim_labels[2])
    }
    
    if (include_marginals & !is.null(facet_variable)) {
        warning("When 'include_marginals' is 'TRUE' and 'facet_variable' isn't 'NULL' only the marginals are plotted.")
    }
    
    if (include_marginals) {
        p <- ggMarginal(p, type = "density",
                        xparams = list(fill = "grey"),
                        yparams = list(fill = "grey"))
    }
    else if (!is.null(facet_variable)) {
        p <- p + facet_wrap(as.formula(paste("~", facet_variable)), nrow = 1)
    }
    
    if (include_current) {
        p_current <- ggplot(ext_res_tibble_, aes(x = I)) + geom_density(aes(y = ..density..), fill = "grey", alpha = 0.8) + 
            xlab("Current [A]") + ylab("Density") + theme_bw() + removeGrid()
        
        if (include_marginals) {
            p_current <- p_current + theme(plot.margin = unit(c(60, 5.5, 60, 10), "points"))
        }
        
        grid.arrange(p, p_current, ncol = 2, nrow = 1, widths = c(0.6, 0.4))
    }
    else {
        return(p)
    }
}

##
rwish_ <- function (v, S) {
    if (!is.matrix(S)) 
        S <- matrix(S)
    if (nrow(S) != ncol(S)) {
        stop(message = "S not square in rwish().\n")
    }
    if (v < nrow(S)) {
        stop(message = "v is less than the dimension of S in rwish().\n")
    }
    p <- nrow(S)
    CC <- chol(S)
    Z <- matrix(0, p, p)
    diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
    if (p > 1) {
        pseq <- 1:(p - 1)
        Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * (p - 1)/2)
    }
    return(crossprod(Z %*% CC))
}

riwish_ <- function (v, S) 
{
    return(solve(rwish_(v, solve(S))))
}

##
is_square_matrix <- function(x) {
    return(nrow(x) == ncol(x))
}

is_symmetric_matrix <- function(x) {
    return(all((x - t(x)) < (2 * (.Machine$double.eps / 2)^(1/3))))
}

is_positive_definite <- function(x) {
    if (!is.matrix(x)) 
        stop("x is not a matrix.")
    if (!is_square_matrix(x)) 
        stop("x is not a square matrix.")
    if (!is_symmetric_matrix(x)) 
        stop("x is not a symmetric matrix.")
    eigs <- eigen(x, symmetric = TRUE)$values
    if (any(is.complex(eigs))) 
        return(FALSE)
    if (all(eigs > 0)) 
        pd <- TRUE
    else pd <- FALSE
    return(pd)
}

rmatrixnorm_ <- function (M, U, V) {
    if (missing(M)) 
        stop("Matrix M is missing.")
    if (!is.matrix(M)) 
        M <- matrix(M)
    if (missing(U)) 
        stop("Matrix U is missing.")
    if (!is.matrix(U)) 
        U <- matrix(U)
    if (!is_positive_definite(U)) 
        stop("Matrix U is not positive-definite.")
    if (missing(V)) 
        stop("Matrix V is missing.")
    if (!is.matrix(V)) 
        V <- matrix(V)
    if (!is_positive_definite(V)) 
        stop("Matrix V is not positive-definite.")
    if (nrow(M) != nrow(U)) 
        stop("Dimensions of M and U are incorrect.")
    if (ncol(M) != ncol(V)) 
        stop("Dimensions of M and V are incorrect.")
    n <- nrow(U)
    k <- ncol(V)
    Z <- matrix(rnorm(n * k), n, k)
    X <- M + t(chol(U)) %*% Z %*% chol(V)
    return(X)
}

##

#' @export
head.resistance_parameters <- function(x, ...) {
    res <- lapply(x[-length(x)], head, ...)
    res[[names(x[length(x)])]] <- x[[length(x)]]
    class(res) <- "resistance_parameters"
    return(res)
}

#' @export
tail.resistance_parameters <- function(x, ...) {
    res <- lapply(x[-length(x)], tail, ...)
    res[[names(x[length(x)])]] <- x[[length(x)]]
    class(res) <- "resistance_parameters"
    return(res)
}

#' @export
as_tibble.resistance_parameters <- function(x, ...) {
    res <- cbind(Week = x$W, x$Parameters) %>% as_tibble() %>% 
        gather(ParameterName, Parameter, -Week)
    
    return(res)
}

#' @export
as_tibble.post_w <- function(x, ...) {
    return(tibble("Week" = x$W, "Posterior" = x$Posterior))
}

#' @export
mean.post_w <- function(x, ...) {
    W = x$W
    post_prob = x$Posterior
    
    res <- sum(W * post_prob)
    return(res)
}

#' @export
cumsum.post_w <- function(x) {
    W = x$W
    post_prob = x$Posterior
    
    res <- cumsum(post_prob)
    names(res) <- W
    return(res)
}

##

#' @title Sample from histogram
#' 
#' @description Sample from the bins of a histogram.
#' 
#' @param x An object of class '\link{histogram}'.
#' @param size The number of samples to draw.
#' 
#' @return Numeric vector.
#' @export
sample_histogram <- function(x, size = 1) {
    bins <- sample(x = length(x$mids), size, p = x$density, replace = TRUE)
    res <- runif(length(bins), x$breaks[bins], x$breaks[bins + 1]) 
    return(res)
}

#' @export
ggplot.histogram <- function(xhist, freq = FALSE, density = TRUE, xlab = "", ylab = NULL) {
    cc <- xhist$density
    if (freq || !density) {
        cc <- xhist$counts
    }
    
    bin_widths <- mean(diff(xhist$mids))
    ggplot(data.frame(x = xhist$mids, y = cc), aes(x = x, y = y)) + 
        geom_bar(stat = "identity", width = bin_widths - 0.002, colour = "black", fill = "white") + 
        #scale_x_continuous(breaks = seq(0, 1, 0.2)) + 
        xlab(xlab) + ylab(ifelse(is.null(ylab), ifelse(freq, "Count", "Density"), ylab)) + 
        theme_bw()
}

##

#' @title Compare windows
#' 
#' @description Compare all windows, with size W, of two vectors, 'I1' and 'I2', with the added requirement that starting index of the windows being compared of two other vectors 'SOC1' and 'SOC2' are within some tolerance 'delta'. The windows are compared by counting the number of indices where they differ more than some tolerance epsilon.
#' 
#' @param I1 Numeric vector of length 'T1'.
#' @param I2 Numeric vector of length 'T2'.
#' @param V1 Numeric vector of length 'T1'.
#' @param V2 Numeric vector of length 'T2'.
#' @param W Scalar: the window size.
#' @param R Scalar: the maximum allowed number of place the two windows exceed epsilon.
#' @param epsilon Tolerance used when comparing windows of 'I1' and 'I2'.
#' @param delta Tolerance used when comparing the initial index of the window in 'SOC1' and 'SOC2'.
#' @param trace TRUE/FLASE: Show trace?
#' @param trace_limit Scalar: Limit the trace shown to 'trace_limit'.
#' @param return_tibble TRUE/FALSE: Should a tibble be returned?
#' 
#' @details No default values are used for 'W', 'epsilon', or 'delta' as the choise of e.g. 'epsilon' and 'delta' are heavily dependent on the scale of the 'Is' and the 'SOCs'. Furthermore, the 'SOCs' are needed to reduce the amount of "relevant" information, otherwise the method will run out of memory quickly.
#'
#' @return Either a list of lists, or a \link{tibble}.
#' @export
compare_windows <- function(I1, I2, V1, V2, W, R, epsilon, delta, trace = TRUE, trace_limit = 10000, return_tibble = TRUE) {
    if (length(I1) != length(V1)) {
        stop("'I1' and 'V1' must have the same length.")
    }
    
    if (length(I2) != length(V2)) {
        stop("'I2' and 'V2' must have the same length.")
    }
    
    if ((W > length(I1)) || (W > length(I2))) {
        stop("The window size, 'W', should be smaller both 'I1'/'V1' and 'I2'/'V2'.")
    }
    
    if (R > W) {
        stop("'R' should be smaller than, or equal to, 'W'.́")
    }
    
    if (length(epsilon) == 1) {
        epsilon <- rep(epsilon, 2)
    }
    else if (length(epsilon) > 2) {
        epsilon <- epsilon[1:2]
    }
    
    system.time({
        res <- BatteryPrognostics:::compare_windows_cpp(I1, I2, V1, V2, W, R, epsilon, delta, trace, trace_limit)
    })
    
    if (return_tibble) {
        res <- tibble() 
    }
    
    return(res)
}

#' @title Compare windows using genetic algorithm
#' 
#' @description Compare windows of two vectors 'I1' and 'I2' using variable window sizes in a multi-objective GA.
#' 
#' @param I1 Numeric vector of length 'T1'.
#' @param I2 Numeric vector of length 'T2'.
#' @param V1 Numeric vector of length 'T1'.
#' @param V2 Numeric vector of length 'T2'.
#' @param Temp1 Numeric vector of length 'T1'.
#' @param Temp2 Numeric vector of length 'T2'.
#' @param MutationWindow Scalar: The size of the mutation window.
#' @param RestrictTemperature Numeric: The maximal allowed initial difference between the temperatures.
#' @param RestrictVoltage Numeric: The maximal allowed difference between the beginning and end between the voltages.
#' @param Wmin Scalar: Minimum the window size.
#' @param Wmax Scalar: Maximum the window size.
#' @param Imin Scalar: Minumum allowed starting current (A).
#' @param N_evolution Scalar: The number of evolutions of the GA.
#' @param N_keep Scalar: The number of individuals to keep form the GA.
#' @param trace TRUE/FLASE: Show trace?
#' @param trace_limit Scalar: Limit the trace shown to 'trace_limit'.
#' @param return_tibble TRUE/FALSE: Should a tibble be returned?
#' 
#' @details No default values are used for 'W', 'epsilon', or 'delta' as the choise of e.g. 'epsilon' and 'delta' are heavily dependent on the scale of the 'Is' and the 'SOCs'. Furthermore, the 'SOCs' are needed to reduce the amount of "relevant" information, otherwise the method will run out of memory quickly.
#'
#' @return Either a list of lists, or a \link{tibble}.
#' @export
compare_windows_ga <- function(I1, I2, V1, V2, Temp1, Temp2, 
                               MutationWindow, RestrictTemperature, RestrictVoltage,
                               Wmin, Wmax, Imin,
                               N_evolution = 2e6, N_keep = 100, 
                               trace = TRUE, trace_limit = 10000, return_tibble = TRUE) {
    if (length(I1) != length(V1)) {
        stop("'I1' and 'V1' must have the same length.")
    }
    
    if (length(I2) != length(V2)) {
        stop("'I2' and 'V2' must have the same length.")
    }
    
    if (length(I1) != length(Temp1)) {
        stop("'I1' and 'Temp1' must have the same length.")
    }
    
    if (length(I2) != length(Temp2)) {
        stop("'I2' and 'Temp2' must have the same length.")
    }
    
    if ((Wmin > length(I1)) || (Wmin > length(I2))) {
        stop("The minimum window size, 'Wmin', should be smaller both 'I1'/'V1' and 'I2'/'V2'.")
    }
    
    if ((Wmax > length(I1)) || (Wmax > length(I2))) {
        stop("The maximum window size, 'Wmax', should be smaller both 'I1'/'V1' and 'I2'/'V2'.")
    }
    
    if (Wmax < Wmin) {
        stop("The maximum window size, 'Wmax', should be larger than minimum window size, 'Wmin'.")
    }
    
    res <- BatteryPrognostics:::compare_windows_ga_cpp(I1, I2, V1, V2, Temp1, Temp2, 
                                                       MutationWindow, RestrictTemperature, RestrictVoltage, 
                                                       Wmin, Wmax, Imin,
                                                       N_evolution, N_keep, trace, trace_limit)
    
    if (return_tibble) {
        res <- tibble() 
    }
    
    return(res)
}

