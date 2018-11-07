##' @title Extract resistance
##' 
##' @description Extracts the resistance from a dynamic current/voltage profile.
##' 
##' @param I Vector containing the current.
##' @param V Vector containing the voltage.
##' @param T_s Vector containing the time step sizes (in seconds). If 'NULL' it is set as '1:length(I)'.
##' @param epsilon The lower limit used for identifying change in current. Default is '1e-3'.
##' @param Q_max The maximal capacity of the battery cell.
##' @param eta The Coulombic effeciency.
##' @param SOC_0 The initial value of the state-of-charge.
##' @param convert_to_tibble Organise in \link{tibble}.
##' 
##' @return A list containing the extracted resistance and their time scales extracted from the profile.
##' @export
extract_resistance <- function(I, V, T_s = NULL, epsilon = 1e-3, Q_max = 2.56, eta = 1.0, SOC_0 = 0.0, 
                               convert_to_tibble = FALSE) {
    if (is.null(I) | is.null(V)) {
        stop("Both 'I' and 'V' has to be provided.")
    }
    
    if (is.null(T_s)) {
        T_s = seq_along(I)
    }
    
    if (is.null(epsilon)) {
        epsilon = 1e-3
    }
    
    res_matrix <- extract_resistance_cpp(I, V, T_s, epsilon, Q_max, eta, SOC_0)
    class(res_matrix) <- "ext_res"
    
    if (convert_to_tibble) {
        res_tibble <- as_tibble(res_matrix)
        data_tibble <- tibble("I" = I, "V" = V, "T_s" = T_s) %>% 
            bind_cols(res_tibble) %>% 
            group_by(ID) %>% 
            mutate("SL" = max(S)) %>% 
            ungroup() %>% 
            mutate("Day" = floor(T_s / 3600 / 24))
        
        return(data_tibble)
    }
    else {
        return(res_matrix)
    }
}

##' @export
as_tibble.ext_res <- function(x, ...) {
    res_tibble <-x[1:2] %>% enframe() %>%
        unnest() %>%
        group_by(name) %>%
        mutate(TT = 1:n()) %>%
        spread(key = "name", value = "value") %>%
        select(-TT) %>%
        mutate(SOC = x$SOC,
               ID = x$ID) 
    
    return(res_tibble)
}

##' @title 2d-density plot of the extracted resistance 
##' 
##' @description A wrapping function calling ggplot with 'stat_density_2d' and 'ggMarginal' (from the 'ggExtra'-package). 
##' 
##' @param ext_res_tibble \link{tibble} made with the \link{extract_resistance}-function. 
##' @param restrict_dim_y Quantiles used to restrict the y-axis. See details.
##' @param include_marginals TRUE/FALSE: Should marginal plots be included? Note: requires 'ggExtra'.
##' @param dims The number of dimensions to be plotted.
##' @param facet_variable String passed to 'facet_wrap'.
##' @param plot_type Character of the type of plot. Takes the values 'density', 'hex', and 'bin'. 
##' @param n_bins The number of bins passed to the 'hex' and 'bin' plot types.
##' 
##' @details If 'restrict_dim_y' has length 1 it is taken as a restriction on the upper quantile. If it is 'NULL' the y-axis is not restricted.
##' 
##' @return ggplot-object.
##' @export
plot_extracted_resistance <- function(ext_res_tibble, restrict_dim_y = NULL, 
                                      include_marginals = FALSE, dims = c("SOC", "R"), 
                                      facet_variable = NULL, 
                                      plot_type = "hex", n_bins = NULL) {
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
    
    return(p)
}

