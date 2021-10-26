#' @title Rainflow 
#' 
#' @description An implementation of the rainflow cycle counting algorithm by Adam Nieslony. 
#' 
#' @param x Vector: input data. 
#' @param flm Numeric: the fixed load mean (default: 0).
#' @param ul Numeric: the ultimate load (default: 1e16).
#' @param pls Numeric: the partial load scaling (default: 0.5).
#' 
#' @return A tibble of the cycle counts
#' @export
rainflow <- function(x, flm = 0, ul = 1e16, pls = 0.5) {
    
    res <- rainflow_cpp(x, flm, ul, pls)
    res <- res %>% 
        as_tibble() %>% 
        rename("" = "V1", "" = "V2", "" = "V3", "" = "V4", "" = "V5",)
    
    return(res)
}
