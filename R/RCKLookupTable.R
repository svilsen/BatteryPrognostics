#' @title Extract look-up tables.
#' 
#' @description Extracts RC-K parameter look-up tables from simple charge/discharge cycling profile.
#' 
#' @param current A charge/discharge cycling current profile. Note: assumes that negative current implies discharging of the battery.
#' @param voltage A charge/discharge cycling voltage profile.
#' @param time_s The sampling time.
#' @param t_k The time intervals for each of the RC-units.
#' @param h R0 reduction fraction. NOTE: THIS SHOULD NOT BE NECESSARY AND, THEREFORE, SET TO 1.
#' @param SOCValues A vector of possible 'SOC' values.
#' @param IValues A vector of possible 'Current' values.
#' 
#' @return A list of look-up table matrices.
#' @export
extractLookupTables <- function(current, voltage, time_s, t_k, 
                                SOCValues = seq(5, 95, 5), 
                                IValues = c(-10, -7.5, -5, -2.5, -1.25, -0.25, 0.25, 1.25, 2.5, 5, 7.5, 10), 
                                h = 1) {
    v_tau <- function(t, V, tau) {
        return(V * exp(-t / tau))
    }
    
    if (is.null(h) || (h < 0)) {
        h = 0
    }
    else if (h > 1) {
        h = 1
    } 
    
    RawTibble = tibble(V = c(voltage), I = c(current), Time_s = c(time_s)) %>% 
        mutate(CD = ifelse(I == 0, 0, ifelse(I < 0, -1, 1))) 
    
    rle_CD = rle(RawTibble$CD)
    RawTibble <- RawTibble %>% 
        mutate(RLECD = rep(1:length(rle_CD$lengths), times = rle_CD$lengths)) %>% 
        group_by(RLECD, CD) %>% 
        mutate(DeltaI = median(I), CRate = round(DeltaI / 2.5, 2)) 
    
    soc_stage <- c(1, cumsum((RawTibble$CRate[-length(RawTibble$CRate)] == -4.0) & (RawTibble$CRate[-1] != -4.0)) + 1)
    RawTibble$SOC = 5 * soc_stage
    
    R0 = matrix(0, nrow = length(SOCValues), ncol = length(IValues))
    
    K = length(t_k)
    Rk = rep(list(matrix(0, nrow = length(SOCValues), ncol = length(IValues))), times = K)
    Ck = rep(list(matrix(0, nrow = length(SOCValues), ncol = length(IValues))), times = K)
    
    OCV = matrix(0, nrow = length(SOCValues), ncol = length(IValues))
    for (i in seq_along(SOCValues)) {
        for (j in seq_along(IValues)) {
            # cat("i =", i, ":: j =", j, "\n")
            ## Set-up
            soc = SOCValues[i]
            cr = IValues[j] / 2.5
            
            ## OCV Calculation
            if (abs(cr + 4) < 1e-6) {
                which_ij = which((abs(RawTibble$SOC - soc) < 1e-6))
                OCV[i, j] = RawTibble$V[which_ij[length(which_ij)]]
            } else  {
                which_ij_ = which((abs(RawTibble$SOC - soc) < 1e-6) & (abs(RawTibble$CRate - cr) < 1e-6))
                
                if (cr == 0.5) {
                    which_ij_ = which_ij_[(which(diff(which_ij_) > 1) + 1):length(which_ij_)]
                }
                
                cr_mod = -cr
                if (cr < 0) {
                    cr_mod = abs(cr) + 1
                    if (cr > -1) {
                        cr_mod = ifelse(cr < -0.2, 1, 0.5)
                    } 
                }
                
                which_ij = which((abs(RawTibble$SOC - soc) < 1e-6) & (abs(RawTibble$CRate - cr_mod) < 1e-6))
                if ((cr_mod == 0.5)) {
                    which_ij = which_ij[(which(diff(which_ij) > 1) + 1):length(which_ij)]
                } 
                
                if (sign(cr) == -1) {
                    OCV[i, j] = max((RawTibble$V[(which_ij_[length(which_ij_)] + 1):(which_ij[1] - 1)]))
                }
                else {
                    OCV[i, j] = min((RawTibble$V[(which_ij_[length(which_ij_)] + 1):(which_ij[1] - 1)]))
                }
            }
            
            ## 
            which_ij = which((abs(RawTibble$SOC - soc) < 1e-6) & (abs(RawTibble$CRate - cr) < 1e-6))
            which_ij_last = which_ij[length(which_ij)]
            if ((cr == 0.5)) {
                which_ij = which_ij[(which(diff(which_ij) > 1) + 1):length(which_ij)]
            }
            
            ## R0
            V_H = RawTibble$V[which_ij[1]] * h + RawTibble$V[which_ij[1] - 1] * (1 - h)
            R0[i, j] = (V_H - RawTibble$V[which_ij[1] - 1]) /  IValues[j]
            
            ## Rk and Ck Lists
            t_1 = numeric(K)
            t_2 = numeric(K)
            tau = numeric(K)
            V = numeric(K)
            if (abs(cr + 4) > 1e-8) {
                for (k in 0:(K - 1)) {
                    t_1[K - k] = RawTibble$Time_s[which_ij_last + t_k[[K - k]][1]] - RawTibble$Time_s[which_ij_last]
                    t_2[K - k] = RawTibble$Time_s[which_ij_last + t_k[[K - k]][2]] - RawTibble$Time_s[which_ij_last]
                    
                    V_k1 = OCV[i, j] - RawTibble$V[which_ij_last + t_k[[K - k]][1]]
                    V_k2 = OCV[i, j] - RawTibble$V[which_ij_last + t_k[[K - k]][2]]
                    
                    if (k > 0) {
                        for (k_ in 0:(k - 1)) {
                            V_k1 = V_k1 - v_tau(t_1[K - k], V[K - k_], tau[K - k_])
                            V_k2 = V_k2 - v_tau(t_2[K - k], V[K - k_], tau[K - k_])
                        }
                    }
                    
                    if (sign(V_k1) != sign(V_k2)) 
                        V_k2 = sign(V_k1) * abs(V_k2)
                    
                    if (abs(V_k1) < abs(V_k2)) {
                        V_temp = V_k1
                        V_k2 = V_k1
                        V_k1 = V_temp
                    }
                    
                    #
                    tau[K - k] = (t_2[K - k] - t_1[K - k]) / log(V_k1 / V_k2)
                    if (is.infinite(tau[K - k]) | is.nan(tau[K - k])) {
                        tau[K - k] = 1L / .Machine$double.eps
                    } 
                    
                    V[K - k] = V_k1 * exp(t_1[K - k] / tau[K - k])
                    
                    # Rk[[K - k]][i, j] = abs((V_k2 - V_k1) / IValues[j])
                    Rk[[K - k]][i, j] = V[K - k] / (-IValues[j] * (1L - exp(-(RawTibble$Time_s[which_ij[length(which_ij)]] - RawTibble$Time_s[which_ij[1]]) / tau[K - k])))
                    if (Rk[[K - k]][i, j] < 0) 
                        Rk[[K - k]][i, j] = abs(Rk[[K - k]][i, j])
                    
                    Ck[[K - k]][i, j] = tau[K - k] / (capacityScalingParameters[[K - k]][i] * Rk[[K - k]][i, j])
                }
            }
        }
    }
    
    ## Temp. CRate = -4 fixes for the R_k, C_k, and OCV tables.
    Rk <- lapply(Rk, function(x) {x[, 1] = x[, ncol(x)]; return(x)})
    Ck <- lapply(Ck, function(x) {x[, 1] = x[, ncol(x)]; return(x)})
    OCV[, 1] = matrix(OCV[, 2] - 4 * (OCV[, 12] - OCV[, 11]), ncol = 1)
    OCV_ = OCV
    OCV = matrix(c(apply(OCV_[, 1:6], 1, mean), apply(OCV_[, 7:12], 1, mean)), ncol = 2, byrow = F)
    
    # CapLimited = c(2.55493, 2.55554, 2.56072, 2.57195, 2.57751, 2.57878, 2.57275, 2.58152, 2.57692, 2.57347, 2.56234, 2.55717)
    # IListLimited = c(-10.00, -7.50, -5.00, -2.50, -1.25, -0.25, 0.25, 1.25, 2.50, 5.00, 7.50, 10.00)
    # lmCap <- lm(CapLimited ~ IListLimited + I(IListLimited^2))
    
    # Cap = matrix(predict(lmCap, newdata = data.frame(IListLimited = IValues)), nrow = 1)
    # Cap = c(2.55493, 2.55554, 2.56072, 2.57195, 2.57751, 2.57878, 2.57275, 2.58152, 2.57692, 2.57347, 2.56234, 2.55717)
    res <- list(Tables = list(R0 = R0, Rk = Rk, Ck = Ck, OCV = OCV), #Cap = Cap), 
                Values = list(SOC = SOCValues, I = IValues))
    return(res)
}

#' @title Simulate RC-K model
#' 
#' @description Simulate battery voltage using look-up tables to approximate the parameters and assuming the Electrical Equivalent Circuit (EEC) model has 'K' RC-units.
#' 
#' @param current A current profile.
#' @param voltage A current profile reflecting the change in current (not the measured current).
#' @param extractedTables A list containing the possible 'SOC' and 'I' values, and the following look-up tables: 'R0', 'Rk', 'Ck', 'Cap', and 'OCV', where 'Rk' and 'Ck' are lists containing a look-up table for each RC-unit. 
#' @param dt The simulation step-size.
#' @param SOCStart The starting SOC value.
#' @param trace TRUE/FALSE: Show trace?
#' @param traceLimit Only show the trace for the iterations of specified 'traceLimit'.
#' 
#' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
#' @export
RCK <- function(current, extractedTables, dt = 1, SOCStart = 0.0, trace = T, traceLimit = 10000) { # voltage, 
    currentDifference = diff(current)
    currentChange = current
    currentFixed = current
    
    last_non_zero_current = current[1]
    last_non_zero_difference = current[1]
    for (k in seq(2, length(current))) {
        if (abs(currentDifference[k - 1]) > 1e-6) {
            last_non_zero_difference = currentDifference[k - 1]
            
            if (abs(current[k - 1]) > 1e-6) {
                last_non_zero_current = current[k - 1]
            }
        }
        
        if (abs(current[k]) < 1e-6) {
            currentChange[k] = last_non_zero_difference
            currentFixed[k] = last_non_zero_current
        }
    }
    
    res <- RCKCpp(current, currentChange, currentFixed, # voltage,
                  extractedTables$Tables$R0, extractedTables$Tables$Rk, extractedTables$Tables$Ck, 
                  matrix(extractedTables$Tables$Cap, nrow = 1), extractedTables$Tables$OCV,
                  extractedTables$Values$SOC / 100, extractedTables$Values$I,
                  dt, SOCStart, trace, traceLimit)
    return(res)
}