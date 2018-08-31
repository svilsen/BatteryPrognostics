#' @title Extract look-up tables.
#' 
#' @description Extracts RC2 parameter look-up tables from simple charge/discharge cycling profile.
#' 
#' @param current A periodic current profile.
#' @param voltage A periodic voltage profile.
#' @param time_s The sampling time.
#' @param t_k The time intervals for each of the RC-units.
#' @param h R0 reduction fraction.
#' @param SOCList A vector of possible 'SOC' values.
#' @param IList A vector of possible 'Current' values.
#' 
#' @return A list of look-up table matrices.
#' @export
extractLookupTables <- function(current, voltage, time_s, t_k, 
                                capacityScalingParameters = NULL,
                                h = 0.65, SOCList = seq(5, 95, 5), 
                                IList = c(-10, -7.5, -5, -2.5, -1.25, -0.25, 0.25, 1.25, 2.5, 5, 7.5, 10)) {
    RawTibble = tibble(V = c(voltage), I = c(current), Time_s = c(time_s)) %>% 
        mutate(CD = ifelse(I == 0, 0, ifelse(I < 0, -1, 1))) 
    
    rle_CD = rle(RawTibble$CD)
    RawTibble <- RawTibble %>% 
        mutate(RLECD = rep(1:length(rle_CD$lengths), times = rle_CD$lengths)) %>% 
        group_by(RLECD, CD) %>% 
        mutate(DeltaI = median(I), CRate = round(DeltaI / 2.5, 2)) 
    
    soc_stage <- c(1, cumsum((RawTibble$CRate[-length(RawTibble$CRate)] == -4.0) & (RawTibble$CRate[-1] != -4.0)) + 1)
    RawTibble$SOC = 5 * soc_stage
    
    R0 = matrix(0, nrow = length(SOCList), ncol = length(IList))
    
    K = length(t_k)
    if (is.null(capacityScalingParameters)) {
        capacityScalingParameters = list(A = 1, K = rep(list(rep(1, length(IList))), times = K))
    }
    
    Rk = rep(list(matrix(0, nrow = length(SOCList), ncol = length(IList))), times = K)
    Ck = rep(list(matrix(0, nrow = length(SOCList), ncol = length(IList))), times = K)
    
    OCV = matrix(0, nrow = length(SOCList), ncol = length(IList))
    for (i in seq_along(SOCList)) {
        for (j in seq_along(IList)) {
            ## Set-up
            soc = SOCList[i]
            cr = IList[j] / 2.5
            
            which_ij = which((abs(RawTibble$SOC - soc) < 1e-6) & (abs(RawTibble$CRate - cr) < 1e-6))
            if ((cr == 0.5)) {
                which_ij = which_ij[(which(diff(which_ij) > 1) + 1):length(which_ij)]
            }
            
            ## R0
            V_H = RawTibble$V[which_ij[1]] * h + RawTibble$V[which_ij[1] - 1] * (1 - h)
            R0[i, j] = (V_H - RawTibble$V[which_ij[1] - 1]) /  IList[j] #RawTibble$DeltaI[which_ij[1]]
            
            ## Rk and Ck Lists
            for (k in 1:K) {
                t_k1 = RawTibble$Time_s[which_ij[t_k[[k]][1]]] - RawTibble$Time_s[which_ij[1]]
                t_k2 = RawTibble$Time_s[which_ij[t_k[[k]][2]]] - RawTibble$Time_s[which_ij[1]]
                
                V_k1 = RawTibble$V[which_ij[t_k[[k]][1]]]
                V_k2 = RawTibble$V[which_ij[t_k[[k]][2]]]
                
                tau_k = abs((t_k2 - t_k1) / log(V_k2 / V_k1))
                Rk[[k]][i, j] = (V_k2 - V_k1) / IList[j]
                
                # U_k = V_k1 * exp(t_k1 / tau_k)
                # Rk[[k]][i, j] = U_k / (IList[j] * (1 - exp(-(RawTibble$Time_s[which_ij[length(which_ij)]] - RawTibble$Time_s[which_ij[1]]) / tau_k)))
                Ck[[k]][i, j] = tau_k / (exp(log(Rk[[k]][i, j]) + log(capacityScalingParameters$A) + log(capacityScalingParameters$K[[k]][j])))
            }
            
            ## OCV Calculation
            if (abs(cr + 4) < 1e-6) {
                which_ij = which((abs(RawTibble$SOC - soc) < 1e-6))
                OCV[i, j] = RawTibble$V[which_ij[length(which_ij)]]
            } else  {
                cr_mod = -cr
                if (cr < 0) {
                    cr_mod = abs(cr) + 1
                    if (cr > -1) {
                        cr_mod = ifelse(cr < -0.2, 1, 0.5)
                    } 
                }
                
                which_ij = which((abs(RawTibble$SOC - soc) < 1e-6) & (abs(RawTibble$CRate - cr_mod) < 1e-6))
                if ((cr_mod == 0.5)) 
                    which_ij = which_ij[(which(diff(which_ij) > 1) + 1):length(which_ij)]
                
                OCV[i, j] = sum((RawTibble$V[(which_ij[1] - 4):(which_ij[1] - 1)] * seq(1, 4) / sum(seq(1, 4))))
            }
        }
    }
    
    CapLimited = c(2.55493, 2.55554, 2.56072, 2.57195, 2.57751, 2.57878, 2.57275, 2.58152, 2.57692, 2.57347, 2.56234, 2.55717)
    IListLimited = c(-10.00, -7.50, -5.00, -2.50, -1.25, -0.25, 0.25, 1.25, 2.50, 5.00, 7.50, 10.00)
    lmCap <- lm(CapLimited ~ IListLimited + I(IListLimited^2))
    
    Cap = matrix(predict(lmCap, newdata = data.frame(IListLimited = IList)), nrow = 1)
    res <- list(R0 = R0, Rk = Rk, Ck = Ck, Cap = Cap, OCV = OCV)
    return(res)
}

#' @title Simulate RC-K model
#' 
#' @description Simulate battery voltage using look-up tables to approximate the parameters and assuming the Electrical Equivalent Circuit (EEC) model has 'K' RC-units.
#' 
#' @param current A current profile.
#' @param voltage A current profile reflecting the change in current (not the measured current).
#' @param listOfTables A list containing the look-up tables: 'R0', 'Rk', 'Ck', 'Cap', and 'OCV', where 'Rk' and 'Ck' are lists containing a look-up table for each RC-unit. 
#' @param listOfDimensions A list containing vectors of the possible 'SOC' and 'Current' values in the look-up tables.
#' @param dt The simulation step-size.
#' @param SOCStart The starting SOC value.
#' @param trace TRUE/FALSE: Show trace?
#' @param traceLimit Only show the trace for the iterations of specified 'traceLimit'.
#' 
#' @return A list containing the voltage, the filtered voltage, the SOC, the polarisation voltage, and the state variance.
#' @export
RCK <- function(current, voltage, listOfTables, listOfDimensions, dt = 1, SOCStart = 0.0, trace = T, traceLimit = 10000) {
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
    
    res <- RCKCpp(current, currentChange, currentFixed, voltage, 
                   listOfTables$R0, listOfTables$Rk, listOfTables$Ck, listOfTables$Cap, listOfTables$OCV, 
                   listOfDimensions$SOCList, listOfDimensions$IList, 
                   dt, SOCStart, trace, traceLimit)
    return(res)
}