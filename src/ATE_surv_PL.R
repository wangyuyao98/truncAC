## PL estimator with IPW weights to handle confounding
# Estimate the survival probabilities at each time point in `t0.list`; use IPW weights to handle confounding

#' @param est_approach_PS The approach used to correct for the truncation bias in the propensity score estimation. Two approaches are implemented: "truncIPW.F" and "truncIPCW".
surv_PL_IPW <- function(dat, time.eval = NULL, X.name, Q.name, event.name, A.name,
                        cov.names.A, model.A = "logistic", est_approach_PS = "truncIPCW",
                        trim = 1e-7){
    
    names = c(X.name, Q.name, event.name, cov.names.A)
    
    if(est_approach_PS == "truncIPCW"){
        # estimate S_D using KM fit
        if(is.na(match("Y", names)) & is.na(match("delta.D", names))){
            dat$Y = dat[,X.name] - dat[,Q.name]    # residual censored event time
            dat$delta.D = 1- dat[,event.name]    # censoring indicator
        }else{
            stop("The names of the variables cannot be 'Y' or 'delta.D'. Conflict with intermediate variables in the function.")
        }
        
        KMfit_Sd = survfit(Surv(Y, delta.D) ~ 1, data = dat)
        KMcurve_Sd = stepfun(KMfit_Sd$time, c(1, KMfit_Sd$surv))
        SdY.vec = KMcurve_Sd(dat$Y)
        
        id.1 = (dat$delta==1)  # index for uncensored subjects
        dat.1 = dat[id.1, ]
        IPCW_weights.1 = 1/pmax(SdY.vec[id.1], trim)
        
        # Estimate G using IPCW weights and uncensored data
        if(sum(names == "Q2")){
            stop("The names of the variables cannot be 'Q2'. It is used in the middle of the computation.")
        }
        if(sum(names == "T2")){
            stop("The names of the variables cannot be 'T2'. It is used in the middle of the computation.")
        }
        dat.1$Q2 = tau - dat.1[,Q.name]
        dat.1$T2 = tau - dat.1[,X.name]
        dat.1$ones = 1
        
        tau = max(dat[,X.name],+1)
        PLfit_G = survfit(Surv(T2,Q2,ones)~1, data = dat.1, weights = IPCW_weights.1)
        PLcurve_G = stepfun(PLfit_G$time, c(1,PLfit_G$surv))
        Gt.vec.1 = PLcurve_G(dat.1$T2)
        trunc_weights.1 = 1/pmax(Gt.vec.1, trim)
        
        truncIPCW_weights.1 = pmax(IPCW_weights.1 * trunc_weights.1, trim)
        
        ## Fit the PS model
        PS_hat = est_PS(dat.1, dat, model = model.A, A.name, cov.names.A, weights = IPCW_weights.1*trunc_weights.1)
        
    }else{
        stop("Currently only the oprtion est_approach_PS = 'truncIPCW' is  implemented.")
    }
    
    
    ## IPTW to account for confounding + PL for estimating survival curve under random left truncation and right censoring
    IPTW_weights_a1 = 1/pmax(PS_hat, trim)
    IPTW_weights_a0 = 1/pmax(1-PS_hat, trim)
    id.a1 = (dat[,A.name] == 1)
    id.a0 = (dat[,A.name] == 0)
    
    # PL fit for a = 1
    formula.T = as.formula(paste("Surv(", Q.name,",",X.name,",", event.name, ") ~ 1", sep = ""))
    PLfit_a1 = survfit(formula.T, data = dat[id.a1,], weights = IPTW_weights_a1[id.a1])
    PLfit_a0 = survfit(formula.T, data = dat[id.a0,], weights = IPTW_weights_a0[id.a0])
    PLcurve_a1 = stepfun(PLfit_a1$time, c(1, PLfit_a1$surv))
    PLcurve_a0 = stepfun(PLfit_a0$time, c(1, PLfit_a0$surv))
    # plot(PLfit_a1$time, PLfit_a1$surv, type = "l")
    # lines(PLfit_a0$time, PLfit_a0$surv, type = "l")
    
    
   ## Compute the survival probabilities at time.eval
    if(is.null(time.eval)){
        time.eval = sort(unique(dat[,X.name]))
    }
    
    surv_a1 = PLcurve_a1(time.eval)
    surv_a0 = PLcurve_a0(time.eval)
    
    
    
    return(list(surv_a1 = surv_a1,
                surv_a0 = surv_a0,
                time.eval = time.eval))
    
}




surv_PL_IPW.boot <- function(dat, id, time.eval = NULL, X.name, Q.name, event.name, A.name,
                        cov.names.A, model.A = "logistic", est_approach_PS = "truncIPCW",
                        trim = 1e-7){
    
    dat = dat[id,]
    
    result = surv_PL_IPW(dat, time.eval, X.name, Q.name, event.name, A.name,
                         cov.names.A, model.A, est_approach_PS,
                         trim)
    
    surv_a1a0 = c(result$surv_a1, result$surv_a0)
    names(surv_a1a0) <- rep(result$time.eval, 2)
    
    return(surv_a1a0)
}



