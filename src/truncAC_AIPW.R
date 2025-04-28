### trunc AIPW, truncC AIPW, truncAC AIPW


### trunc AIPW ----------------------------------------------------------------------------

# function that compute different parts involved in the estimating function (EF) for truncAIPW

# Fuz.mx: the matrix representing the CDF of T|(A,Z) for each subject in dat
# Gvz.mx: the matrix representing the CDF of Q|(A,Z) for each subject in dat

truncAIPW_transMean_EF <- function(dat, nu, Fuz.mx, Gvz.mx,
                                   T.name, Q.name, trim = 1e-7, u = NULL, v = NULL) {
    
    if (nrow(dat) != nrow(Fuz.mx) || nrow(dat) != nrow(Gvz.mx)) {
        stop("The number of rows of dat, Fuz.mx, and Gvz.mx must be the same.")
    }
    
    nn <- nrow(dat)
    Ttime <- as.numeric(dat[, T.name])
    Q <- as.numeric(dat[, Q.name])
    
    if(is.null(u)){
        u <- as.numeric(colnames(Fuz.mx))  # jumps.T
    }
    if(is.null(v)){
        v <- as.numeric(colnames(Gvz.mx))  # jumps.Q
    }
    
    nuu <- nu(u)            # used in matrix form
    nu_time <- nu(Ttime)     # pointwise evaluation
    
    result <- truncAIPW_transMean_EF_cpp(
        nn = nn,
        time = Ttime,
        Q = Q,
        Fuz_mx = Fuz.mx,
        Gvz_mx = Gvz.mx,
        u = u,
        v = v,
        nu_time = nu_time,
        nuu = nuu,
        trim = trim
    )
    
    result <- lapply(result, as.vector)
    
    return(result)
}


### truncC AIPW ----------------------------------------------------------------------------

# Return applying truncC_AIPW to the function \nu(T) and the constant function 1
#### ! Need testing - already tested?
#### Double check plus or minus augmentation term - should be correct, but be aware if the results are not as expected.
truncC_AIPW_transMean <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                                  X.name, Q.name, event.name, trim = 1e-7,
                                  u = NULL, v = NULL, d = NULL){
    
    ## compute the augmentation term
    nn <- nrow(dat)
    X <- as.numeric(dat[, X.name])
    Q <- as.numeric(dat[, Q.name])
    Delta <- as.integer(dat[, event.name])
    
    if(is.null(u)){
        u <- as.numeric(colnames(Fuz.mx))
    }
    if(is.null(v)){
        v <- as.numeric(colnames(Gvz.mx))
    }
    if(is.null(d)){
        d <- as.numeric(colnames(Sdz.mx))
    }
    
    # compute the IPCW weights
    X.res = X - Q
    Sdyz.vec <- diag(1 - CDF_eval_mx_cpp(X.res, 1 - Sdz.mx, d))
    id1 = (dat[, event.name] == 1)
    w_IPCW = rep(0, nrow(dat))
    w_IPCW[id1] = 1/pmax(Sdyz.vec[id1], trim)
    
    nuu <- nu(u)
    
    # apply truncAIPW
    efs = truncAIPW_transMean_EF(dat = dat, nu = nu, Fuz.mx = Fuz.mx, Gvz.mx = Gvz.mx, 
                                 T.name = X.name, Q.name = Q.name, trim = trim)
    Num_AIPW = as.numeric(efs$Num_AIPW)
    Den_AIPW = as.numeric(efs$Den_AIPW)
    
    Num_IPW.Q = as.numeric(efs$Num_IPW.Q)
    Den_IPW.Q = as.numeric(efs$Den_IPW.Q)
    
    result_Aug_QD = aug_QD(nn, nuu, X, Q, Delta, Fuz.mx, Gvz.mx, Sdz.mx,
                           u, v, d, Sdyz.vec, trim)
    
    Aug_QD_nu = result_Aug_QD$Aug_QD_nu
    Aug_QD_const1 = result_Aug_QD$Aug_QD_const1
    
    truncC_nu = w_IPCW * Num_AIPW + Aug_QD_nu 
    truncC_const1 = w_IPCW * Den_AIPW + Aug_QD_const1
    
    
    return(list(truncC_nu = truncC_nu, 
                truncC_const1 = truncC_const1))
}






### truncAC AIPW ----------------------------------------------------------------------------

## dr estimator for triple biases: confounding, left truncation, right censoring

#' @param mBootstrap whether multiplier boostrap procedure is used. If `mBootstrap = TRUE`, random weights are simulated and attached to the each observation, and the estimator is computed with this weighted sample. If `mBootstrap = FALSE`, the estimators are computed with the original data 
dr_truncAC_AIPW <- function(dat, nu, X.name, Q.name, event.name, A.name, 
                            model.T, model.Q, model.A, model.D, 
                            cov.names.T, cov.names.Q, cov.names.A, cov.names.D,
                            trim = 1e-7, trim.est = trim,
                            # mBootstrap = FALSE,
                            simulation = TRUE){
    
    names = c(X.name, Q.name, event.name, A.name, cov.names.T, cov.names.Q, cov.names.A, cov.names.D)
    if(sum(names == "delta.1")){
        stop("The names of the variables cannot be 'delta.1'. It is used in the middle of the computation.")
    }
    
    dat$delta.1 = rep(1, nrow(dat))
    event.name.Q = "delta.1"
    event.name.T = event.name
    
    dat.est = dat
    dat.fit = dat
    
    
    ### Estimate the nuisance parameters
    jumps.X = sort(dat[,X.name])
    jumps.Q = sort(dat[,Q.name])
    jumps.Y = sort(dat[,X.name] - dat[,Q.name])
    
    tau1 = min(jumps.X)
    tau2 = max(jumps.Q)
    tau3 = max(jumps.Y)
    u = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)
    v = c(tau1-1e-10, jumps.Q[jumps.Q>=tau1], max(jumps.Q)+1e-10)
    # v = c(tau1-1e-10, jumps.Q, max(jumps.Q)+1e-10)
    d = c(min(jumps.Y)-1e-10, jumps.Y, max(jumps.Y)+1e-10)
    
    
    
    # Estimate F and Sd
    Fuz.mx = est_F(dat.fit, dat.est, u, model.T, X.name, Q.name, event.name.T, cov.names.T, trim = trim.est)
    Sdz.mx = est_Sd(dat.fit, dat.est, d, model.D, X.name, Q.name, event.name.T, cov.names.D, trim = trim.est)   # d taking values on the residual time scale
    
    
    # Estimate G using IPCW weights
    id1.fit = (dat.fit[,event.name.T]==1)
    dat.fit.1 = dat.fit[id1.fit, ]  # uncensored data used to fit
    
    X.res = dat.fit[,X.name] - dat.fit[,Q.name]
    X.res.1 = X.res[id1.fit]
    Sd.vec.1 = diag(est_Sd(dat.fit, dat.fit.1, X.res.1, model.D, X.name, Q.name, event.name.T, cov.names.D, trim = trim.est) )
    w_IPCW.1 = 1/pmax(Sd.vec.1, trim)
    
    Gvz.mx = est_G(dat.fit.1, dat.est, v, model.Q, X.name, Q.name, event.name.Q, cov.names.Q,
                   weights = w_IPCW.1, trim = trim.est)
    
    
    # Estimate PS using IPCW and IPQW weights
    Gxz.vec.1 = pmax(diag(est_G(dat.fit.1, dat.fit.1, dat.fit.1$X, model.Q, X.name, Q.name, event.name.Q, cov.names.Q,
                                weights = w_IPCW.1, trim = trim.est)), trim)
    w_IPCW_IPQW.1 = 1/(pmax(Sd.vec.1, trim) * pmax(Gxz.vec.1, trim))
    # summary(w_IPCW_IPQW.1)
    PS = est_PS(dat.fit.1, dat.est, model.A, A.name, cov.names.A, weights = w_IPCW_IPQW.1, trim = trim.est)   # a vector of estimated PS for each individual
    
    # compute the mu_1 and mu_0 for the treatment augmentation
    u.T = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)  # need the values of F in all jumps of T to compute mu()
    
    dat_A1.est = dat.est
    dat_A1.est[,A.name] <- 1
    dat_A0.est = dat.est
    dat_A0.est[,A.name] <- 0
    if(simulation){
        dat_A1.est[,"AZ1"] = dat_A1.est[,A.name] * dat_A1.est[,"Z1"]    # need to change these two lines if the name of the covariates changes
        dat_A0.est[,"AZ1"] = dat_A0.est[,A.name] * dat_A0.est[,"Z1"]
    }
    
    Fuz_A1.mx = est_F(dat.fit, dat_A1.est, u.T, model.T, X.name, Q.name, event.name = event.name.T,
                      cov.names = cov.names.T, trim = trim.est)
    Fuz_A0.mx = est_F(dat.fit, dat_A0.est, u.T, model.T, X.name, Q.name, event.name = event.name.T,
                      cov.names = cov.names.T, trim = trim.est)
    
    nuuT.mx = matrix(rep(nu(u.T), nrow(dat.est)), nrow = nrow(dat.est), byrow = TRUE)
    tau.Tmax = max(u.T) + 1e-5
    mu_A1 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A1.mx))
    mu_A0 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A0.mx))
    
    
    ### Compute the estimators ----------------------------------------------------
    truncC_AIPW_result <- truncC_AIPW_transMean(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                                                X.name, Q.name, event.name, trim)
    
    truncC_AIPW_1 = truncC_AIPW_result$truncC_const1
    truncC_AIPW_nuT = truncC_AIPW_result$truncC_nu

    ## truncAC_AIPW
    A = dat[, A.name]
    PS_hat = PS
    mu1_hat = mu_A1
    mu0_hat = mu_A0
    Num.a1 = A/pmax(PS_hat, trim) * truncC_AIPW_nuT - (A-PS_hat)/pmax(PS_hat, trim) * mu1_hat * truncC_AIPW_1
    Num.a0 = (1-A)/pmax(1-PS_hat, trim) * truncC_AIPW_nuT + (A-PS_hat)/pmax(1-PS_hat, trim) * mu0_hat * truncC_AIPW_1
    Den = truncC_AIPW_1

    est_truncAC = c(est.a1 = mean(Num.a1)/mean(Den),
                    est.a0 = mean(Num.a0)/mean(Den))
    
    # The IF-based SE for the 3-bias truncAC_AIPW estimator
    n = nrow(dat)
    beta_hat = 1/mean(Den)
    IF.a1 = beta_hat * (Num.a1 - est_truncAC[1]*Den)
    IF.a0 = beta_hat * (Num.a0 - est_truncAC[2]*Den)
    SE_IF.a1 = sqrt(mean(IF.a1^2))/sqrt(n)
    SE_IF.a0 = sqrt(mean(IF.a0^2))/sqrt(n)
    SE_IF_ATE = sqrt(mean((IF.a1-IF.a0)^2))/sqrt(n)
    SE_IF_truncAC = c(SE.a1 = SE_IF.a1, SE.a0 = SE_IF.a0, SE_IF_ATE = SE_IF_ATE)
    
    
    
    ## IPW
    v <- as.numeric(colnames(Gvz.mx))
    d <- as.numeric(colnames(Sdz.mx))
    
    Delta = dat[, event.name]
    GXZ_hat = diag(CDF_eval_mx_cpp(dat.est[,X.name], Gvz.mx, v))
    SdXZ_hat = 1 - diag(CDF_eval_mx_cpp(dat.est[,X.name] - dat.est[,Q.name], 1 - Sdz.mx, d))
    
    Num_IPW.a1 = A/pmax(PS_hat, trim) * 1/pmax(GXZ_hat, trim) * Delta/pmax(SdXZ_hat, trim) * nu(dat[, X.name])
    Den_IPW.a1 = A/pmax(PS_hat, trim) * 1/pmax(GXZ_hat, trim) * Delta/pmax(SdXZ_hat, trim)
    
    Num_IPW.a0 = (1-A)/pmax(1-PS_hat, trim) * 1/pmax(GXZ_hat, trim) * Delta/pmax(SdXZ_hat, trim) * nu(dat[, X.name])
    Den_IPW.a0 = (1-A)/pmax(1-PS_hat, trim) * 1/pmax(GXZ_hat, trim) * Delta/pmax(SdXZ_hat, trim)
    
    est_IPW = c(est.a1 = mean(Num_IPW.a1)/mean(Den_IPW.a1),
                est.a0 = mean(Num_IPW.a0)/mean(Den_IPW.a0))
    
    # The robust sandwich SE estimator for IPW
    beta_hat_IPW.a1 = 1/mean(Den_IPW.a1)
    beta_hat_IPW.a0 = 1/mean(Den_IPW.a0)
    SF_IPW.a1 = beta_hat_IPW.a1 *(Num_IPW.a1 - est_IPW[1]*Den_IPW.a1)
    SF_IPW.a0 = beta_hat_IPW.a0 *(Num_IPW.a0 - est_IPW[1]*Den_IPW.a0)
    SE_IPW.a1 = sqrt(mean(SF_IPW.a1^2))/sqrt(n)
    SE_IPW.a0 = sqrt(mean(SF_IPW.a0^2))/sqrt(n)
    SE_IPW.ATE = sqrt(mean((SF_IPW.a1-SF_IPW.a0)^2))/sqrt(n)
    SE_IPW = c(SE.a1 = SE_IPW.a1, SE.a0 = SE_IPW.a0, SE_IPW.ATE = SE_IPW.ATE)
    
    
    return(list(est_truncAC_AIPW = est_truncAC, SE_truncAC_AIPW = SE_IF_truncAC,
                est_IPW = est_IPW, SE_IPW = SE_IPW,
                truncC_AIPW_1 = truncC_AIPW_1,
                truncC_AIPW_nuT = truncC_AIPW_nuT,
                mu1_hat = mu1_hat,
                mu0_hat = mu0_hat,
                PS_hat = PS_hat,
                GXZ_hat = GXZ_hat,
                SdXZ_hat = SdXZ_hat,
                beta_hat_truncAC_AIPW = beta_hat,
                beta_hat_IPW.a1 = beta_hat_IPW.a1, beta_hat_IPW.a0 = beta_hat_IPW.a0))
    
}



dr_truncAC_AIPW_boot <- function(dat, id, nu, X.name, Q.name, event.name, A.name, 
                                 model.T, model.Q, model.A, model.D, 
                                 cov.names.T, cov.names.Q, cov.names.A, cov.names.D,
                                 trim, trim.est = trim, simulation = TRUE){
    
    dat = dat[id,]
    
    est_truncAC_AIPW = dr_truncAC_AIPW(dat, nu, X.name, Q.name, event.name, A.name,
                                       model.T, model.Q, model.A, model.D, 
                                       cov.names.T, cov.names.Q, cov.names.A, cov.names.D,
                                       trim = trim, trim.est = trim.est, simulation = simulation)
    
    est = rbind(truncAC_AIPW = est_truncAC_AIPW$est_truncAC_AIPW,
                IPW = est_truncAC_AIPW$est_IPW)
    est.a1 = est[,1]
    est.a0 = est[,2]
    
    est.a1a0 = c(est.a1, est.a0)
    names(est.a1a0) = c(paste(rownames(est), ".a1", sep = ""),
                        paste(rownames(est), ".a0", sep = ""))
    
    
    SE = rbind(SE.truncAC_AIPW = est_truncAC_AIPW$SE_truncAC_AIPW,
               SE.IPW = est_truncAC_AIPW$SE_IPW)
    SE.a1 = SE[,1]
    SE.a0 = SE[,2]
    SE.ATE = SE[,3]
    
    SE.a1a0 = c(SE.a1, SE.a0, SE.ATE)
    names(SE.a1a0) = c(paste(rownames(SE), ".a1", sep = ""),
                       paste(rownames(SE), ".a0", sep = ""),
                       paste(rownames(SE), ".ATE", sep = ""))
    
    
    
    return(c(est.a1a0, SE.a1a0))
    
}


### cf estimator for triple biases: confounding, left truncation, right censoring

## the cross fitted truncAC_AIPW estimator
## Two options are implemented for the 
## Use another sample splitting for estimating the weights involved in nuisance estimation and the nuisance estimation
## the cross fitted truncAC_AIPW estimator

#' @param cov.names.binary.T vector of names for binary covariates if \texttt{model = "pCox"} is used. For other models, \texttt{c(cov.names, cov.names,binary)} is used as the covariates.
#' @param options.F a list of arguments in est_F() that are used for specific methods, including trim, formula.survPen, df, nfolds, s, alpha, lambda  
#' @param options.G a list of arguments in est_G() that are used for specific methods, including trim, formula.survPen, df, nfolds, s, alpha, lambda, trunc.weighting, where 'trunc.weighting' is a logical value indicating whether truncation weights are use to estimate G; if FALSE, G is estimated from uncensored subjects with IPCW weights 
#' @param options.PS a list of arguments in est_PS() that are used for specific methods, including trim, df. Currently PS is always estimated using truncating weighting.
#' @param options.Sd a list of arguments in est_Sd() that are used for specific methods, including trim, mtry, ntree, formula.survPen, df, nfolds, s, alpha, lambda  

cf_truncAC_AIPW <- function(dat, K, nu, X.name, Q.name, event.name, A.name, 
                            model.T, model.Q, model.A, model.D, 
                            cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim = 1e-7,
                            cov.names.binary.T = NULL, cov.names.binary.Q = NULL, cov.names.binary.A = NULL, cov.names.binary.D = NULL,
                            options.F = NULL, options.G = NULL, 
                            options.PS = NULL, options.Sd = NULL,
                            est_approach_G = "truncIPW.F", est_approach_PS = "truncIPW.F",
                            Fuz.mx = NULL, Fuz_A1.mx = NULL, Fuz_A0.mx = NULL,
                            Gvz.mx = NULL, Sdz.mx = NULL, PS = NULL,
                            simulation = FALSE){
    
    dat_A1 = dat
    dat_A1[,A.name] <- 1
    dat_A0 = dat
    dat_A0[,A.name] <- 0
    if(simulation){
        dat_A1[,"AZ1"] = dat_A1[,A.name] * dat_A1[,"Z1"]    # need to change these two lines if the name of the covariates changes
        dat_A0[,"AZ1"] = dat_A0[,A.name] * dat_A0[,"Z1"]
    }
    
    if(is.null(Gvz.mx)){
        names = c(X.name, Q.name, event.name, A.name, cov.names.T, cov.names.Q, cov.names.A, cov.names.D)
        if(sum(names == "delta.1")){
            stop("The names of the variables cannot be 'delta.1'. It is used in the middle of the computation.")
        }
        dat$delta.1 = rep(1, nrow(dat))
        event.name.Q = "delta.1"
    }
    
    event.name.T = event.name
    
    n = nrow(dat)
    folds <- cvFolds(n, K)
    
    jumps.X = sort(dat[,X.name])
    jumps.Q = sort(dat[,Q.name])
    jumps.Y = sort(dat[,X.name] - dat[,Q.name])
    
    tau1 = min(jumps.X)
    tau2 = max(jumps.Q)
    tau3 = max(jumps.Y)
    tau.max = max(jumps.X)
    
    # Compute \V{\nu} and \V{1} using out of fold estimate of the nuisance parameters
    truncC_AIPW_1 = rep(NA, n)  # \V{1}
    truncC_AIPW_nuT = rep(NA, n)  # \V{\nu(T)}
    
    # Store some intermediate results that are useful for constructing IPW and naive estimators
    mu0_hat = rep(NA, n)
    mu1_hat = rep(NA, n)
    GXZ_hat = rep(NA, n)
    SdXZ_hat = rep(NA, n)
    PS_hat = rep(NA, n)
    
    
    # ct = 0
    folds <- cvFolds(n, K)
    for(k in 1:K){
        id.est = folds$subsets[folds$which == k]
        id.fit = folds$subsets[folds$which != k]
        dat.est = dat[id.est, ]
        dat.fit = dat[id.fit, ]
        dat_A1.est = dat_A1[id.est, ]
        dat_A0.est = dat_A0[id.est, ]
        
        
        ### Estimate the nuisance parameters
        
        ## Estimate F and F_A1 F_A0
        
        if(is.null(Fuz.mx)){
            u = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)
            
            # ### Begin debug -- line-by-line -----------------------------------------------------
            # time.eval = u
            # model = model.T
            # time.name = X.name
            # Q.name = Q.name
            # event.name = event.name.T
            # cov.names = cov.names.T
            # cov.names.binary = cov.names.binary.T
            # trim = options.F$trim
            # mtry = options.F$mtry
            # ntree = options.F$ntree
            # formula.survPen = options.F$formula.survPen
            # nfolds = options.F$nfolds
            # s = options.F$s
            # alpha = options.F$alpha
            # lambda = options.F$lambda
            # df = options.F$df
            # Fuz.mx = NULL
            # Fuz_A1.mx = NULL
            # Fuz_A0.mx = NULL
            # Gvz.mx = NULL
            # Sdz.mx = NULL
            # PS = NULL
            # weights = NULL
            # ### End debug -----------------------------------------------------------------------
            Fuz.mx.est = est_F(dat.fit, dat.est, u, model.T, X.name, Q.name, event.name.T, 
                               cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                               trim = options.F$trim, 
                               mtry = options.F$mtry, ntree = options.F$ntree, 
                               formula.survPen = options.F$formula.survPen, 
                               nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                               lambda = options.F$lambda, df = options.F$df)
        }else{
            Fuz.mx.est = Fuz.mx[id.est, ]
            u = as.numeric(colnames(Fuz.mx))
        }
        
        if(is.null(Fuz_A1.mx)){
            u_A1 = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)
            Fuz_A1.mx.est = est_F(dat.fit, dat_A1.est, u_A1, model.T, X.name, Q.name, event.name.T, 
                                  cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                                  trim = options.F$trim, 
                                  mtry = options.F$mtry, ntree = options.F$ntree, 
                                  formula.survPen = options.F$formula.survPen, 
                                  nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                  lambda = options.F$lambda, df = options.F$df)
        }else{
            Fuz_A1.mx.est = Fuz_A1.mx[id.est, ]
            u_A1 = as.numeric(colnames(Fuz_A1.mx))
        }
        
        if(is.null(Fuz_A0.mx)){
            u_A0 = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)
            Fuz_A0.mx.est = est_F(dat.fit, dat_A0.est, u_A0, model.T, X.name, Q.name, event.name.T, 
                                  cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                                  trim = options.F$trim, 
                                  mtry = options.F$mtry, ntree = options.F$ntree, 
                                  formula.survPen = options.F$formula.survPen, 
                                  nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                  lambda = options.F$lambda, df = options.F$df)
        }else{
            Fuz_A0.mx.est = Fuz_A0.mx[id.est, ]
            u_A0 = as.numeric(colnames(Fuz_A0.mx))
        }
        
        
        
        ## Estimate S_D
        if(is.null(Sdz.mx)){
            d = c(min(jumps.Y)-1e-10, jumps.Y, max(jumps.Y)+1e-10)
            Sdz.mx.est = est_Sd(dat.fit, dat.est, d, model.D, X.name, Q.name, event.name.T, 
                                cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                trim = options.Sd$trim, 
                                mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                formula.survPen = options.Sd$formula.survPen, 
                                nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                lambda = options.Sd$lambda, df = options.Sd$df)   # d taking values on the residual time scale
            
        }else{
            Sdz.mx.est = Sdz.mx[id.est, ]
            d <- as.numeric(colnames(Sdz.mx))
        }
        
        
        
        ## Estimate G - using truncation weights 1/{1-F(Q|A,Z)} estimated using fit.si, and use fit.sj to estimate G, i\neq j \in \{1,2\} 
        if(is.null(Gvz.mx)){
            v = c(tau1-1e-10, jumps.Q[jumps.Q>=tau1], max(jumps.Q)+1e-10)
            # v = c(tau1-1e-10, jumps.Q, max(jumps.Q)+1e-10)
            
            if(est_approach_G == "truncIPW.F"){
                Fq.vec.fit = NULL
                
                # Compute the truncation weights 1/{1-F(Q|A,Z)} 
                if(is.null(Fuz.mx)){
                    Fq.vec.fit = diag(est_F(dat.fit, dat.fit, dat.fit[,Q.name], 
                                            model.T, X.name, Q.name, event.name.T, 
                                            cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                                            trim = options.F$trim, 
                                            mtry = options.F$mtry, ntree = options.F$ntree, 
                                            formula.survPen = options.F$formula.survPen, 
                                            nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                            lambda = options.F$lambda, df = options.F$df))
                }else{
                    Fq.vec.fit = diag(CDF_eval_mx_cpp(dat.fit[ ,Q.name], Fuz.mx[id.fit,], u))
                }
                
                
                w_truncF.fit = 1/pmax(1-Fq.vec.fit, trim)
                
                # Estimate G using fit.s2
                Gvz.mx.est = est_G(dat.fit, dat.est, v, model.Q, X.name, Q.name, event.name.Q, 
                                   cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                   weights = w_truncF.fit, trunc = FALSE, 
                                   tau = tau.max, trim = options.G$trim,
                                   formula.survPen = options.G$formula.survPen, 
                                   nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                   lambda = options.G$lambda, df = options.G$df)
                
                
                
            }else if(est_approach_G == "IPCW"){
                
                id.fit_1 = (dat.fit[, event.name.T] == 1)
                dat.fit_1 = dat.fit[id.fit_1, ]
                X.res_1 = dat.fit_1[,X.name] - dat.fit_1[,Q.name]
                
                # Estimate the censoring weights
                if(is.null(Sdz.mx)){
                    Sdy.vec.fit_1 = diag(est_Sd(dat.fit, dat.fit_1, X.res_1, 
                                                model.D, X.name, Q.name, event.name.T,
                                                cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                                trim = options.Sd$trim, 
                                                mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                                formula.survPen = options.Sd$formula.survPen, 
                                                nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                                lambda = options.Sd$lambda, df = options.Sd$df))
                }else{
                    Sdy.vec.fit_1 = 1- diag(CDF_eval_mx_cpp(X.res_1, 1-Sdz.mx[id.fit,][id.fit_1,], d))
                }
                
                w_IPCW_1 = 1/pmax(Sdy.vec.fit_1, trim)
                
                Gvz.mx.est = est_G(dat.fit_1, dat.est, v, model.Q, X.name, Q.name, event.name.Q, 
                                   cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                   weights = w_IPCW_1, tau = tau.max,
                                   trim = options.G$trim,
                                   formula.survPen = options.G$formula.survPen, 
                                   nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                   lambda = options.G$lambda, df = options.G$df)
                
                
            }else{
                stop("'est_approach_G' Should be either 'truncIPW.F' or 'IPCW'. ")
            }
            
            
        }else{
            Gvz.mx.est = Gvz.mx[id.est, ]
            v <- as.numeric(colnames(Gvz.mx))
        }
        
        
        
        ## Estimate PS
        if(is.null(PS)){
            
            if(est_approach_PS == "truncIPW.F"){
                
                if(is.null(Fq.vec.fit)){ # The truncation weights from F haven't been computed yet
                    
                    # Compute the truncation weights 1/{1-F(Q|A,Z)} estimated using fit.s1
                    if(is.null(Fuz.mx)){
                        Fq.vec.fit = diag(est_F(dat.fit, dat.fit, dat.fit[,Q.name], 
                                                model.T, X.name, Q.name, event.name.T, 
                                                cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                                                trim = options.F$trim, 
                                                mtry = options.F$mtry, ntree = options.F$ntree, 
                                                formula.survPen = options.F$formula.survPen, 
                                                nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                                lambda = options.F$lambda, df = options.F$df))
                    }else{
                        Fq.vec.fit = diag(CDF_eval_mx_cpp(dat.fit[ ,Q.name], Fuz.mx[id.fit,], u))
                    }
                    
                    w_truncF.fit = 1/pmax(1-Fq.vec.fit, trim)
                    
                } # End if(is.null(Fq.vec.fit))
                
                # Estimate PS using truncation weights 1/(1-F(Q|A,Z))
                PS.est = est_PS(dat.fit, dat.est, model.A, A.name, 
                                cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                weights = w_truncF.fit,
                                trim = options.PS$trim, df = options.PS$df,
                                ntree = options.PS$ntree)   # a vector of estimated PS for each individual
                
                
            }else if(est_approach_PS == "truncIPCW"){
                
                id.fit_1 = (dat.fit[, event.name.T] == 1)
                dat.fit_1 = dat.fit[id.fit_1, ]
                X.res_1 = dat.fit_1[,X.name] - dat.fit_1[,Q.name]
                X_1 = dat.fit_1[,X.name] 
                
                # Compute the IPCW weights
                if(is.null(Sdz.mx)){
                    Sdy.vec.fit = diag(est_Sd(dat.fit, dat.fit_1, X.res_1, 
                                              model.D, X.name, Q.name, event.name.T,
                                              cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                              trim = options.Sd$trim, 
                                              mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                              formula.survPen = options.Sd$formula.survPen, 
                                              nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                              lambda = options.Sd$lambda, df = options.Sd$df))
                }else{
                    Sdy.vec.fit = 1- diag(CDF_eval_mx_cpp(X.res_1, 1-Sdz.mx[id.fit,][id.fit_1,], d))
                }
                
                # IPCW weights
                w_IPCW.fit_1 = 1/pmax(Sdy.vec.fit, trim)   # for the uncensored subjects
                
                if(is.null(Gvz.mx)){
                    # For uncensored subjects
                    Gx.vec.fit_1 = diag( est_G(dat.fit_1, dat.fit_1, X_1, model.Q, X.name, Q.name, event.name.Q, 
                                               cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                               weights = w_IPCW.fit_1, tau = tau.max,
                                               trim = options.G$trim,
                                               formula.survPen = options.G$formula.survPen, 
                                               nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                               lambda = options.G$lambda, df = options.G$df) )
                    
                }else{
                    Gx.vec.fit_1 = diag(CDF_eval_mx_cpp(X_1, Gvz.mx[id.fit,][id.fit_1,], v)) # For uncensored subjects
                }
                
                w_trunc.fit_1 = 1/pmax(Gx.vec.fit_1, trim)
                
                
                # Estimate PS using uncensored subjects with weights 1/{G(X|A,Z)S_D(X-Q|A,Z)}
                PS.est = est_PS(dat.fit_1, dat.est, model.A, A.name, 
                                cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                weights = w_trunc.fit_1 * w_IPCW.fit_1,
                                trim = options.PS$trim, df = options.PS$df,
                                ntree = options.PS$ntree) 
                
            }else{
                stop("'est_approach_PS' Should be either 'truncIPW.F' or 'truncIPCW'. ")
            }
            
            
        }else{
            PS.est = PS[id.est]
        }
        
        
        # Compute the mu_1(Z), mu_0(Z) for the treatment augmentation, and m(Z) invoved in the R-learner
        nuuT_A1.mx = matrix(rep(nu(u_A1), nrow(dat.est)), nrow = nrow(dat.est), byrow = TRUE)
        nuuT_A0.mx = matrix(rep(nu(u_A0), nrow(dat.est)), nrow = nrow(dat.est), byrow = TRUE)
        tau.Tmax_A1 = max(u_A1) + 1e-5
        tau.Tmax_A0 = max(u_A0) + 1e-5
        mu_A1.est = as.vector(int_fmx_dF_cpp(tau.Tmax_A1, nuuT_A1.mx, Fuz_A1.mx.est, u_A1))
        mu_A0.est = as.vector(int_fmx_dF_cpp(tau.Tmax_A0, nuuT_A0.mx, Fuz_A0.mx.est, u_A0))
        
        
        
        # # Debug line-by-line ------------------------------------------------------
        # dat = dat.est
        # nu = nu 
        # Fuz.mx = Fuz.mx.est
        # Gvz.mx = Gvz.mx.est
        # Sdz.mx = Sdz.mx.est
        # X.name = X.name
        # Q.name = Q.name
        # event.name = event.name
        # trim = trim
        # # End debug ---------------------------------------------------------------
        # Compute \V{1}, \V{\nu}
        truncC_AIPW_result <- truncC_AIPW_transMean(dat.est, nu, Fuz.mx.est, Gvz.mx.est, Sdz.mx.est,
                                                    X.name, Q.name, event.name, trim)
        
        # truncCAIPW_const1 = truncC_AIPW_result$truncC_const1
        truncC_AIPW_1[id.est] = truncC_AIPW_result$truncC_const1
        truncC_AIPW_nuT[id.est] = truncC_AIPW_result$truncC_nu
        
        
        ### Record some intermediate quantities
        mu1_hat[id.est] = mu_A1.est
        mu0_hat[id.est] = mu_A0.est
        PS_hat[id.est] = PS.est
        GXZ_hat[id.est] = diag(CDF_eval_mx_cpp(dat.est[,X.name], Gvz.mx.est, v))
        SdXZ_hat[id.est] = 1 - diag(CDF_eval_mx_cpp(dat.est[,X.name] - dat.est[,Q.name], 1-Sdz.mx.est, d))
        
    }
    
    # ### !!! Debug -- save the temp results -----------------------------------------
    # cf_truncC_AIPW_result = list(truncC_AIPW_1 = truncC_AIPW_1, truncC_AIPW_nuT = truncC_AIPW_nuT)
    # cf_truncC_details = list(dat = dat, nu = nu, # Fuz.mx = Fuz.mx, Gvz.mx = Gvz.mx, Sdz.mx = Sdz.mx, 
    #                          PS = PS_hat, mu_A1 = mu1_hat, mu_A0 = mu0_hat) 
    # save(cf_truncC_AIPW_result, cf_truncC_details, 
    #      file = "debug_results/cf_truncC_AIPW_results.rda")
    # ### End Debug ------------------------------------------------------------------
    # 
    
    
    A = dat[,A.name]
    Delta = dat[, event.name]
    
    ### Compute the estimator -----------
    # truncAC_AIPW
    Num.a1 = A/pmax(PS_hat, trim) * truncC_AIPW_nuT - (A-PS_hat)/pmax(PS_hat, trim) * mu1_hat * truncC_AIPW_1
    Num.a0 = (1-A)/pmax(1-PS_hat, trim) * truncC_AIPW_nuT + (A-PS_hat)/pmax(1-PS_hat, trim) * mu0_hat * truncC_AIPW_1
    Den = truncC_AIPW_1
    
    est_truncAC = c(est.a1 = mean(Num.a1)/mean(Den),
                    est.a0 = mean(Num.a0)/mean(Den))
    
    ## Add the IF-based SE for the 3-bias truncAC_AIPW estimator
    n = nrow(dat)
    beta_hat = 1/mean(Den)
    IF.a1 = beta_hat * (Num.a1 - est_truncAC[1]*Den)
    IF.a0 = beta_hat * (Num.a0 - est_truncAC[2]*Den)
    SE_IF.a1 = sqrt(mean(IF.a1^2))/sqrt(n)
    SE_IF.a0 = sqrt(mean(IF.a0^2))/sqrt(n)
    SE_IF_ATE = sqrt(mean((IF.a1-IF.a0)^2))/sqrt(n)
    SE_IF_truncAC = c(SE.a1 = SE_IF.a1, SE.a0 = SE_IF.a0, SE_IF_ATE = SE_IF_ATE)
    
    
    # IPW
    Num_IPW.a1 = A/pmax(PS_hat, trim) * 1/pmax(GXZ_hat, trim) * Delta/pmax(SdXZ_hat, trim) * nu(dat[, X.name])
    Den_IPW.a1 = A/pmax(PS_hat, trim) * 1/pmax(GXZ_hat, trim) * Delta/pmax(SdXZ_hat, trim)
    
    Num_IPW.a0 = (1-A)/pmax(1-PS_hat, trim) * 1/pmax(GXZ_hat, trim) * Delta/pmax(SdXZ_hat, trim) * nu(dat[, X.name])
    Den_IPW.a0 = (1-A)/pmax(1-PS_hat, trim) * 1/pmax(GXZ_hat, trim) * Delta/pmax(SdXZ_hat, trim)
    
    est_IPW = c(est.a1 = mean(Num_IPW.a1)/mean(Den_IPW.a1),
                est.a0 = mean(Num_IPW.a0)/mean(Den_IPW.a0))
    
    # Add the robust sandwich SE estimator for IPW
    beta_hat_IPW.a1 = 1/mean(Den_IPW.a1)
    beta_hat_IPW.a0 = 1/mean(Den_IPW.a0)
    SF_IPW.a1 = beta_hat_IPW.a1 *(Num_IPW.a1 - est_IPW[1]*Den_IPW.a1)
    SF_IPW.a0 = beta_hat_IPW.a0 *(Num_IPW.a0 - est_IPW[1]*Den_IPW.a0)
    SE_IPW.a1 = sqrt(mean(SF_IPW.a1^2))/sqrt(n)
    SE_IPW.a0 = sqrt(mean(SF_IPW.a0^2))/sqrt(n)
    SE_IPW.ATE = sqrt(mean((SF_IPW.a1-SF_IPW.a0)^2))/sqrt(n)
    SE_IPW = c(SE.a1 = SE_IPW.a1, SE.a0 = SE_IPW.a0, SE_IPW.ATE = SE_IPW.ATE)
    
    
    return(list(est_truncAC_AIPW = est_truncAC, SE_truncAC_AIPW = SE_IF_truncAC,
                est_IPW = est_IPW, SE_IPW = SE_IPW,
                truncC_AIPW_1 = truncC_AIPW_1,
                truncC_AIPW_nuT = truncC_AIPW_nuT,
                mu1_hat = mu1_hat,
                mu0_hat = mu0_hat,
                PS_hat = PS_hat,
                GXZ_hat = GXZ_hat,
                SdXZ_hat = SdXZ_hat,
                beta_hat_truncAC_AIPW = beta_hat,
                beta_hat_IPW.a1 = beta_hat_IPW.a1, beta_hat_IPW.a0 = beta_hat_IPW.a0))
    
}





cf_truncAC_AIPW_boot <- function(dat, id, K, nu, X.name, Q.name, event.name, A.name, 
                                 model.T, model.Q, model.A, model.D, 
                                 cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim = 1e-7,
                                 cov.names.binary.T = NULL, cov.names.binary.Q = NULL, cov.names.binary.A = NULL, cov.names.binary.D = NULL,
                                 options.F = NULL, options.G = NULL, 
                                 options.PS = NULL, options.Sd = NULL,
                                 est_approach_G = "truncIPW.F", est_approach_PS = "truncIPW.F",
                                 Fuz.mx = NULL, Fuz_A1.mx = NULL, Fuz_A0.mx = NULL,
                                 Gvz.mx = NULL, Sdz.mx = NULL, PS = NULL,
                                 simulation = FALSE){
    
    dat = dat[id, ]
    
    est_truncAC_AIPW = cf_truncAC_AIPW(dat, K, nu, X.name, Q.name, event.name, A.name, 
                                       model.T, model.Q, model.A, model.D, 
                                       cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim,
                                       cov.names.binary.T = cov.names.binary.T, cov.names.binary.Q = cov.names.binary.Q, cov.names.binary.A = cov.names.binary.A, cov.names.binary.D = cov.names.binary.D,
                                       options.F = options.F, options.G = options.G, 
                                       options.PS = options.PS, options.Sd = options.Sd,
                                       est_approach_G = est_approach_G, est_approach_PS = est_approach_PS,
                                       Fuz.mx = Fuz.mx, Fuz_A1.mx = Fuz_A1.mx, Fuz_A0.mx = Fuz_A0.mx,
                                       Gvz.mx = Gvz.mx , Sdz.mx = Sdz.mx, PS = PS,
                                       simulation = simulation)
    
    
    est = rbind(truncAC_AIPW = est_truncAC_AIPW$est_truncAC_AIPW,
                IPW = est_truncAC_AIPW$est_IPW)
    est.a1 = est[,1]
    est.a0 = est[,2]
    
    est.a1a0 = c(est.a1, est.a0)
    names(est.a1a0) = c(paste(rownames(est), ".a1", sep = ""),
                        paste(rownames(est), ".a0", sep = ""))
    
    
    
    SE = rbind(SE.truncAC_AIPW = est_truncAC_AIPW$SE_truncAC_AIPW,
               SE.IPW = est_truncAC_AIPW$SE_IPW)
    SE.a1 = SE[,1]
    SE.a0 = SE[,2]
    SE.ATE = SE[,3]
    
    SE.a1a0 = c(SE.a1, SE.a0, SE.ATE)
    names(SE.a1a0) = c(paste(rownames(SE), ".a1", sep = ""),
                       paste(rownames(SE), ".a0", sep = ""),
                       paste(rownames(SE), ".ATE", sep = ""))
    
    
    
    return(c(est.a1a0, SE.a1a0))
}




