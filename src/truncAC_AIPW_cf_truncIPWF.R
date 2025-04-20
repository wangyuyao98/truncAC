## the cross fitted truncAC_AIPW estimator

#' @param cov.names.binary.T vector of names for binary covariates if \texttt{model = "pCox"} is used. For other models, \texttt{c(cov.names, cov.names,binary)} is used as the covariates.
#' @param options.F a list of arguments in est_F() that are used for specific methods, including trim, formula.survPen, df, nfolds, s, alpha, lambda  
#' @param options.G a list of arguments in est_G() that are used for specific methods, including trim, formula.survPen, df, nfolds, s, alpha, lambda  
#' @param options.PS a list of arguments in est_PS() that are used for specific methods, including trim, df
#' @param options.Sd a list of arguments in est_Sd() that are used for specific methods, including trim, mtry, ntree, formula.survPen, df, nfolds, s, alpha, lambda  

cf_truncAC_AIPW <- function(dat, K, nu, X.name, Q.name, event.name, A.name, 
                            model.T, model.Q, model.A, model.D, 
                            cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim = 1e-7,
                            cov.names.binary.T = NULL, cov.names.binary.Q = NULL, cov.names.binary.A = NULL, cov.names.binary.D = NULL,
                            options.F = NULL, options.G = NULL, 
                            options.PS = NULL, options.Sd = NULL){
    
    
    dat_A1 = dat
    dat_A1[,A.name] <- 1
    dat_A0 = dat
    dat_A0[,A.name] <- 0
    
    names = c(X.name, Q.name, event.name, A.name, cov.names.T, cov.names.Q, cov.names.A, cov.names.D)
    if(sum(names == "delta.1")){
        stop("NThe names of the variables cannot be 'delta.1'. It is used in the middle of the computation.")
    }
    
    dat$delta.1 = rep(1, nrow(dat))
    event.name.Q = "delta.1"
    event.name.T = event.name
    
    n = nrow(dat)
    folds <- cvFolds(n, K)
    
    jumps.X = sort(dat[,X.name])
    jumps.Q = sort(dat[,Q.name])
    jumps.Y = sort(dat[,X.name] - dat[,Q.name])
    
    tau1 = min(jumps.X)
    tau2 = max(jumps.Q)
    tau3 = max(jumps.Y)
    u = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)
    v = c(tau1-1e-10, jumps.Q[jumps.Q>=tau1], max(jumps.Q)+1e-10)
    d = c(min(jumps.Y)-1e-10, jumps.Y, max(jumps.Y)+1e-10)
    
    tau.max = max(jumps.X)
    
    # For the truncAC_AIPW estimator
    Num.a1 = rep(NA, n)
    Num.a0 = rep(NA, n)
    Den = rep(NA, n)
    
    Num_truncA_IPCW.a1 = rep(NA, n)
    Num_truncA_IPCW.a0 = rep(NA, n)
    Den_truncA_IPCW = rep(NA, n)
    Num_IPW.a1 = rep(NA, n)
    Num_IPW.a0 = rep(NA, n)
    Den_IPW.a1 = rep(NA, n)
    Den_IPW.a0 = rep(NA, n)
    
    ct = 0
    for(k in 1:K){
        dat.est = dat[folds$subsets[folds$which == k], ]
        dat.fit = dat[folds$subsets[folds$which != k], ]
        dat_A1.est = dat_A1[folds$subsets[folds$which == k], ]
        dat_A0.est = dat_A0[folds$subsets[folds$which == k], ]
        
        ### Estimate the nuisance parameters
        # Estimate F and Sd
        Fuz.mx = est_F(dat.fit, dat.est, u, model.T, X.name, Q.name, event.name.T, 
                       cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                       trim = options.F$trim, 
                       mtry = options.F$mtry, ntree = options.F$ntree, 
                       formula.survPen = options.F$formula.survPen, 
                       nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                       lambda = options.F$lambda, df = options.F$df)

        Sdz.mx = est_Sd(dat.fit, dat.est, d, model.D, X.name, Q.name, event.name.T, 
                        cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                        trim = options.Sd$trim, 
                        mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                        formula.survPen = options.Sd$formula.survPen, 
                        nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                        lambda = options.Sd$lambda, df = options.Sd$df)   # d taking values on the residual time scale
        
        # Estimate G using IPCW weights
        id1.fit = (dat.fit[,event.name.T]==1)
        dat.fit.1 = dat.fit[id1.fit, ]  # uncensored data used to fit
        
        X.res = dat.fit[,X.name] - dat.fit[,Q.name]
        X.res.1 = X.res[id1.fit]
        
        Sd.vec.1 = pmax(diag(est_Sd(dat.fit, dat.fit, X.res.1, model.D, 
                                    X.name, Q.name, event.name.T,
                                    cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                    trim = options.Sd$trim, 
                                    mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                    formula.survPen = options.Sd$formula.survPen, 
                                    nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                    lambda = options.Sd$lambda, df = options.Sd$df,
                                    OOF = TRUE, nfolds.OOF = options.Sd$nfolds.OOF)[id1.fit,] ), trim)
        w_IPCW.1 = 1/Sd.vec.1
        
        Gvz.mx = est_G(dat.fit.1, dat.est, v, model.Q, X.name, Q.name, event.name.Q, 
                       cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                       weights = w_IPCW.1, tau = tau.max,
                       trim = options.G$trim,
                       formula.survPen = options.G$formula.survPen, 
                       nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                       lambda = options.G$lambda, df = options.G$df)
        
        
        # Estimate PS using IPCW and IPQW weights
        Gxz.vec.1 = pmax(diag(est_G(dat.fit.1, dat.fit.1, dat.fit.1$X, model.Q, X.name, Q.name, event.name.Q, 
                                    cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                    weights = w_IPCW.1, tau = tau.max,
                                    trim = options.G$trim,  
                                    formula.survPen = options.G$formula.survPen, 
                                    nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                    lambda = options.G$lambda, df = options.G$df)), trim)
        w_IPCW_IPQW.1 = 1/(Sd.vec.1*Gxz.vec.1)
        # summary(w_IPCW_IPQW.1)
        PS = est_PS(dat.fit.1, dat.est, model.A, A.name, 
                    cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                    weights = w_IPCW_IPQW.1,
                    trim = options.PS$trim, df = options.PS$df,
                    ntree = options.PS$ntree)   # a vector of estimated PS for each individual
        
        # compute the mu_1 and mu_0 for the treatment augmentation
        u.T = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)  # need the values of F in all jumps of T to compute mu()
        
        Fuz_A1.mx = est_F(dat.fit, dat_A1.est, u.T, model.T, X.name, Q.name, event.name.T, 
                          cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                          trim = options.F$trim, 
                          mtry = options.F$mtry, ntree = options.F$ntree, 
                          formula.survPen = options.F$formula.survPen, 
                          nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                          lambda = options.F$lambda, df = options.F$df)
        Fuz_A0.mx = est_F(dat.fit, dat_A0.est, u.T, model.T, X.name, Q.name, event.name.T, 
                          cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                          trim = options.F$trim, 
                          mtry = options.F$mtry, ntree = options.F$ntree, 
                          formula.survPen = options.F$formula.survPen, 
                          nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                          lambda = options.F$lambda, df = options.F$df)
        
        nuuT.mx = matrix(rep(nu(u.T), nrow(dat.est)), nrow = nrow(dat.est), byrow = TRUE)
        tau.Tmax = max(u.T) + 1e-5
        mu_A1 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A1.mx))
        mu_A0 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A0.mx))
        
        truncAC_EF = truncAC_AIPW_transMean_EF(dat.est, nu, Fuz.mx, Gvz.mx, Sdz.mx, PS, mu_A1, mu_A0, 
                                               X.name, Q.name, event.name.T, A.name, trim)
        
        nk = nrow(dat.est)
        Num.a1[(ct+1):(ct+nk)] = truncAC_EF$Num_truncAC_AIPW.a1
        Num.a0[(ct+1):(ct+nk)] = truncAC_EF$Num_truncAC_AIPW.a0
        Den[(ct+1):(ct+nk)] = truncAC_EF$Den_truncAC_AIPW
        
        Num_truncA_IPCW.a1[(ct+1):(ct+nk)] = truncAC_EF$Num_truncA_AIPW_IPCW.a1
        Num_truncA_IPCW.a0[(ct+1):(ct+nk)] = truncAC_EF$Num_truncA_AIPW_IPCW.a0
        Den_truncA_IPCW[(ct+1):(ct+nk)] = truncAC_EF$Den_truncA_AIPW_IPCW
        
        Num_IPW.a1[(ct+1):(ct+nk)] = truncAC_EF$Num_IPW.a1
        Num_IPW.a0[(ct+1):(ct+nk)] = truncAC_EF$Num_IPW.a0
        Den_IPW.a1[(ct+1):(ct+nk)] = truncAC_EF$Den_IPW.a1
        Den_IPW.a0[(ct+1):(ct+nk)] = truncAC_EF$Den_IPW.a0
        
        ct = ct + nk
    }
    
    
    est_truncAC = c(est.a1 = mean(Num.a1)/mean(Den),
                    est.a0 = mean(Num.a0)/mean(Den))
    
    est_truncA_IPCW = c(est.a1 = mean(Num_truncA_IPCW.a1)/mean(Den_truncA_IPCW),
                        est.a0 = mean(Num_truncA_IPCW.a0)/mean(Den_truncA_IPCW))
    
    est_IPW = c(est.a1 = mean(Num_IPW.a1)/mean(Den_IPW.a1),
                est.a0 = mean(Num_IPW.a0)/mean(Den_IPW.a0))
    
    return(list(truncAC_AIPW = est_truncAC,
                truncA_AIPW_IPCW = est_truncA_IPCW,
                IPW = est_IPW))
    
}





cf_truncAC_AIPW_boot <- function(dat, id, K, nu, X.name, Q.name, event.name, A.name, 
                            model.T, model.Q, model.A, model.D, 
                            cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim,
                            cov.names.binary.T, cov.names.binary.QL, cov.names.binary.A, cov.names.binary.D,
                            options.F, options.G, 
                            options.PS, options.Sd){
    
    dat = dat[id, ]
    
    est_truncAC_AIPW = cf_truncAC_AIPW(dat, K, nu, X.name, Q.name, event.name, A.name, 
                                       model.T, model.Q, model.A, model.D, 
                                       cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim,
                                       cov.names.binary.T, cov.names.binary.Q, cov.names.binary.A, cov.names.binary.D,
                                       options.F, options.G, 
                                       options.PS, options.Sd)
    
    est = rbind(truncAC_AIPW = est_truncAC_AIPW$truncAC_AIPW,
                truncA_AIPW_IPCW = est_truncAC_AIPW$truncA_AIPW_IPCW,
                IPW = est_truncAC_AIPW$IPW)
    est.a1 = est[,1]
    est.a0 = est[,2]
    
    est.a1a0 = c(est.a1, est.a0)
    names(est.a1a0) = c(paste(rownames(est), ".a1", sep = ""),
                        paste(rownames(est), ".a0", sep = ""))
    
    return(est.a1a0)
}




## Use another sample splitting for estimating the weights involved in nuisance estimation and the nuisance estimation
## the cross fitted truncAC_AIPW estimator

#' @param cov.names.binary.T vector of names for binary covariates if \texttt{model = "pCox"} is used. For other models, \texttt{c(cov.names, cov.names,binary)} is used as the covariates.
#' @param options.F a list of arguments in est_F() that are used for specific methods, including trim, formula.survPen, df, nfolds, s, alpha, lambda  
#' @param options.G a list of arguments in est_G() that are used for specific methods, including trim, formula.survPen, df, nfolds, s, alpha, lambda, trunc.weighting, where 'trunc.weighting' is a logical value indicating whether truncation weights are use to estimate G; if FALSE, G is estimated from uncensored subjects with IPCW weights 
#' @param options.PS a list of arguments in est_PS() that are used for specific methods, including trim, df. Currently PS is always estimated using truncating weighting.
#' @param options.Sd a list of arguments in est_Sd() that are used for specific methods, including trim, mtry, ntree, formula.survPen, df, nfolds, s, alpha, lambda  

cf2_truncAC_AIPW <- function(dat, K, nu, X.name, Q.name, event.name, A.name, 
                             model.T, model.Q, model.A, model.D, 
                             cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim = 1e-7,
                             cov.names.binary.T = NULL, cov.names.binary.Q = NULL, cov.names.binary.A = NULL, cov.names.binary.D = NULL,
                             options.F = NULL, options.G = NULL, 
                             options.PS = NULL, options.Sd = NULL){
    
    
    dat_A1 = dat
    dat_A1[,A.name] <- 1
    dat_A0 = dat
    dat_A0[,A.name] <- 0
    
    names = c(X.name, Q.name, event.name, A.name, cov.names.T, cov.names.Q, cov.names.A, cov.names.D)
    if(sum(names == "delta.1")){
        stop("The names of the variables cannot be 'delta.1'. It is used in the middle of the computation.")
    }
    
    dat$delta.1 = rep(1, nrow(dat))
    event.name.Q = "delta.1"
    event.name.T = event.name
    
    n = nrow(dat)
    folds <- cvFolds(n, K)
    
    jumps.X = sort(dat[,X.name])
    jumps.Q = sort(dat[,Q.name])
    jumps.Y = sort(dat[,X.name] - dat[,Q.name])
    
    tau1 = min(jumps.X)
    tau2 = max(jumps.Q)
    tau3 = max(jumps.Y)
    u = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)
    v = c(tau1-1e-10, jumps.Q[jumps.Q>=tau1], max(jumps.Q)+1e-10)
    d = c(min(jumps.Y)-1e-10, jumps.Y, max(jumps.Y)+1e-10)
    
    tau.max = max(jumps.X)
    
    # For the truncAC_AIPW estimator
    Num.a1 = rep(NA, n)
    Num.a0 = rep(NA, n)
    Den = rep(NA, n)
    
    Num_truncA_IPCW.a1 = rep(NA, n)
    Num_truncA_IPCW.a0 = rep(NA, n)
    Den_truncA_IPCW = rep(NA, n)
    Num_IPW.a1 = rep(NA, n)
    Num_IPW.a0 = rep(NA, n)
    Den_IPW.a1 = rep(NA, n)
    Den_IPW.a0 = rep(NA, n)
    
    ct = 0
    for(k in 1:K){
        dat.est = dat[folds$subsets[folds$which == k], ]
        dat.fit = dat[folds$subsets[folds$which != k], ]
        dat_A1.est = dat_A1[folds$subsets[folds$which == k], ]
        dat_A0.est = dat_A0[folds$subsets[folds$which == k], ]
        
        
        # split the fit data into 4 folds
        n.fit = nrow(dat.fit)
        folds.fit <- cvFolds(n.fit, 2)
        dat.fit.s1 = dat.fit[folds.fit$subsets[folds.fit$which == 1], ]
        dat.fit.s2 = dat.fit[folds.fit$subsets[folds.fit$which == 2], ]
        
        ### Estimate the nuisance parameters
        
        ## truncation weights 1/(1-F(Q|A,Z))
        Fq.vec.PS = diag(est_F(dat.fit.s1, dat.fit.s2, dat.fit.s2[,Q.name], 
                               model.T, X.name, Q.name, event.name.T, 
                               cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                               trim = options.F$trim, 
                               mtry = options.F$mtry, ntree = options.F$ntree, 
                               formula.survPen = options.F$formula.survPen, 
                               nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                               lambda = options.F$lambda, df = options.F$df))
        w_trunc.F = 1/pmax(1-Fq.vec.PS, trim)
        # summary(w_trunc.F)
        
        
        ## Estimate F and Sd
        Fuz.mx = est_F(dat.fit.s2, dat.est, u, model.T, X.name, Q.name, event.name.T, 
                       cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                       trim = options.F$trim, 
                       mtry = options.F$mtry, ntree = options.F$ntree, 
                       formula.survPen = options.F$formula.survPen, 
                       nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                       lambda = options.F$lambda, df = options.F$df)
        
        Sdz.mx = est_Sd(dat.fit.s2, dat.est, d, model.D, X.name, Q.name, event.name.T, 
                        cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                        trim = options.Sd$trim, 
                        mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                        formula.survPen = options.Sd$formula.survPen, 
                        nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                        lambda = options.Sd$lambda, df = options.Sd$df)   # d taking values on the residual time scale
        
        ## Estimate G
        if(is.null(options.G$trunc.weighting)|options.G$trunc.weighting){
            
            # Estimate G using truncation weights 1/(1-F(Q|A,Z))
            Gvz.mx = est_G(dat.fit.s2, dat.est, v, model.Q, X.name, Q.name, event.name.Q, 
                           cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                           weights = w_trunc.F, trunc = FALSE, 
                           tau = tau.max, trim = options.G$trim,
                           formula.survPen = options.G$formula.survPen, 
                           nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                           lambda = options.G$lambda, df = options.G$df)
            
            
        }else{
            # Estimate G using uncensored subjects with IPCW weights
            id1.fit.s2 = (dat.fit.s2[,event.name.T]==1)
            dat.fit.s2_1 = dat.fit.s2[id1.fit.s2, ]  # uncensored data used to fit G
            X.res.s2_1 = dat.fit.s2_1[,X.name] - dat.fit.s2_1[,Q.name]
            
            Sd.vec.G1 = diag(est_Sd(dat.fit.s1, dat.fit.s2_1, X.res.s2_1, 
                                    model.D, X.name, Q.name, event.name.T,
                                    cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                    trim = options.Sd$trim, 
                                    mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                    formula.survPen = options.Sd$formula.survPen, 
                                    nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                    lambda = options.Sd$lambda, df = options.Sd$df))
            
            w_IPCW.G1 = 1/pmax(Sd.vec.G1, trim)
            # summary(w_IPCW.G1)
            # length(w_IPCW.G1)
            Gvz.mx = est_G(dat.fit.s2_1, dat.est, v, model.Q, X.name, Q.name, event.name.Q, 
                           cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                           weights = w_IPCW.G1, tau = tau.max,
                           trim = options.G$trim,
                           formula.survPen = options.G$formula.survPen, 
                           nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                           lambda = options.G$lambda, df = options.G$df)
            
            # warning: algorithm does not converge with n = 999 
            # - check it converges when n = 1999 - yes for most, but still has a few that does not converge
        }
        
        # Estimate PS using truncation weights 1/(1-F(Q|A,Z))
        PS = est_PS(dat.fit.s2, dat.est, model.A, A.name, 
                    cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                    weights = w_trunc.F,
                    trim = options.PS$trim, df = options.PS$df,
                    ntree = options.PS$ntree)   # a vector of estimated PS for each individual
        
        # compute the mu_1 and mu_0 for the treatment augmentation
        u.T = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)  # need the values of F in all jumps of T to compute mu()
        
        Fuz_A1.mx = est_F(dat.fit.s2, dat_A1.est, u.T, model.T, X.name, Q.name, event.name.T, 
                          cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                          trim = options.F$trim, 
                          mtry = options.F$mtry, ntree = options.F$ntree, 
                          formula.survPen = options.F$formula.survPen, 
                          nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                          lambda = options.F$lambda, df = options.F$df)
        Fuz_A0.mx = est_F(dat.fit.s2, dat_A0.est, u.T, model.T, X.name, Q.name, event.name.T, 
                          cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                          trim = options.F$trim, 
                          mtry = options.F$mtry, ntree = options.F$ntree, 
                          formula.survPen = options.F$formula.survPen, 
                          nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                          lambda = options.F$lambda, df = options.F$df)
        
        nuuT.mx = matrix(rep(nu(u.T), nrow(dat.est)), nrow = nrow(dat.est), byrow = TRUE)
        tau.Tmax = max(u.T) + 1e-5
        mu_A1 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A1.mx))
        mu_A0 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A0.mx))
        
        truncAC_EF = truncAC_AIPW_transMean_EF(dat.est, nu, Fuz.mx, Gvz.mx, Sdz.mx, PS, mu_A1, mu_A0, 
                                               X.name, Q.name, event.name.T, A.name, trim)
        
        nk = nrow(dat.est)
        Num.a1[(ct+1):(ct+nk)] = truncAC_EF$Num_truncAC_AIPW.a1
        Num.a0[(ct+1):(ct+nk)] = truncAC_EF$Num_truncAC_AIPW.a0
        Den[(ct+1):(ct+nk)] = truncAC_EF$Den_truncAC_AIPW
        
        Num_truncA_IPCW.a1[(ct+1):(ct+nk)] = truncAC_EF$Num_truncA_AIPW_IPCW.a1
        Num_truncA_IPCW.a0[(ct+1):(ct+nk)] = truncAC_EF$Num_truncA_AIPW_IPCW.a0
        Den_truncA_IPCW[(ct+1):(ct+nk)] = truncAC_EF$Den_truncA_AIPW_IPCW
        
        Num_IPW.a1[(ct+1):(ct+nk)] = truncAC_EF$Num_IPW.a1
        Num_IPW.a0[(ct+1):(ct+nk)] = truncAC_EF$Num_IPW.a0
        Den_IPW.a1[(ct+1):(ct+nk)] = truncAC_EF$Den_IPW.a1
        Den_IPW.a0[(ct+1):(ct+nk)] = truncAC_EF$Den_IPW.a0
        
        ct = ct + nk
    }
    
    
    est_truncAC = c(est.a1 = mean(Num.a1)/mean(Den),
                    est.a0 = mean(Num.a0)/mean(Den))
    
    est_truncA_IPCW = c(est.a1 = mean(Num_truncA_IPCW.a1)/mean(Den_truncA_IPCW),
                        est.a0 = mean(Num_truncA_IPCW.a0)/mean(Den_truncA_IPCW))
    
    est_IPW = c(est.a1 = mean(Num_IPW.a1)/mean(Den_IPW.a1),
                est.a0 = mean(Num_IPW.a0)/mean(Den_IPW.a0))
    
    return(list(truncAC_AIPW = est_truncAC,
                truncA_AIPW_IPCW = est_truncA_IPCW,
                IPW = est_IPW))
    
}





cf2_truncAC_AIPW_boot <- function(dat, id, K, nu, X.name, Q.name, event.name, A.name, 
                                  model.T, model.Q, model.A, model.D, 
                                  cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim,
                                  cov.names.binary.T, cov.names.binary.QL, cov.names.binary.A, cov.names.binary.D,
                                  options.F, options.G, 
                                  options.PS, options.Sd){
    
    dat = dat[id, ]
    
    est_truncAC_AIPW = cf2_truncAC_AIPW(dat, K, nu, X.name, Q.name, event.name, A.name, 
                                        model.T, model.Q, model.A, model.D, 
                                        cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim,
                                        cov.names.binary.T, cov.names.binary.Q, cov.names.binary.A, cov.names.binary.D,
                                        options.F, options.G, 
                                        options.PS, options.Sd)
    
    est = rbind(truncAC_AIPW = est_truncAC_AIPW$truncAC_AIPW,
                truncA_AIPW_IPCW = est_truncAC_AIPW$truncA_AIPW_IPCW,
                IPW = est_truncAC_AIPW$IPW)
    est.a1 = est[,1]
    est.a0 = est[,2]
    
    est.a1a0 = c(est.a1, est.a0)
    names(est.a1a0) = c(paste(rownames(est), ".a1", sep = ""),
                        paste(rownames(est), ".a0", sep = ""))
    
    return(est.a1a0)
}










