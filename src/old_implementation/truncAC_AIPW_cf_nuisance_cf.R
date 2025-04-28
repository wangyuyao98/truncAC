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
        
        # split the fit data into 2 folds
        n.fit = nrow(dat.fit)
        folds.fit <- cvFolds(n.fit, 2)
        id.fit.s1  = folds.fit$subsets[folds.fit$which == 1]
        id.fit.s2  = folds.fit$subsets[folds.fit$which == 2]
        dat.fit.s1 = dat.fit[id.fit.s1, ]
        dat.fit.s2 = dat.fit[id.fit.s2, ]
        
        # split the fit data into 3 folds - for one approach of PS estimation
        folds.fit.PS <- cvFolds(n.fit, 3)
        id.fit.ss1  = folds.fit.PS$subsets[folds.fit.PS$which == 1]
        id.fit.ss2  = folds.fit.PS$subsets[folds.fit.PS$which == 2]
        id.fit.ss3  = folds.fit.PS$subsets[folds.fit.PS$which == 3]
        dat.fit.ss1 = dat.fit[id.fit.ss1, ]
        dat.fit.ss2 = dat.fit[id.fit.ss2, ]
        dat.fit.ss3 = dat.fit[id.fit.ss3, ]
        
        
        
        ### Estimate the nuisance parameters
        
        ## Estimate F and F_A1 F_A0
        if(is.null(Fuz.mx)){
            u = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)
            Fuz.mx.est = est_F(dat.fit, dat.est, u, model.T, X.name, Q.name, event.name.T, 
                               cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                               trim = options.F$trim, 
                               mtry = options.F$mtry, ntree = options.F$ntree, 
                               formula.survPen = options.F$formula.survPen, 
                               nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                               lambda = options.F$lambda, df = options.F$df)
        }else{
            Fuz.mx.est = Fuz.mx[id.est, ]
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
        }
        
        
        
        ## Estimate G - using truncation weights 1/{1-F(Q|A,Z)} estimated using fit.si, and use fit.sj to estimate G, i\neq j \in \{1,2\} 
        v = c(tau1-1e-10, jumps.Q[jumps.Q>=tau1], max(jumps.Q)+1e-10)
        if(is.null(Gvz.mx)){
            if(est_approach_G == "truncIPW.F"){
                Fq.vec.fit.s1 = NULL
                Fq.vec.fit.s2 = NULL
                
                # Compute the truncation weights 1/{1-F(Q|A,Z)} 
                if(is.null(Fuz.mx)){
                    Fq.vec.fit.s2 = diag(est_F(dat.fit.s1, dat.fit.s2, dat.fit.s2[,Q.name], 
                                               model.T, X.name, Q.name, event.name.T, 
                                               cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                                               trim = options.F$trim, 
                                               mtry = options.F$mtry, ntree = options.F$ntree, 
                                               formula.survPen = options.F$formula.survPen, 
                                               nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                               lambda = options.F$lambda, df = options.F$df))
                    
                    Fq.vec.fit.s1 = diag(est_F(dat.fit.s2, dat.fit.s1, dat.fit.s1[,Q.name], 
                                               model.T, X.name, Q.name, event.name.T, 
                                               cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                                               trim = options.F$trim, 
                                               mtry = options.F$mtry, ntree = options.F$ntree, 
                                               formula.survPen = options.F$formula.survPen, 
                                               nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                               lambda = options.F$lambda, df = options.F$df))
                    
                }else{
                    Fq.vec.fit.s2 = CDF_eval(dat.fit.s2[ ,Q.name], Fuz.mx[id.fit,][id.fit.s2, ])
                    Fq.vec.fit.s1 = CDF_eval(dat.fit.s1[ ,Q.name], Fuz.mx[id.fit,][id.fit.s1, ])
                }
                
                
                w_truncF.fit.s2 = 1/pmax(1-Fq.vec.fit.s2, trim)
                w_truncF.fit.s1 = 1/pmax(1-Fq.vec.fit.s1, trim)
                
                # Estimate G using fit.s2
                Gvz.mx.est.s2 = est_G(dat.fit.s2, dat.est, v, model.Q, X.name, Q.name, event.name.Q, 
                                      cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                      weights = w_truncF.fit.s2, trunc = FALSE, 
                                      tau = tau.max, trim = options.G$trim,
                                      formula.survPen = options.G$formula.survPen, 
                                      nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                      lambda = options.G$lambda, df = options.G$df)
                # Estimate G using fit.s1
                Gvz.mx.est.s1 = est_G(dat.fit.s1, dat.est, v, model.Q, X.name, Q.name, event.name.Q, 
                                      cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                      weights = w_truncF.fit.s1, trunc = FALSE, 
                                      tau = tau.max, trim = options.G$trim,
                                      formula.survPen = options.G$formula.survPen, 
                                      nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                      lambda = options.G$lambda, df = options.G$df)
                # Average the two estimates
                Gvz.mx.est = (Gvz.mx.est.s2 + Gvz.mx.est.s1)/2
                
            }else if(est_approach_G == "IPCW"){
                
                id.fit.s1_1 = (dat.fit.s1[, event.name.T] == 1)
                id.fit.s2_1 = (dat.fit.s2[, event.name.T] == 1)
                dat.fit.s1_1 = dat.fit.s1[id.fit.s1_1, ]
                dat.fit.s2_1 = dat.fit.s2[id.fit.s2_1, ]
                X.res.s1_1 = dat.fit.s1_1[,X.name] - dat.fit.s1_1[,Q.name]
                X.res.s2_1 = dat.fit.s2_1[,X.name] - dat.fit.s2_1[,Q.name]
                
                # Estimate the censoring weights
                if(is.null(Sdz.mx)){
                    Sdy.vec.fit.s2_1 = diag(est_Sd(dat.fit.s1, dat.fit.s2_1, X.res.s2_1, 
                                                   model.D, X.name, Q.name, event.name.T,
                                                   cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                                   trim = options.Sd$trim, 
                                                   mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                                   formula.survPen = options.Sd$formula.survPen, 
                                                   nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                                   lambda = options.Sd$lambda, df = options.Sd$df))
                    Sdy.vec.fit.s1_1 = diag(est_Sd(dat.fit.s2, dat.fit.s1_1, X.res.s1_1, 
                                                   model.D, X.name, Q.name, event.name.T,
                                                   cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                                   trim = options.Sd$trim, 
                                                   mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                                   formula.survPen = options.Sd$formula.survPen, 
                                                   nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                                   lambda = options.Sd$lambda, df = options.Sd$df))
                }else{
                    Sdy.vec.fit.s2_1 = 1- CDF_eval(X.res.s2_1, 1-Sdz.mx[id.fit,][id.fit.s2,][id.fit.s2_1,])
                    Sdy.vec.fit.s1_1 = 1- CDF_eval(X.res.s1_1, 1-Sdz.mx[id.fit,][id.fit.s1,][id.fit.s1_1,])
                }
                
                w_IPCW_s1_1 = 1/pmax(Sdy.vec.fit.s1_1, trim)
                w_IPCW_s2_1 = 1/pmax(Sdy.vec.fit.s2_1, trim)
                
                Gvz.mx.est.s2 = est_G(dat.fit.s2_1, dat.est, v, model.Q, X.name, Q.name, event.name.Q, 
                                      cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                      weights = w_IPCW_s2_1, tau = tau.max,
                                      trim = options.G$trim,
                                      formula.survPen = options.G$formula.survPen, 
                                      nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                      lambda = options.G$lambda, df = options.G$df)
                Gvz.mx.est.s1 = est_G(dat.fit.s1_1, dat.est, v, model.Q, X.name, Q.name, event.name.Q, 
                                      cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                      weights = w_IPCW_s1_1, tau = tau.max,
                                      trim = options.G$trim,
                                      formula.survPen = options.G$formula.survPen, 
                                      nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                      lambda = options.G$lambda, df = options.G$df)
                
                # Average the two estimates
                Gvz.mx.est = (Gvz.mx.est.s2 + Gvz.mx.est.s1)/2
                
                
            }else{
                stop("'est_approach_G' Should be either 'truncIPW.F' or 'IPCW'. ")
            }
            
            
        }else{
            Gvz.mx.est = Gvz.mx[id.est, ]
        }
        
        
        
        ## Estimate PS
        if(is.null(PS)){
            
            if(est_approach_PS == "truncIPW.F"){
                
                if(is.null(Fq.vec.fit.s1)|is.null(Fq.vec.fit.s2)){ # The truncation weights from F haven't been computed yet
                    
                    # Compute the truncation weights 1/{1-F(Q|A,Z)} estimated using fit.s1
                    if(is.null(Fuz.mx)){
                        Fq.vec.fit.s2 = diag(est_F(dat.fit.s1, dat.fit.s2, dat.fit.s2[,Q.name], 
                                                   model.T, X.name, Q.name, event.name.T, 
                                                   cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                                                   trim = options.F$trim, 
                                                   mtry = options.F$mtry, ntree = options.F$ntree, 
                                                   formula.survPen = options.F$formula.survPen, 
                                                   nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                                   lambda = options.F$lambda, df = options.F$df))
                        
                        Fq.vec.fit.s1 = diag(est_F(dat.fit.s2, dat.fit.s1, dat.fit.s1[,Q.name], 
                                                   model.T, X.name, Q.name, event.name.T, 
                                                   cov.names = cov.names.T, cov.names.binary = cov.names.binary.T, 
                                                   trim = options.F$trim, 
                                                   mtry = options.F$mtry, ntree = options.F$ntree, 
                                                   formula.survPen = options.F$formula.survPen, 
                                                   nfolds = options.F$nfolds, s = options.F$s, alpha = options.F$alpha,
                                                   lambda = options.F$lambda, df = options.F$df))
                    }else{
                        Fq.vec.fit.s2 = CDF_eval(dat.fit.s2[ ,Q.name], Fuz.mx[id.fit,][id.fit.s2, ])
                        Fq.vec.fit.s1 = CDF_eval(dat.fit.s1[ ,Q.name], Fuz.mx[id.fit,][id.fit.s1, ])
                    }
                    
                    w_truncF.fit.s2 = 1/pmax(1-Fq.vec.fit.s2, trim)
                    w_truncF.fit.s1 = 1/pmax(1-Fq.vec.fit.s1, trim)
                } # End if(is.null(Fq.vec.fit.s1)|is.null(Fq.vec.fit.s2)
                
                # Estimate PS using truncation weights 1/(1-F(Q|A,Z))
                PS.est.s2 = est_PS(dat.fit.s2, dat.est, model.A, A.name, 
                                   cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                   weights = w_truncF.fit.s2,
                                   trim = options.PS$trim, df = options.PS$df,
                                   ntree = options.PS$ntree)   # a vector of estimated PS for each individual
                PS.est.s1 = est_PS(dat.fit.s1, dat.est, model.A, A.name, 
                                   cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                   weights = w_truncF.fit.s1,
                                   trim = options.PS$trim, df = options.PS$df,
                                   ntree = options.PS$ntree)   # a vector of estimated PS for each individual
                PS.est = (PS.est.s1 + PS.est.s2)/2
                
                
                
            }else if(est_approach_PS == "truncIPCW"){
                
                id.fit.ss1_1 = (dat.fit.ss1[, event.name.T] == 1)
                id.fit.ss2_1 = (dat.fit.ss2[, event.name.T] == 1)
                id.fit.ss3_1 = (dat.fit.ss3[, event.name.T] == 1)
                dat.fit.ss1_1 = dat.fit.ss1[id.fit.ss1_1, ]
                dat.fit.ss2_1 = dat.fit.ss2[id.fit.ss2_1, ]
                dat.fit.ss3_1 = dat.fit.ss3[id.fit.ss3_1, ]
                X.res.ss1_1 = dat.fit.ss1_1[,X.name] - dat.fit.ss1_1[,Q.name]
                X.res.ss2_1 = dat.fit.ss2_1[,X.name] - dat.fit.ss2_1[,Q.name]
                X.res.ss3_1 = dat.fit.ss3_1[,X.name] - dat.fit.ss3_1[,Q.name]
                X.ss1_1 = dat.fit.ss1_1[,X.name] 
                X.ss2_1 = dat.fit.ss2_1[,X.name] 
                X.ss3_1 = dat.fit.ss3_1[,X.name] 
                
                # 6 combinations in total for how dat.fit is used to estimate the nuisance.
                # Compute the IPCW weights
                if(is.null(Sdz.mx)){
                    Sdy.vec.fit.ss1_2 = diag(est_Sd(dat.fit.ss2, dat.fit.ss1_1, X.res.ss1_1, 
                                                    model.D, X.name, Q.name, event.name.T,
                                                    cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                                    trim = options.Sd$trim, 
                                                    mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                                    formula.survPen = options.Sd$formula.survPen, 
                                                    nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                                    lambda = options.Sd$lambda, df = options.Sd$df))
                    Sdy.vec.fit.ss1_3 = diag(est_Sd(dat.fit.ss3, dat.fit.ss1_1, X.res.ss1_1, 
                                                    model.D, X.name, Q.name, event.name.T,
                                                    cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                                    trim = options.Sd$trim, 
                                                    mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                                    formula.survPen = options.Sd$formula.survPen, 
                                                    nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                                    lambda = options.Sd$lambda, df = options.Sd$df))
                    
                    Sdy.vec.fit.ss2_1 = diag(est_Sd(dat.fit.ss1, dat.fit.ss2_1, X.res.ss2_1, 
                                                    model.D, X.name, Q.name, event.name.T,
                                                    cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                                    trim = options.Sd$trim, 
                                                    mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                                    formula.survPen = options.Sd$formula.survPen, 
                                                    nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                                    lambda = options.Sd$lambda, df = options.Sd$df))
                    Sdy.vec.fit.ss2_3 = diag(est_Sd(dat.fit.ss3, dat.fit.ss2_1, X.res.ss2_1, 
                                                    model.D, X.name, Q.name, event.name.T,
                                                    cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                                    trim = options.Sd$trim, 
                                                    mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                                    formula.survPen = options.Sd$formula.survPen, 
                                                    nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                                    lambda = options.Sd$lambda, df = options.Sd$df))
                    
                    Sdy.vec.fit.ss3_1 = diag(est_Sd(dat.fit.ss1, dat.fit.ss3_1, X.res.ss3_1, 
                                                    model.D, X.name, Q.name, event.name.T,
                                                    cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                                    trim = options.Sd$trim, 
                                                    mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                                    formula.survPen = options.Sd$formula.survPen, 
                                                    nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                                    lambda = options.Sd$lambda, df = options.Sd$df))
                    Sdy.vec.fit.ss3_2 = diag(est_Sd(dat.fit.ss2, dat.fit.ss3_1, X.res.ss3_1, 
                                                    model.D, X.name, Q.name, event.name.T,
                                                    cov.names = cov.names.D, cov.names.binary = cov.names.binary.D, 
                                                    trim = options.Sd$trim, 
                                                    mtry = options.Sd$mtry, ntree = options.Sd$ntree, 
                                                    formula.survPen = options.Sd$formula.survPen, 
                                                    nfolds = options.Sd$nfolds, s = options.Sd$s, alpha = options.Sd$alpha,
                                                    lambda = options.Sd$lambda, df = options.Sd$df))
                    
                }else{
                    Sdy.vec.fit.ss1 = 1- CDF_eval(X.res.ss1_1, 1-Sdz.mx[id.fit,][id.fit.ss1,][id.fit.ss1_1,])
                    Sdy.vec.fit.ss2 = 1- CDF_eval(X.res.ss2_1, 1-Sdz.mx[id.fit,][id.fit.ss2,][id.fit.ss2_1,])
                    Sdy.vec.fit.ss3 = 1- CDF_eval(X.res.ss3_1, 1-Sdz.mx[id.fit,][id.fit.ss3,][id.fit.ss3_1,])
                    
                    Sdy.vec.fit.ss1_2 = Sdy.vec.fit.ss1
                    Sdy.vec.fit.ss1_3 = Sdy.vec.fit.ss1
                    
                    Sdy.vec.fit.ss2_1 = Sdy.vec.fit.ss2
                    Sdy.vec.fit.ss2_3 = Sdy.vec.fit.ss2
                    
                    Sdy.vec.fit.ss3_1 = Sdy.vec.fit.ss3
                    Sdy.vec.fit.ss3_2 = Sdy.vec.fit.ss3
                }
                
                # IPCW weights
                w_IPCW.ss1_2 = 1/pmax(Sdy.vec.fit.ss1_2, trim)   # the IPCW weights for dat.fit.ss1_1 with censoring model fitted using dat.fit.ss2
                w_IPCW.ss1_3 = 1/pmax(Sdy.vec.fit.ss1_3, trim)
                
                w_IPCW.ss2_1 = 1/pmax(Sdy.vec.fit.ss2_1, trim)
                w_IPCW.ss2_3 = 1/pmax(Sdy.vec.fit.ss2_3, trim)
                
                w_IPCW.ss3_1 = 1/pmax(Sdy.vec.fit.ss3_1, trim)
                w_IPCW.ss3_2 = 1/pmax(Sdy.vec.fit.ss3_2, trim)
                
                if(is.null(Gvz.mx)){
                    Gx.vec.fit.ss1_23 = diag( est_G(dat.fit.ss2_1, dat.fit.ss1_1, X.ss1_1, model.Q, X.name, Q.name, event.name.Q, 
                                                    cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                                    weights = w_IPCW.ss2_3, tau = tau.max,
                                                    trim = options.G$trim,
                                                    formula.survPen = options.G$formula.survPen, 
                                                    nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                                    lambda = options.G$lambda, df = options.G$df) )
                    Gx.vec.fit.ss1_32 = diag( est_G(dat.fit.ss3_1, dat.fit.ss1_1, X.ss1_1, model.Q, X.name, Q.name, event.name.Q, 
                                                    cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                                    weights = w_IPCW.ss3_2, tau = tau.max,
                                                    trim = options.G$trim,
                                                    formula.survPen = options.G$formula.survPen, 
                                                    nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                                    lambda = options.G$lambda, df = options.G$df) )
                    
                    Gx.vec.fit.ss2_13 = diag( est_G(dat.fit.ss1_1, dat.fit.ss2_1, X.ss2_1, model.Q, X.name, Q.name, event.name.Q, 
                                                    cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                                    weights = w_IPCW.ss1_3, tau = tau.max,
                                                    trim = options.G$trim,
                                                    formula.survPen = options.G$formula.survPen, 
                                                    nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                                    lambda = options.G$lambda, df = options.G$df) )
                    Gx.vec.fit.ss2_31 = diag( est_G(dat.fit.ss3_1, dat.fit.ss2_1, X.ss2_1, model.Q, X.name, Q.name, event.name.Q, 
                                                    cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                                    weights = w_IPCW.ss3_1, tau = tau.max,
                                                    trim = options.G$trim,
                                                    formula.survPen = options.G$formula.survPen, 
                                                    nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                                    lambda = options.G$lambda, df = options.G$df) )
                    
                    Gx.vec.fit.ss3_12 = diag( est_G(dat.fit.ss1_1, dat.fit.ss3_1, X.ss3_1, model.Q, X.name, Q.name, event.name.Q, 
                                                    cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                                    weights = w_IPCW.ss1_2, tau = tau.max,
                                                    trim = options.G$trim,
                                                    formula.survPen = options.G$formula.survPen, 
                                                    nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                                    lambda = options.G$lambda, df = options.G$df) )
                    Gx.vec.fit.ss3_21 = diag( est_G(dat.fit.ss2_1, dat.fit.ss3_1, X.ss3_1, model.Q, X.name, Q.name, event.name.Q, 
                                                    cov.names = cov.names.Q, cov.names.binary = cov.names.binary.Q,
                                                    weights = w_IPCW.ss2_1, tau = tau.max,
                                                    trim = options.G$trim,
                                                    formula.survPen = options.G$formula.survPen, 
                                                    nfolds = options.G$nfolds, s = options.G$s, alpha = options.G$alpha,
                                                    lambda = options.G$lambda, df = options.G$df) )
                    
                    
                }else{
                    Gx.vec.fit.ss1 = CDF_eval(X.ss1_1, Gvz.mx[id.fit,][id.fit.ss1,][id.fit.ss1_1,])
                    Gx.vec.fit.ss2 = CDF_eval(X.ss2_1, Gvz.mx[id.fit,][id.fit.ss2,][id.fit.ss2_1,])
                    Gx.vec.fit.ss3 = CDF_eval(X.ss3_1, Gvz.mx[id.fit,][id.fit.ss3,][id.fit.ss3_1,])
                    
                    Gx.vec.fit.ss1_23 = Gx.vec.fit.ss1
                    Gx.vec.fit.ss1_32 = Gx.vec.fit.ss1
                    
                    Gx.vec.fit.ss2_13 = Gx.vec.fit.ss2
                    Gx.vec.fit.ss2_31 = Gx.vec.fit.ss2
                    
                    Gx.vec.fit.ss3_12 = Gx.vec.fit.ss3
                    Gx.vec.fit.ss3_21 = Gx.vec.fit.ss3
                }
                
                w_trunc.ss1_23 = 1/pmax(Gx.vec.fit.ss1_23, trim)
                w_trunc.ss1_32 = 1/pmax(Gx.vec.fit.ss1_32, trim)
                
                w_trunc.ss2_13 = 1/pmax(Gx.vec.fit.ss2_13, trim)
                w_trunc.ss2_31 = 1/pmax(Gx.vec.fit.ss2_31, trim)
                
                w_trunc.ss3_12 = 1/pmax(Gx.vec.fit.ss3_12, trim)
                w_trunc.ss3_21 = 1/pmax(Gx.vec.fit.ss3_21, trim)
                
                
                # Estimate PS using truncIPCW weights \Delta/{G(X|A,Z)S_D(X-Q|A,Z)}
                PS.ss1_23 = est_PS(dat.fit.ss1_1, dat.est, model.A, A.name, 
                                   cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                   weights = w_trunc.ss1_23 * w_IPCW.ss1_3,
                                   trim = options.PS$trim, df = options.PS$df,
                                   ntree = options.PS$ntree) 
                PS.ss1_32 = est_PS(dat.fit.ss1_1, dat.est, model.A, A.name, 
                                   cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                   weights = w_trunc.ss1_32 * w_IPCW.ss1_2,
                                   trim = options.PS$trim, df = options.PS$df,
                                   ntree = options.PS$ntree)
                
                PS.ss2_13 = est_PS(dat.fit.ss2_1, dat.est, model.A, A.name, 
                                   cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                   weights = w_trunc.ss2_13 * w_IPCW.ss2_3,
                                   trim = options.PS$trim, df = options.PS$df,
                                   ntree = options.PS$ntree)
                PS.ss2_31 = est_PS(dat.fit.ss2_1, dat.est, model.A, A.name, 
                                   cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                   weights = w_trunc.ss2_31 * w_IPCW.ss2_1,
                                   trim = options.PS$trim, df = options.PS$df,
                                   ntree = options.PS$ntree)
                
                PS.ss3_12 = est_PS(dat.fit.ss3_1, dat.est, model.A, A.name, 
                                   cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                   weights = w_trunc.ss3_12 * w_IPCW.ss3_2,
                                   trim = options.PS$trim, df = options.PS$df,
                                   ntree = options.PS$ntree)
                PS.ss3_21 = est_PS(dat.fit.ss3_1, dat.est, model.A, A.name, 
                                   cov.names = cov.names.A, cov.names.binary = cov.names.binary.A, 
                                   weights = w_trunc.ss3_21 * w_IPCW.ss3_1,
                                   trim = options.PS$trim, df = options.PS$df,
                                   ntree = options.PS$ntree)
                
                PS.est = (PS.ss1_23 + PS.ss1_32 + PS.ss2_13 + PS.ss2_31 + PS.ss3_12 + PS.ss3_21)/6
                
                
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
        mu_A1.est = as.vector(int_fmx_dF(tau.Tmax_A1, nuuT_A1.mx, Fuz_A1.mx.est))
        mu_A0.est = as.vector(int_fmx_dF(tau.Tmax_A0, nuuT_A0.mx, Fuz_A0.mx.est))
        
        
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
        GXZ_hat[id.est] = CDF_eval(dat.est[,X.name], Gvz.mx.est)
        SdXZ_hat[id.est] = 1 - CDF_eval(dat.est[,X.name] - dat.est[,Q.name], 1-Sdz.mx.est)
        
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






cf_truncAC_AIPW_bootSE <- function(seed.boot, n.boot, dat, K, nu, X.name, Q.name, event.name, A.name, 
                                   model.T, model.Q, model.A, model.D, 
                                   cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim = 1e-7,
                                   cov.names.binary.T = NULL, cov.names.binary.Q = NULL, cov.names.binary.A = NULL, cov.names.binary.D = NULL,
                                   options.F = NULL, options.G = NULL, 
                                   options.PS = NULL, options.Sd = NULL,
                                   est_approach_G = "truncIPW.F", est_approach_PS = "truncIPW.F",
                                   Fuz.mx = NULL, Fuz_A1.mx = NULL, Fuz_A0.mx = NULL,
                                   Gvz.mx = NULL, Sdz.mx = NULL, PS = NULL){
    
    set.seed(seed.boot)
    bootresult = boot(dat, cf_truncAC_AIPW_boot, R = n.boot,
                      K = K, nu = nu, X.name = X.name, Q.name = Q.name, event.name = event.name, A.name = A.name,
                      model.T = model.T, model.Q = model.Q, model.D = model.D, model.A = model.A, 
                      cov.names.T = cov.names.T, cov.names.Q = cov.names.Q, 
                      cov.names.A = cov.names.A, cov.names.D = cov.names.D,
                      trim = trim,
                      cov.names.binary.T = cov.names.binary.T, cov.names.binary.Q = cov.names.binary.Q, 
                      cov.names.binary.A = cov.names.binary.A, cov.names.binary.D = cov.names.binary.D,
                      options.F = options.F, options.G = options.G, 
                      options.PS = options.PS, options.Sd = options.Sd)
    
    return(bootresult)
}






