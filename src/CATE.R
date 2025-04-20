## Functions for estimating CATE

CATE_truth <- function(dat, nn = 1e5, nu,
                      T.min, beta.T, shape.T){
    n = nrow(dat)
    
    ## The names of the covariates are Z1, Z2
    CATE.truth = rep(NA, n)
    for(i in 1:n){
        
        # if(floor(i/10) == i/10){
        #     print(i)
        # }
        
        dat_i = dat[i,]
        zz1 = dat_i$Z1
        zz2 = dat_i$Z2
        
        TT1 = T.min + rweibull(nn, shape = shape.T, scale = exp(-1/shape.T * (cbind(1, 1, zz1, zz2, zz1, zz2)%*%beta.T)))
        TT0 = T.min + rweibull(nn, shape = shape.T, scale = exp(-1/shape.T * (cbind(1, 0, zz1, zz2, -zz1, -zz2)%*%beta.T)))
        
        CATE.truth[i] = mean( nu(TT1) - nu(TT0) ) 
    }
    
    return(CATE.truth)
}

# plot(dat$Z2, CATE.truth)
# plot(dat$Z1, CATE.truth)



##!!!! There are probably bugs in the following code.

## function for xgboost

## Use another sample splitting for estimating the weights involved in nuisance estimation and the nuisance estimation
## the cross fitted truncAC_AIPW estimator

#' @param cov.names.binary.T vector of names for binary covariates if \texttt{model = "pCox"} is used. For other models, \texttt{c(cov.names, cov.names,binary)} is used as the covariates.
#' @param options.F a list of arguments in est_F() that are used for specific methods, including trim, formula.survPen, df, nfolds, s, alpha, lambda  
#' @param options.G a list of arguments in est_G() that are used for specific methods, including trim, formula.survPen, df, nfolds, s, alpha, lambda, trunc.weighting, where 'trunc.weighting' is a logical value indicating whether truncation weights are use to estimate G; if FALSE, G is estimated from uncensored subjects with IPCW weights 
#' @param options.PS a list of arguments in est_PS() that are used for specific methods, including trim, df. Currently PS is always estimated using truncating weighting.
#' @param options.Sd a list of arguments in est_Sd() that are used for specific methods, including trim, mtry, ntree, formula.survPen, df, nfolds, s, alpha, lambda  

cf2_CATE_truncAC_AIPW <- function(dat, CATE.alg = "xgb", K, nu, cov.names.CATE,
                                  X.name, Q.name, event.name, A.name, 
                             model.T, model.Q, model.A, model.D, 
                             cov.names.T, cov.names.Q, cov.names.A, cov.names.D, trim = 1e-7,
                             cov.names.binary.T = NULL, cov.names.binary.Q = NULL, cov.names.binary.A = NULL, cov.names.binary.D = NULL,
                             options.F = NULL, options.G = NULL, 
                             options.PS = NULL, options.Sd = NULL,
                             options.CATE.xgb = NULL){
    
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
    
    # Compute \V{\nu} and \V{1} using out of fold estimate of the nuisance parameters
    truncC_AIPW_weights = rep(NA, n)  # \V{1}
    y_tilde = rep(NA,n)  # [\V{\nu(T)\} - m(Z)] / [\V{1}(A-\pi(Z))]
    PS_weights = rep(NA,n)  # {A-\pi(Z)}^2
    
    ct = 0
    for(k in 1:K){
        dat.est = dat[folds$subsets[folds$which == k], ]
        dat.fit = dat[folds$subsets[folds$which != k], ]
        dat_A1.est = dat_A1[folds$subsets[folds$which == k], ]
        dat_A0.est = dat_A0[folds$subsets[folds$which == k], ]
        
        
        # split the fit data into 2 folds
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
        
        
        # Compute m(Z) from F and \pi
        nuuT.mx = matrix(rep(nu(u.T), nrow(dat.est)), nrow = nrow(dat.est), byrow = TRUE)
        tau.Tmax = max(u.T) + 1e-5
        mu_A1 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A1.mx))
        mu_A0 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A0.mx))
        m_Z = mu_A1 * PS + mu_A0 * (1-PS)
        
        
        # Compute \V{1}, y_tilde, and (A-\pi(Z))^2
        
        truncC_AIPW_result <- truncC_AIPW_transMean(dat.est, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                                                    X.name, Q.name, event.name, trim)
        
        nk = nrow(dat.est)
        truncC_AIPW_const1 = truncC_AIPW_result$truncC_const1
        truncC_AIPW_weights[(ct+1):(ct+nk)] = truncC_AIPW_const1
        
        truncC_AIPW_nuT = truncC_AIPW_result$truncC_nu
        A.est = dat.est[,A.name]
        y_tilde[(ct+1):(ct+nk)] = (truncC_AIPW_nuT/truncC_AIPW_const1 - m_Z) / (A.est-PS)
        
        PS_weights[(ct+1):(ct+nk)] = (A.est-PS)^2
        
        
        ct = ct + nk
    }
    
    
    if(CATE.alg == "xgb"){
        
        Z_CATE = as.matrix(dat[,cov.names.CATE, drop = FALSE])
        dtrain <- xgb.DMatrix(data = Z_CATE, label = y_tilde)
        
        custom_weights = truncC_AIPW_weights * PS_weights
        
        # Define the custom weighted squared loss objective function
        # preds: model's predictions
        # dtrain: xgb.DMatrix containing the training data and labels
        custom_weighted_squared_loss <- function(preds, dtrain) {
            labels <- getinfo(dtrain, "label")
            residuals <- preds - labels
            # Gradient of the weighted squared loss with respect to predictions
            grad <- 2 * residuals * custom_weights
            # Hessian (second derivative) of the loss with respect to predictions
            hess <- 2 * custom_weights
            
            return(list(grad = grad, hess = hess))
        }
        
        # setinfo(dtrain, "custom_weights",  truncC_AIPW_weights * PS_weights )  # This sets equal weight to each instance, replace with actual weights
        # Model parameters
        params <- list(booster = options.CATE.xgb$booster, objective = custom_weighted_squared_loss)
        
        # Training the model with the custom objective
        bst <- xgb.train(params = params, data = dtrain, nrounds = options.CATE.xgb$nrounds)
        
        # make predictions
        dtest <- xgb.DMatrix(data = Z_CATE)
        CATE.pred <- predict(bst, dtest)
        
    }else if(CATE.alg == "xgb.BW"){  ## bound weights to be positive
        
        Z_CATE = as.matrix(dat[,cov.names.CATE, drop = FALSE])
        dtrain <- xgb.DMatrix(data = Z_CATE, label = y_tilde)
        
        custom_weights = truncC_AIPW_weights * PS_weights
        
        # Define the custom weighted squared loss objective function
        # preds: model's predictions
        # dtrain: xgb.DMatrix containing the training data and labels
        custom_weighted_squared_loss <- function(preds, dtrain) {
            labels <- getinfo(dtrain, "label")
            residuals <- preds - labels
            # Gradient of the weighted squared loss with respect to predictions
            grad <- 2 * residuals * custom_weights
            # Hessian (second derivative) of the loss with respect to predictions
            hess <- 2 * custom_weights
            
            return(list(grad = grad, hess = hess))
        }
        
        # setinfo(dtrain, "custom_weights",  truncC_AIPW_weights * PS_weights )  # This sets equal weight to each instance, replace with actual weights
        # Model parameters
        params <- list(booster = options.CATE.xgb$booster, objective = custom_weighted_squared_loss)
        
        # Training the model with the custom objective
        bst <- xgb.train(params = params, data = dtrain, nrounds = options.CATE.xgb$nrounds)
        
        # make predictions
        dtest <- xgb.DMatrix(data = Z_CATE)
        CATE.pred <- predict(bst, dtest)
        
    }else{
        stop(paste("The algorithm '", CATE.alg, "' is not implemented in this function.", sep = ""))
    }
    
    
    return(list(CATE.pred = CATE.pred, Z.CATE = Z_CATE))
    
}




