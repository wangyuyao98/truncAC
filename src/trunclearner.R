## The function for implementing the truncrlearner with K-fold cross-fitting
## For estimating the nuisance parameters, we further split the data into 2 folds when one nuisance parameter estimate depend on the other nuisance estimate. 
## We further switch the role of the two folds of data and take the average of the 2 nuisance estimates.


#' @param K the number of folds for K-fold cross-fitting for computing the loss function
#' @param est_approach_G The Approach of estimating the nuisance parameter G. There are two options: "truncIPW.F" uses truncation weights $1/{1-F(Q|A,Z)}$; "IPCW" uses IPCW weights $\Detla/S_D(X-Q|A,Z)$.
#' @param est_approach_PS The Approach of estimating the nuisance parameter propensity score. There are two options: "truncIPW.F" uses truncation weights $1/\{1-F(Q|A,Z)\}$; "truncIPCW" uses IPW weights $\Detla/\{G(X|A,Z)S_D(X-Q|A,Z)\}$.
trunclearner <- function(dat, nu, cov.names.CATE, cov.names.CATE.binary = NULL, K,
                          X.name, Q.name, event.name, A.name, trim = 1e-7,
                          model.T = NULL, model.Q = NULL, model.A = NULL, model.D = NULL, 
                          cov.names.T = NULL, cov.names.Q = NULL, cov.names.A = NULL, cov.names.D = NULL, 
                          cov.names.binary.T = NULL, cov.names.binary.Q = NULL, cov.names.binary.A = NULL, cov.names.binary.D = NULL,
                          options.F = NULL, options.G = NULL, 
                          options.PS = NULL, options.Sd = NULL,
                          est_approach_G = "truncIPW.F", est_approach_PS = "truncIPW.F", 
                          Fuz.mx = NULL, Fuz_A1.mx = NULL, Fuz_A0.mx = NULL,
                          Gvz.mx = NULL, Sdz.mx = NULL, PS = NULL,
                          df = NULL, nfolds = 10, # alpha = NULL, lambda = NULL,
                          metrics = "rmse",
                          booster = "gbtree",
                          k_folds=NULL,
                          objective= "reg:squarederror",
                          ntrees_max=500,
                          num_search_rounds=200,
                          print_every_n=100,
                          early_stopping_rounds=10,
                          nthread=NULL,
                          verbose=FALSE, 
                          simulation = TRUE,
                          trim1 = 0.05,
                          trim2 = 0,
                          save_cvfit = FALSE){
    # Remove the df and alpha option in the above imputs of the function
    # lambda also seems to be un-used
    
    if(simulation){
        # !! The following need to be manually modified, is not automated yet
        dat_A1 = dat
        dat_A1[,A.name] <- 1
        dat_A1$AZ1 = dat$Z1
        dat_A1$AZ2 = dat$Z2
        dat_A1$Z2sqrt = sqrt(abs(dat$Z2))
        dat_A1$A_sinpiZ1 = sin(pi*dat$Z1)
        dat_A1$A_Z2sqrt = sqrt(abs(dat$Z2))
        
        dat_A0 = dat
        dat_A0[,A.name] <- 0
        dat_A0$AZ1 = 0
        dat_A0$AZ2 = 0
        dat_A0$Z2sqrt = sqrt(abs(dat$Z2))
        dat_A0$A_sinpiZ1 = 0
        dat_A1$A_Z2sqrt = 0
        
    }else{
        dat_A1 = dat
        dat_A1[,A.name] <- 1
        dat_A0 = dat
        dat_A0[,A.name] <- 0
    }
    
    
    if(is.null(Gvz.mx)){
        names = c(X.name, Q.name, event.name, A.name, cov.names.T, cov.names.Q, cov.names.A, cov.names.D)
        if(sum(names == "delta.1")){
            stop("The names of the variables cannot be 'delta.1'. It is used in the middle of the computation.")
        }
        dat$delta.1 = rep(1, nrow(dat))
        event.name.Q = "delta.1"
    }
    
    
    n = nrow(dat)
    event.name.T = event.name
    
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
    # y_R = rep(NA,n)  # [\V{\nu(T)\} - m(Z)] / [\V{1}(A-\pi(Z))]
    # PS_weights = rep(NA,n)  # {A-\pi(Z)}^2
    
    # Store some intermediate results that are useful for constructing IPW and naive estimators
    mZ_hat = rep(NA, n)
    mu0_hat = rep(NA, n)
    mu1_hat = rep(NA, n)
    GXZ_hat = rep(NA, n)
    SdXZ_hat = rep(NA, n)
    PS_hat = rep(NA, n)
    
    t1 <- Sys.time()  # record time
    
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
            # # Look into the estimation line by line --------------------------------
            # dat.fit = dat.fit
            # dat.est = dat.est
            # time.eval = d
            # model = model.D
            # time.name = X.name
            # Q.name = Q.name
            # event.name = event.name.T
            # cov.names = cov.names.D
            # cov.names.binary = cov.names.binary.D
            # trim = options.Sd$trim 
            # mtry = options.Sd$mtry
            # ntree = options.Sd$ntree 
            # formula.survPen = options.Sd$formula.survPen 
            # nfolds = options.Sd$nfolds
            # s = options.Sd$s
            # alpha = options.Sd$alpha
            # lambda = options.Sd$lambda
            # df = options.Sd$df
            # # End ----------------------------------------------------------------
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
                
                # # Look into the estimate - line-by-line -----------------------------------
                # dat.fit = dat.fit_1
                # dat.est = dat.est
                # time.eval = v
                # model = model.Q
                # time.name = X.name
                # event.name = event.name.Q
                # cov.names = cov.names.Q
                # cov.names.binary = cov.names.binary.Q
                # weights = w_IPCW_1
                # tau = tau.max
                # trim = options.G$trim
                # formula.survPen = options.G$formula.survPen
                # nfolds = options.G$nfolds
                # s = options.G$s
                # alpha = options.G$alpha
                # lambda = options.G$lambda
                # df = options.G$df
                # trunc = TRUE  # trunc = TRUE for IPCW approach
                # # End ----------------------------------------------------------------------
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
                    Gx.vec.fit_1 = diag(CDF_eval_mx_cpp(X_1, Gvz.mx[id.fit,][id.fit_1,], v))  # For uncensored subjects
                }
                
                w_trunc.fit_1 = 1/pmax(Gx.vec.fit_1, trim)
                
                
                # Estimate PS using uncensored subjects with weights 1/{G(X|A,Z)S_D(X-Q|A,Z)}
                # # Look into - line-by-line ---------------------------------------------------
                # dat.fit = dat.fit_1
                # dat.est = dat.est
                # model = model.A
                # A.name = A.name
                # cov.names = cov.names.A
                # cov.names.binary = cov.names.binary.A 
                # weights = w_trunc.fit_1 * w_IPCW.fit_1
                # trim = options.PS$trim
                # df = options.PS$df
                # ntree = options.PS$ntree
                # # End ----------------------------------------------------------------------
                
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
        m_Z.est = mu_A1.est * PS.est + mu_A0.est * (1-PS.est)
        
        
        # Compute \V{1}, y_R, and (A-\pi(Z))^2
        truncC_AIPW_result <- truncC_AIPW_transMean(dat.est, nu, Fuz.mx.est, Gvz.mx.est, Sdz.mx.est,
                                                    X.name, Q.name, event.name, trim)
        
        # truncCAIPW_const1 = truncC_AIPW_result$truncC_const1
        truncC_AIPW_1[id.est] = truncC_AIPW_result$truncC_const1
        truncC_AIPW_nuT[id.est] = truncC_AIPW_result$truncC_nu
        
        
        ### Record some intermediate quantities
        mZ_hat[id.est] = m_Z.est
        mu1_hat[id.est] = mu_A1.est
        mu0_hat[id.est] = mu_A0.est
        PS_hat[id.est] = PS.est
        GXZ_hat[id.est] = diag(CDF_eval_mx_cpp(dat.est[,X.name], Gvz.mx.est, v))
        SdXZ_hat[id.est] = diag(CDF_eval_mx_cpp(dat.est[,X.name] - dat.est[,Q.name], Sdz.mx.est, d))
        
    }
    
    t2 <- Sys.time()  # record time
    
    print(paste("Finished computation for cf part. Time elapsed (minutes):", 
                round(as.numeric(difftime(t2, t1, units = "mins")), 2)))
    

    A = dat[,A.name]
    nuT_tilde = truncC_AIPW_nuT/bound_away_zero(truncC_AIPW_1, trim1)
    
    
    
    
    ############### Estimate CATE ################
    
    cov.names.CATE = c(cov.names.CATE, cov.names.CATE.binary)
    Z_CATE = as.matrix(dat[,cov.names.CATE, drop = FALSE])  # `drop = FALSE` will make sure that the output is a matrix even when there is only one covariate

    # truncR-learner with xgboost
    PS_weights = (bound_away_zero(A-PS_hat, 0))^2  # used to be trim1
    weights_R = bound_away_zero(bound_away_zero(truncC_AIPW_1, 0)*PS_weights, trim2)  # used to be trim1
    y_R = (nuT_tilde - mZ_hat) / bound_away_zero(A-PS_hat, trim1)
    # # Begin test step-by-step for cvboost_wsq() ---------------------------------------------------
    # x = Z_CATE
    # y = y_R
    # weights = weights_R
    # k_folds = k_folds
    # ntrees_max = ntrees_max
    # num_search_rounds = num_search_rounds
    # print_every_n = print_every_n
    # nthread  = nthread
    # verbose = verbose
    # # End test ------------------------------------------------------------------
    cvfit_wsq_R  = cvboost_wsq(Z_CATE, y_R, weights = weights_R,
                             k_folds = k_folds,
                             ntrees_max = ntrees_max, 
                             num_search_rounds = num_search_rounds,
                             print_every_n = print_every_n, 
                             nthread  = nthread,
                             verbose = verbose)
    CATE_est_R = predict(cvfit_wsq_R, newx = Z_CATE)
    
    t3 <- Sys.time()  # record time
    
    print(paste("Finished ltrcR-learner tuning and final prediction. Time elapsed (minutes):", 
                round(as.numeric(difftime(t3, t2, units = "mins")), 2)))
    
    
    # truncDR-learner with xgboost
    mu_A.vec = A*mu1_hat + (1-A)*mu0_hat
    y_DR = (A-PS_hat)/(pmax(PS_hat, trim)*pmax(1-PS_hat, trim)) * (nuT_tilde - mu_A.vec) + mu1_hat - mu0_hat
    weights_DR = bound_away_zero(truncC_AIPW_1, 0)  # used to be trim1
    cvfit_wsq_DR  = cvboost_wsq(Z_CATE, y_DR, weights = weights_DR, 
                               k_folds = k_folds,
                               ntrees_max = ntrees_max, 
                               num_search_rounds = num_search_rounds,
                               print_every_n = print_every_n, 
                               nthread  = nthread,
                               verbose = verbose)
    CATE_est_DR = predict(cvfit_wsq_DR, newx = Z_CATE)
    
    t4 <- Sys.time()  # record time
    
    print(paste("Finished ltrcDR-learner tuning and final prediction. Time elapsed (minutes):", 
                round(as.numeric(difftime(t4, t3, units = "mins")), 2)))
    
    
    # S-learner with xgboost
    Delta = dat[, event.name]
    id.Delta1 = (Delta == 1)
    Z_CATE_S = Z_CATE[id.Delta1, , drop = FALSE]
    weights_S = 1/(pmax(GXZ_hat, trim)*pmax(SdXZ_hat, trim))[id.Delta1]
    y_S =  (mu1_hat - mu0_hat)[id.Delta1]
    # ## Begin debug -- step-by-step ------------------------------------------------
    # x = Z_CATE_S
    # y = y_S
    # objective = "reg:squarederror"
    # weights = weights_S
    # k_folds = k_folds
    # ntrees_max = ntrees_max
    # num_search_rounds = num_search_rounds
    # print_every_n = print_every_n
    # early_stopping_rounds = early_stopping_rounds
    # nthread = nthread
    # verbose = verbose
    # ## End debug ------------------------------------------------------------------
    cvfit_S = cvboost2(Z_CATE_S,
                       y_S,
                       objective = "reg:squarederror",
                       weights = weights_S,
                       k_folds = k_folds,
                       ntrees_max = ntrees_max,
                       num_search_rounds = num_search_rounds,
                       print_every_n = print_every_n,
                       early_stopping_rounds = early_stopping_rounds,
                       nthread = nthread,
                       verbose = verbose)
    CATE_est_S = predict(cvfit_S, newx = Z_CATE)
    
    t5 <- Sys.time()  # record time
    
    print(paste("Finished IPW.S-learner tuning and final prediction. Time elapsed (minutes):", 
                round(as.numeric(difftime(t5, t4, units = "mins")), 2)))
    
    print(paste("Total time elapsed (minutes):", 
                round(as.numeric(difftime(t5, t1, units = "mins")), 2)))
    
    
    # ## Visualize the estimated CATE when there is only one covariate and save the results -----------------
    # pdf(paste("HAAS_analysis/plots/CATE/nonparametric/CATE_surv90_pCox_", cov.names.CATE,
    #           "_", options.F$s, ".pdf", sep = ""), 
    #     width = 6, height = 6)
    # plot(dat[,cov.names.CATE], CATE_est_DR, col = 2, 
    #      xlab = cov.names.CATE, ylab = "CATE")
    # points(dat[,cov.names.CATE], CATE_est_R, col = 3)
    # points(dat[,cov.names.CATE], CATE_est_S, col = 4)
    # legend("topright", legend = c("ltrcR", "ltrcDR", "IPW.S"), lty = 1, col = c(3,2,4), bty = "n")
    # dev.off()
    # 
    # # Smoothing using lowess
    # f = 2/3
    # pdf(paste("HAAS_analysis/plots/CATE/nonparametric/CATE_surv90_lowess_pCox_", cov.names.CATE,
    #           "_", options.F$s, ".pdf", sep = ""), 
    #     width = 6, height = 6)
    # # Create a smoothing curve for each estimate using lowess()
    # plot(dat[,cov.names.CATE], CATE_est_DR, col = 2, 
    #      xlab = cov.names.CATE, ylab = "CATE", type = "n")
    # # Add smoothing curves using lowess()
    # lines(lowess(dat[,cov.names.CATE], CATE_est_R, f = f), col = 3, lwd = 2, lty = 1)
    # lines(lowess(dat[,cov.names.CATE], CATE_est_DR, f = f), col = 2, lwd = 2, lty = 2)
    # lines(lowess(dat[,cov.names.CATE], CATE_est_S, f = f), col = 4, lwd = 3, lty = 3)
    # # Add legend
    # legend("topright", legend = c("ltrcR", "ltrcDR", "IPW.S"), 
    #        lty = 1:3, col = c(3, 2, 4), bty = "n", lwd = c(2,2,3))
    # dev.off()
    # 
    # results = list(est_DR = CATE_est_DR,
    #                est_R = CATE_est_R,
    #                est_S = CATE_est_S,
    #                truncC_AIPW_1 = truncC_AIPW_1,
    #                truncC_AIPW_nuT = truncC_AIPW_nuT,
    #                mZ_hat = mZ_hat,
    #                mu1_hat = mu1_hat,
    #                mu0_hat = mu0_hat,
    #                PS_hat = PS_hat,
    #                GXZ_hat = GXZ_hat,
    #                SdXZ_hat = SdXZ_hat)
    # Z_CATE = dat[,cov.names.CATE]
    # 
    # save(results, cov.names.CATE, Z_CATE,
    #      CATE_est_R, CATE_est_DR, CATE_est_S,
    #      truncC_AIPW_1, truncC_AIPW_nuT,
    #      y_R, weights_R, y_DR, weights_DR, y_S, weights_S,
    #      mZ_hat, mu1_hat, mu0_hat, PS_hat, GXZ_hat, SdXZ_hat,
    #      file = paste("HAAS_analysis/plots/CATE/nonparametric/surv90_Cox_3cov.rda", sep = "")
    #      # file = paste("HAAS_analysis/plots/CATE/nonparametric/surv90_Cox_", cov.names.CATE, ".rda", sep = "")
    #      )
    # ## End visualize and save -----------------------------------------------------------------
    # 
    
    if(save_cvfit){
        result_list = list(est_DR = CATE_est_DR,
                    est_R = CATE_est_R,
                    est_S = CATE_est_S,
                    cvfit_wsq_R = cvfit_wsq_R,
                    cvfit_wsq_DR = cvfit_wsq_DR,
                    cvfit_S = cvfit_S,
                    truncC_AIPW_1 = truncC_AIPW_1,
                    truncC_AIPW_nuT = truncC_AIPW_nuT,
                    mZ_hat = mZ_hat,
                    mu1_hat = mu1_hat,
                    mu0_hat = mu0_hat,
                    PS_hat = PS_hat,
                    GXZ_hat = GXZ_hat,
                    SdXZ_hat = SdXZ_hat)
    }else{
        result_list = list(est_DR = CATE_est_DR,
                    est_R = CATE_est_R,
                    est_S = CATE_est_S,
                    truncC_AIPW_1 = truncC_AIPW_1,
                    truncC_AIPW_nuT = truncC_AIPW_nuT,
                    mZ_hat = mZ_hat,
                    mu1_hat = mu1_hat,
                    mu0_hat = mu0_hat,
                    PS_hat = PS_hat,
                    GXZ_hat = GXZ_hat,
                    SdXZ_hat = SdXZ_hat)
    }
    
    return(result_list)
}






# x: a vector
# return max(x_i, trim) if x_i>0, and min(x_i, -trim) if x_i<0
bound_away_zero <- function(x, trim){
    sign_x = sign(x)
    sign_x[x == 0] <- 1
    y = sign_x * pmax(abs(x), trim)
    
    return(y)
}
