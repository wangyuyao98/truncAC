## 

truncrlearner.oracle_HTE3 <- function(dat, nu, cov.names.CATE, cov.names.CATE.binary = NULL,
                                      X.name, Q.name, event.name, A.name, 
                                      cov.names.T.oracle, beta.T, 
                                      cov.names.C.oracle, D.min, beta.C, shape.C,
                                      cov.names.Q.oracle, Q.min, Q.max, tau, beta.Q,
                                      cov.names.A.oracle, gamma.A,
                                      err.dist = c("weibull", "lognormal"),
                                      err.sigma = NULL, err.shape = NULL, err.scale = NULL,
                                      df = 3, nfolds = 10, alpha = 1, lambda = NULL,
                                      metrics = "rmse",
                                      booster = "gbtree",
                                      k_folds=NULL,
                                      objective= "reg:squarederror",
                                      ntrees_max=1000,
                                      num_search_rounds=100,
                                      print_every_n=100,
                                      early_stopping_rounds=10,
                                      nthread=NULL,
                                      verbose=FALSE,
                                      lambda_choice = c("lambda.min","lambda.1se")){
    
    
    cov.names.CATE = c(cov.names.CATE, cov.names.CATE.binary)
    Z_CATE = as.matrix(dat[,cov.names.CATE, drop = FALSE])  # `drop = FALSE` will make sure that the output is a matrix even when there is only one covariate
    
    
    # !! The following need to be manually modified, is not automated yet
    dat_A1 = dat
    dat_A1[,A.name] <- 1
    dat_A1$AZ1 = dat$Z1
    dat_A1$AZ2 = dat$Z2
    dat_A1$Z2sqrt = sqrt(abs(dat$Z2))
    
    dat_A0 = dat
    dat_A0[,A.name] <- 0
    dat_A0$AZ1 = 0
    dat_A0$AZ2 = 0
    dat_A0$Z2sqrt = sqrt(abs(dat$Z2))
    
    # names = c(X.name, Q.name, event.name, A.name, cov.names.T.oracle, cov.names.Q.oracle, cov.names.A.oracle, cov.names.C.oracle)
    # if(sum(names == "delta.1")){
    #     stop("The names of the variables cannot be 'delta.1'. It is used in the middle of the computation.")
    # }
    # 
    # dat$delta.1 = rep(1, nrow(dat))
    # event.name.Q = "delta.1"
    # event.name.T = event.name
    
    n = nrow(dat)
    
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
    
    
    ### Compute the truth of the nuisance parameters
    # cov.names.T.oracle = c("A","Z1","Z2sqrt")
    Fuz.mx = truth_F.AFT(dat, u, cov.names = cov.names.T.oracle, beta.T, 
                         err.dist = err.dist, err.sigma = err.sigma, 
                         err.shape = err.shape, err.scale = err.scale)
    Fuz_A1.mx = truth_F.AFT(dat_A1, u, cov.names = cov.names.T.oracle, beta.T, 
                            err.dist = err.dist, err.sigma = err.sigma, 
                            err.shape = err.shape, err.scale = err.scale)
    Fuz_A0.mx = truth_F.AFT(dat_A0, u, cov.names = cov.names.T.oracle, beta.T, 
                            err.dist = err.dist, err.sigma = err.sigma, 
                            err.shape = err.shape, err.scale = err.scale)
    
    # cov.names.C.oracle = c("Q","A","Z1","Z2","AZ1","AZ2")
    Sdz.mx = truth_Sd(dat, d, 
                      cov.names = cov.names.C.oracle, D.min, beta.C, shape.C)
    
    # cov.names.Q.oracle = c("A","Z1","Z2","AZ1","AZ2")
    Gvz.mx = truth_G(dat, v, 
                     cov.names = cov.names.Q.oracle, Q.min, Q.max, tau, beta.Q)
    
    # cov.names.A.oracle = c("Z1","Z2")
    PS = truth_PS(dat, cov.names.A.oracle, gamma.A)
    
    
    # compute the mu_1 and mu_0 for the treatment augmentation
    u.T = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)  # need the values of F in all jumps of T to compute mu()
    
    
    # Compute m(Z) from F and \pi
    nuuT.mx = matrix(rep(nu(u.T), nrow(dat)), nrow = nrow(dat), byrow = TRUE)
    tau.Tmax = max(u.T) + 1e-5
    mu_A1 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A1.mx))
    mu_A0 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A0.mx))
    m_Z = mu_A1 * PS + mu_A0 * (1-PS)
    
    
    # Compute \V{1}, y_tilde, and (A-\pi(Z))^2
    truncC_AIPW_result <- truncC_AIPW_transMean(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                                                X.name, Q.name, event.name, trim)
    
    truncC_AIPW_weights= truncC_AIPW_result$truncC_const1
    truncC_AIPW_nuT = truncC_AIPW_result$truncC_nu
    A = dat[,A.name]
    
    y_tilde = (truncC_AIPW_nuT/bound_away_zero(truncC_AIPW_weights, 1e-7) - m_Z) / (A-PS)
    PS_weights = (A-PS)^2
    truncPS_weights = truncC_AIPW_weights * (A-PS)^2
    
    
    ############### Estimate CATE ################
    
    ### lm() with the correctly specified model for tau() + bounded "weights"
    weights = pmax(truncPS_weights, 1e-7)
    
    dat$y_tilde = y_tilde
    lm_fit = lm(y_tilde ~ Z1, data = dat, weights = weights)
    CATE_est.lm = predict(lm_fit, newdata = dat)    
    coef_lm = coef(lm_fit)
    
    
    ### xgboost with cvboost2() to do tuning parameter selection + bounded "weights"
    fit = cvboost2(Z_CATE,
                  y_tilde,
                  objective = "reg:squarederror",
                  weights = weights,
                  k_folds = k_folds,
                  ntrees_max = ntrees_max,
                  num_search_rounds = num_search_rounds,
                  print_every_n = print_every_n,
                  early_stopping_rounds = early_stopping_rounds,
                  nthread = nthread,
                  verbose = verbose)
    X_test = Z_CATE
    CATE_est.cvboost = predict(fit, newx = X_test)
    
    
    
    return(list(CATE_est.lm = CATE_est.lm, coef_lm = coef_lm,
                CATE_est.cvboost = CATE_est.cvboost,
                y_tilde = y_tilde, truncPS_weights = truncPS_weights,
                PS_weights = PS_weights, truncC_AIPW_result = truncC_AIPW_result,
                Fuz.mx = Fuz.mx, Fuz_A1.mx = Fuz_A1.mx, Fuz_A0.mx = Fuz_A0.mx))
}
