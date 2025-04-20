## function for CATE estimation using R-learner without left truncation or right censoring
## In this case, X = T 
## The original name of this Rscript is `developement_rlearner.R`

# r-learner for continous outcome without left truncation or right censoring

# # For simu2 simulation scheme
# oracle_CATE_AIPW <- function(dat, 
#                              T.name, A.name,
#                              cov.names.T.oracle, T.min, beta.T, shape.T,
#                              cov.names.A.oracle, gamma.A,
#                              cov.names.CATE,
#                              booster = "gbtree",
#                              eta = 0.3,
#                              gamma = 0.5,
#                              max_depth = 6,
#                              subsample = 0.7,
#                              lambda = 0.5,
#                              alpha = 0,
#                              nfold = 10,
#                              nrounds = 100,
#                              metrics = "rmse",
#                              k_folds = NULL,
#                              # ntrees_max = 1000,
#                              # num_search_rounds = 10,
#                              # print_every_n = 100,
#                              # early_stopping_rounds = 10,
#                              nthread = NULL,
#                              verbose = FALSE,
#                              lambda_choice = c("lambda.min","lambda.1se")){
#     
#     # !! The following need to be manually modified, is not automated yet
#     dat_A1 = dat
#     dat_A1[,A.name] <- 1
#     dat_A1$AZ1 = dat$Z1
#     dat_A1$AZ2 = dat$Z2
#     
#     dat_A0 = dat
#     dat_A0[,A.name] <- 0
#     dat_A0$AZ1 = 0
#     dat_A0$AZ2 = 0
#     
#     n = nrow(dat)
#     
#     jumps.T = sort(dat[,T.name])
#     
#     u = c(min(jumps.T)-1e-10, jumps.T, max(jumps.T)+1e-10)
# 
#     tau.max = max(jumps.T)
#     
#     # Compute the truth of the nuisance parameters
#     Fuz_A1.mx = truth_F(dat_A1, u, 
#                         cov.names = cov.names.T.oracle, T.min, beta.T, shape.T)
#     Fuz_A0.mx = truth_F(dat_A0, u, 
#                         cov.names = cov.names.T.oracle, T.min, beta.T, shape.T)
#     
#     PS = truth_PS(dat, cov.names.A.oracle, gamma.A)
#     
#     # compute the mu_1 and mu_0 for the treatment augmentation
#     u.T = c(min(jumps.T)-1e-10, jumps.T, max(jumps.T)+1e-10)  # need the values of F in all jumps of T to compute mu()
#     
#     
#     # Compute m(Z) from F and \pi
#     nuuT.mx = matrix(rep(nu(u.T), nrow(dat)), nrow = nrow(dat), byrow = TRUE)
#     tau.Tmax = max(u.T) + 1e-5
#     mu_A1 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A1.mx))
#     mu_A0 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A0.mx))
#     m_Z = mu_A1 * PS + mu_A0 * (1-PS)
#     
#     nuT.vec = nu(dat[,T.name])
#     
#     A = dat[,A.name]
#     
#     yy = nuT.vec - m_Z
#     A_tilde = A - PS
#     y_tilde = yy/A_tilde  # pseudo_outcome
#     
#     weights = A_tilde^2
#     
#     Z_CATE = as.matrix(dat[,cov.names.CATE, drop = FALSE])  # `drop = FALSE` will make sure that the output is a matrix even when there is only one covariate
#     
#     
#     ## estimate the CATE using xgbost
#     # Using cvboot() or rboost from the rlearner package
#     
#     # tau_fit.rboost = rboost(x = Z_CATE, w = A, y = nuT.vec)
#     # 
#     # 
#     # tau_fit = cvboost(Z_CATE,
#     #                   pseudo_outcome,
#     #                   objective = "reg:squarederror",
#     #                   weights = weights,
#     #                   k_folds = k_folds,
#     #                   ntrees_max = ntrees_max,
#     #                   num_search_rounds = num_search_rounds,
#     #                   print_every_n = print_every_n,
#     #                   early_stopping_rounds = early_stopping_rounds,
#     #                   nthread = nthread,
#     #                   verbose = verbose)
#     
#     #!!!! The above commented code is not working
# 
#     
#     ### Using the xgb.cv() function
#     
#     # Create DMatrix including the weights
#     dtrain <- xgb.DMatrix(data = Z_CATE, label = y_tilde, weight = weights)
#     
#     ## Cross-validation parameters
#     # Previous setting
#     params <- list(objective = "reg:squarederror", 
#                    booster = booster,  # "gblinear", "gbtree", 
#                    eta = eta,
#                    gamma = gamma,
#                    max_depth = max_depth,
#                    subsample = subsample,
#                    lambda = lambda,
#                    alpha = alpha)
#     
#     # Perform cross-validation
#     cv_results <- xgb.cv(params = params, 
#                          data = dtrain, 
#                          nrounds = nrounds, 
#                          nfold = nfold, 
#                          metrics = metrics, 
#                          showsd = TRUE, 
#                          stratified = FALSE, 
#                          print_every_n = 10, 
#                          early_stopping_rounds = 10, 
#                          maximize = FALSE)
#     
#     # Extract the Best Number of Rounds
#     best_nrounds <- cv_results$best_iteration
#     best_nrounds
#     
#     e_log = cv_results$evaluation_log
#     plot(e_log$iter, e_log$test_rmse_mean)
#     
#     # Train the Final Model
#     final_model <- xgb.train(params = list(objective = "reg:squarederror", booster = "gbtree"),
#                              data = dtrain,
#                              nrounds = best_nrounds)
#     
#     
#     # make predictions
#     dtest <- xgb.DMatrix(data = Z_CATE)
#     CATE_est.weightSquare_cv <- predict(final_model, dtest)
#     summary(CATE_est.weightSquare_cv)
#     summary(CATE_est.weightSquare_cv - CATE.truth.full)
#     
#     
#     
#     ### Uisng customized loss function - square loss
#     dtrain <- xgb.DMatrix(data = Z_CATE, label = y_tilde)
#     
#     custom_weights = weights
#     
#     # Define the custom weighted squared loss objective function
#     # preds: model's predictions
#     # dtrain: xgb.DMatrix containing the training data and labels
#     custom_weighted_squared_loss <- function(preds, dtrain) {
#         labels <- getinfo(dtrain, "label")
#         residuals <- preds - labels
#         # Gradient of the weighted squared loss with respect to predictions
#         grad <- 2 * residuals * custom_weights
#         # Hessian (second derivative) of the loss with respect to predictions
#         hess <- 2 * custom_weights
#         
#         return(list(grad = grad, hess = hess))
#     }
#     
#     # setinfo(dtrain, "custom_weights",  truncC_AIPW_weights * PS_weights )  # This sets equal weight to each instance, replace with actual weights
#     # Model parameters
#     params <- list(booster = options.CATE.xgb$booster, objective = custom_weighted_squared_loss)
#     
#     # Training the model with the custom objective
#     # bst <- xgb.train(params = params, data = dtrain, nrounds = options.CATE.xgb$nrounds)   # using the default
#     bst <- xgb.train(params = params, data = dtrain, nrounds =  best_nrounds)
#     
#     # make predictions
#     dtest <- xgb.DMatrix(data = Z_CATE)
#     CATE_est.cutomizedLoss <- predict(bst, dtest)
#     
#     
#     
#     
#     
#     
#     
#     
#     ### Try glmnet with CV
#     
#     
#     
#     
#     return(list(CATE_est.cutomizedLoss = CATE_est.cutomizedLoss,
#                 CATE_est.weightSquare_cv = CATE_est.weightSquare_cv))
#     
#     
# }




#' @param lambda_choice for fitting using lasso

oracle_HTE2_rlearner <- function(dat.full, 
                             T.name, A.name,
                             cov.names.T.oracle, beta.T, sigma_eps,
                             cov.names.A.oracle, gamma.A,
                             cov.names.CATE,
                             booster = "gbtree",
                             eta = 0.3,
                             gamma = 0.5,
                             max_depth = 6,
                             subsample = 0.7,
                             lambda = 0.5,
                             alpha = 0,
                             nfold = 10,
                             nrounds = 100,
                             metrics = "rmse",
                             k_folds=NULL,
                             objective= "reg:squarederror",
                             ntrees_max=1000,
                             num_search_rounds=10,
                             print_every_n=100,
                             early_stopping_rounds=10,
                             nthread=NULL,
                             verbose=FALSE,
                             lambda_choice = c("lambda.min","lambda.1se")){
    
    dat = dat.full
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
    
    n = nrow(dat)
    
    jumps.T = sort(dat[,T.name])
    
    u = c(min(jumps.T)-1e-10, jumps.T, max(jumps.T)+1e-10)
    
    tau.max = max(jumps.T)
    
    # Compute the truth of the nuisance parameters
    Fuz_A1.mx = truth_F.AFT(dat_A1, u, cov.names = cov.names.T.oracle, beta.T, 
                            err.dist = "lognormal", err.sigma = sigma_eps)
    Fuz_A0.mx = truth_F.AFT(dat_A0, u, cov.names = cov.names.T.oracle, beta.T, 
                            err.dist = "lognormal", err.sigma = sigma_eps)
    
    PS = truth_PS(dat, cov.names.A.oracle, gamma.A)
    
    # compute the mu_1 and mu_0 for the treatment augmentation
    u.T = c(min(jumps.T)-1e-10, jumps.T, max(jumps.T)+1e-10)  # need the values of F in all jumps of T to compute mu()
    
    
    # Compute m(Z) from F and \pi
    nuuT.mx = matrix(rep(nu(u.T), nrow(dat)), nrow = nrow(dat), byrow = TRUE)
    tau.Tmax = max(u.T) + 1e-5
    mu_A1 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A1.mx))
    mu_A0 = as.vector(int_fmx_dF(tau.Tmax, nuuT.mx, Fuz_A0.mx))
    m_Z = mu_A1 * PS + mu_A0 * (1-PS)
    
    nuT.vec = nu(dat[,T.name])
    
    A = dat[,A.name]
    
    yy = nuT.vec - m_Z
    A_tilde = A - PS
    y_tilde = yy/A_tilde  # pseudo_outcome
    
    weights = A_tilde^2
    
    Z_CATE = as.matrix(dat[,cov.names.CATE, drop = FALSE])  # `drop = FALSE` will make sure that the output is a matrix even when there is only one covariate
    
    
    ##### estimate the CATE using xgbost - Using rlearner::cvboot() to select the best tuning parameters
    
    ### Using the xgb.cv() function
    
    # Create DMatrix including the weights
    dtrain <- xgb.DMatrix(data = Z_CATE, label = y_tilde, weight = weights)
    
    ## Cross-validation parameters
    # Previous setting
    params <- list(objective = "reg:squarederror", 
                   booster = booster,  # "gblinear", "gbtree", 
                   eta = eta,
                   gamma = gamma,
                   max_depth = max_depth,
                   subsample = subsample,
                   lambda = lambda,
                   alpha = alpha)

    # Perform cross-validation
    cv_results <- xgb.cv(params = params, 
                         data = dtrain, 
                         nrounds = nrounds, 
                         nfold = nfold, 
                         metrics = metrics, 
                         showsd = TRUE, 
                         stratified = FALSE, 
                         print_every_n = 10, 
                         early_stopping_rounds = 10, 
                         maximize = FALSE)
    
    
    # Extract the Best Number of Rounds
    best_nrounds <- cv_results$best_iteration
    # best_nrounds
    
    e_log = cv_results$evaluation_log
    # plot(e_log$iter, e_log$test_rmse_mean)
    
    # Train the Final Model
    final_model <- xgb.train(params = list(objective = "reg:squarederror", booster = "gbtree"),
                             data = dtrain,
                             nrounds = best_nrounds)
    
    
    # make predictions
    dtest <- xgb.DMatrix(data = Z_CATE)
    CATE_est.weightSquare_cv <- predict(final_model, dtest)
    # summary(CATE_est.weightSquare_cv)
    # summary(CATE_est.weightSquare_cv - CATE.truth.full)
    
    
    
    
    
    
    ### Uisng customized loss function - square loss
    dtrain <- xgb.DMatrix(data = Z_CATE, label = y_tilde)
    
    custom_weights = weights
    
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
    # bst <- xgb.train(params = params, data = dtrain, nrounds = options.CATE.xgb$nrounds)   # using the default
    bst <- xgb.train(params = params, data = dtrain, nrounds =  best_nrounds)
    
    # make predictions
    dtest <- xgb.DMatrix(data = Z_CATE)
    CATE_est.cutomizedLoss <- predict(bst, dtest)
    
    
    
    ### xgboost with tuning using rlearner::cvboost()
    
    fit = cvboost(Z_CATE,
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
    
    
    
    
    
    
    ##### Fit the true linear model
    dat$y_tilde = y_tilde
    lm_fit = lm(y_tilde ~ Z1, data = dat, weights = weights)
    CATE_est.lm = predict(lm_fit, newdata = dat)    
    coef_lm = coef(lm_fit)
    
    
    ##### Fit lasso
    ###!!! The following design matrix is only for testing purposes. In  practice need to replace by natural spline basis
    Z_CATE = as.matrix(cbind(Z1 = dat$Z1, Z2 = dat$Z2, Z2sqrt = dat$Z2sqrt))
    standardization = caret::preProcess(Z_CATE, method=c("center", "scale")) # get the standardization params
    x_scl = predict(standardization, Z_CATE)							 # standardize the input
    x_scl = x_scl[,!is.na(colSums(x_scl)), drop = FALSE]
    
    s_scl_tilde = cbind(as.numeric(A-PS) * cbind(1, x_scl))
    x_scl_pred = cbind(1, x_scl)
    
    tau_fit = glmnet::cv.glmnet(s_scl_tilde, yy, alpha = 1, standardize = FALSE)
    tau_beta = as.vector(t(coef(tau_fit, s = lambda_choice)[-1]))
    
    CATE_est.glmnet = x_scl_pred %*% tau_beta
    
    
    
    return(list(CATE_est.cutomizedLoss = CATE_est.cutomizedLoss,
                CATE_est.weightSquare_cv = CATE_est.weightSquare_cv,
                best_nrounds = best_nrounds,
                e_log = e_log,
                CATE_est.lm = CATE_est.lm, coef_lm = coef_lm,
                CATE_est.glmnet =  CATE_est.glmnet, 
                CATE_est.cvboost = CATE_est.cvboost,
                y_tilde = y_tilde,
                weights = weights,
                Z_CATE = Z_CATE))
}





