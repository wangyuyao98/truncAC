### function that output the conditional survival probabilities for D at given time points ----------------
#' @title Estimate the Conditional survival probabilities of the residual censoring time given Covariates
#' @description Estimate the conditional csurvival probabilties of the event time given covariates evaluated at given time points. The options implemented in this function are: Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen'.
#' @param dat.fit data frame that is used to fit the model for the full data conditional distribution of the event time given the covariates.
#' @param dat.est data frame that contains the subjects for which the estimated conditional CDF is computed.
#' @param time.eval vector of time points at which the conditional CDF is evaluated.
#' @param model method used to estimate the conditional CDF. The options available are "Cox" and "survPen", corresponding to Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen', respectively.
#' @param time.name name of the censored event time variable.
#' @param Q.name name of the left truncation time variable.
#' @param event.name name of the event indicator.
#' @param cov.names vector of the names of covariates.
#' @param trim constant for bounding the estimated conditional CDF from 1.
#' @param formula.survPen the formula when applying the hazard model with penalized splines implemented in \code{survPen::survPen}.
#' @param nfolds number of folds for cross-validation
#' @param OOF a logical value indicating whether out of fold prediction will be used when dat.fit = dat.est
#' @param nfolds.OOF number of folds for doing out-of-fold prediction
#' @return \code{Sd_est()} returns a matrix of the estimated conditional conditional survival functions for subjects in `\code{data.est}' evaluated at the time points in the vector `\code{time.eval}'. Each row corresponds to a subject and each column corresponds to a time point. The column names of the matrix are the times in `\code{time.eval}'.
#' @export
#' @import survival stats survPen
#' @seealso  \code{\link{G_est}}, \code{\link{F_est}}
#' @examples
#' data("simu_truncAC")
#' d = 0:5
#' Sdz.mx = est_Sd(simu_truncAC, simu_truncAC, d, 'Cox', 'X', 'Q', 'delta', c("A", "Z1"))
est_Sd <- function(dat.fit, dat.est = dat.fit, time.eval, model, 
                   time.name, Q.name, event.name, cov.names = NULL, trim = 0, 
                   mtry = NA, ntree = NA, formula.survPen = NA, 
                   x.fit = NULL, x.est = NULL, cov.names.binary = NULL,
                   nfolds = 10, s = "lambda.min", alpha = 1, lambda = NULL, 
                   df = 5, OOF = FALSE, nfolds.OOF = 10){
    
    # u = time.eval
    
    names = c(time.name, Q.name, event.name, cov.names)
    if(is.na(match("Y", names)) & is.na(match("delta.D", names))){
        dat.fit$Y = dat.fit[,time.name] - dat.fit[,Q.name]    # residual censored event time
        dat.fit$delta.D = 1- dat.fit[,event.name]    # censoring indicator
    }else{
        stop("The names of the variables cannot be 'Y' or 'delta.D'. Conflict with intermediate variables in the function.")
    }
    
    
    if(model == "RF"){
        
        if(OOF){
            
            if(nrow(dat.fit)!=nrow(dat.est)){
                stop("The out-of-fold prediction implemented here is onlt for the case that dat.fit = dat.est.")
            }
            if(is.null(nfolds.OOF)){ nfolds.OOF = 10 }
            
            cov.names = c(cov.names, cov.names.binary)
            formula.D = formula(paste("Surv(Y, delta.D) ~ ", 
                                      paste(cov.names, collapse = " + "), collapse = ""))
            
            n.fit = nrow(dat.fit)
            folds.OOF <- cvFolds(n.fit, K)
            
            Sdz.mx = matrix(nrow = nrow(dat.fit), ncol = length(time.eval))
            ct = 0
            for(k in 1:nfolds.OOF){
                
                ddat.est = dat.fit[folds.OOF$subsets[folds.OOF$which == k], ]
                ddat.fit = dat.fit[folds.OOF$subsets[folds.OOF$which != k], ]
                
                rf_fit = rfsrc(formula.D, data = ddat.fit, ntree = ntree, mtry = mtry, 
                               splitrule = 'bs.gradient')
                newdata = ddat.est
                time.interest = predict(rf_fit, newdata = newdata)$time.interest
                surv_rf.mx = predict(rf_fit, newdata = newdata)$survival
                colnames(surv_rf.mx) = time.interest
                
                Fdz.mx0 = 1 - surv_rf.mx
                Sdz.mx_k = 1 - CDF_eval.mx(time.eval, Fdz.mx0)
                
                nk = nrow(ddat.est)
                Sdz.mx[(ct+1):(ct+nk),] = Sdz.mx_k
                ct = ct + nk
            }
            
            colnames(Sdz.mx) = time.eval
            
            
        }else{
            cov.names = c(cov.names, cov.names.binary)
            formula.D = formula(paste("Surv(Y, delta.D) ~ ", 
                                      paste(cov.names, collapse = " + "), collapse = ""))
            rf_fit = rfsrc(formula.D, data = dat.fit, ntree = ntree, mtry = mtry, 
                           splitrule = 'bs.gradient')
            
            newdata = dat.est
            time.interest = predict(rf_fit, newdata = newdata)$time.interest
            surv_rf.mx = predict(rf_fit, newdata = newdata)$survival
            colnames(surv_rf.mx) = time.interest
            
            Fdz.mx0 = 1 - surv_rf.mx
            Sdz.mx = 1 - CDF_eval.mx(time.eval, Fdz.mx0)
        }
        
        

    }else if(model == "Cox"){
        cov.names = c(cov.names, cov.names.binary)
        dat.cox = dat.fit
        formula.D = formula(paste("Surv(Y, delta.D) ~ ", 
                                  paste(cov.names, collapse = " + "), collapse = ""))
        # fit Cox-PH model for D|Z
        fit.D = coxph(formula.D, data = dat.cox)
        basehaz.D = basehaz(fit.D, centered = FALSE)
        beta.D = coef(fit.D)
        
        Z.D = as.matrix(dat.est[,cov.names])
        Sdz.mx = 1-Fz(time.eval, Z.D, basehaz.D, beta.D)
        # # Visualize the estimates --------------------------------------------
        # nplot = 8
        # plot(d, Sdz.mx[1,], type = "l", col = 1)
        # for(i in 2:nplot){
        #     lines(d, Sdz.mx[i,], col = 4+i)
        # }
        # # End visualize -------------------------------------------------------
        
    }else if(model == "pCox"){
        
        if(is.null(cov.names)){
            if(is.null(x.fit)|is.null(x.est)){
                stop("Need to input 'cov.names' if 'x.fit' or 'x.est' is NULL. ")
            }
        }else{
            dat.cmb = rbind(dat.fit[, c(cov.names, cov.names.binary)], dat.est[, c(cov.names, cov.names.binary)])
            XX = ns_mx(dat.cmb, cov.names, cov.names.binary, df)
            x.fit = XX[1:nrow(dat.fit), ]
            x.est = XX[-(1:nrow(dat.fit)), ]
        }
        
        yss = Surv(dat.fit$Y, dat.fit$delta.D)
        cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, alpha = alpha, lambda = lambda)
        fit.pred = survival::survfit(cv.fit, s = s, x = x.fit, y = yss, newx = x.est)
        Sdz.mx = t(fit.pred$surv)
        colnames(Sdz.mx) <- fit.pred$time
        
        Fdz.mx0 = 1 - Sdz.mx
        Sdz.mx = 1 - CDF_eval.mx(time.eval, Fdz.mx0)
        # # Visualize the estimates --------------------------------------------
        # nplot = 8
        # plot(d, Sdz.mx[1,], type = "l", col = 1)
        # for(i in 2:nplot){
        #     lines(d, Sdz.mx[i,], col = 4+i)
        # }
        # # End visualize -------------------------------------------------------

    }else if(model == "survPen"){
        
        cov.names = c(cov.names, cov.names.binary)
        n.u = length(time.eval)
        n.est = nrow(dat.est)
        
        # May want to modify how the residual time is coded. Y is involed in the formula, which need to be specified by user
        mod.D = survPen(formula.survPen, data = dat.fit, t1 = Y, event = delta.D)
        dat.est.u = data.frame(Y = rep(time.eval, n.est))
        for(j in cov.names){
            dat.est.u = cbind(dat.est.u, rep(dat.est[,j], each = n.u))
        }
        colnames(dat.est.u) = c("Y", cov.names)
        
        Sdz.mx = matrix(predict(mod.D, dat.est.u)$surv, ncol = n.u, byrow = TRUE)
        Sdz.mx = pmax(Sdz.mx, trim)
        colnames(Sdz.mx) = time.eval
        
    }else{
        
        stop("This S_D model is not implemented in this function!")
        
    }
    
    Sdz.mx = pmax(Sdz.mx, trim)
    return(Sdz.mx)
}


# return a matrix of true conditional survival probabilities for D given Z
truth_Sd <- function(dat.est, time.eval, 
                     cov.names, D.min, beta.C, shape.C){
    
    Sd.est <- function(t){
        surv = pweibull(rep(t - D.min, nrow(dat.est)), shape = shape.C, 
                        scale = exp(-1/shape.C * (cbind(1, as.matrix(dat.est[,cov.names]))%*%beta.C)),
                        lower = F)
        return(surv)
    }
    
    Sdz.mx = sapply(time.eval, Sd.est)
    colnames(Sdz.mx) <- time.eval
    
    return(Sdz.mx)
}




### function that output the estimated propensity scores for each subject ----------------


est_PS <- function(dat.fit, dat.est, model, A.name, cov.names, cov.names.binary = NULL, 
                   weights = rep(1, nrow(dat.fit)), 
                   trim = 0, df = 7, ntree = 2000){
    if(model == "logistic"){
        cov.names = c(cov.names, cov.names.binary)
        formula.A = formula(paste(A.name, "~", paste(cov.names, collapse = " + ")))
        
        # Use IPW.Q weights
        psfit <- glm(formula.A, family = "binomial", data = dat.fit, weights = weights)
        gamma = coef(psfit)
        PS = predict(psfit, newdata = dat.est, type = "response")
        
    }else if(model == "ns.logistic"){
        # logistic regression with natural splines basis
        if(length(cov.names)>1 | length(cov.names.binary)>0){
            stop("The logistic multivariate regression with natural splines is not implemented yet.")
        }
        
        Z_ns <- ns(dat.fit[, cov.names], df = df)
        colnames(Z_ns) = paste("Z", colnames(Z_ns), sep = "")
        dat.glm = as.data.frame(cbind(A = dat.fit[, A.name], Z_ns))
        model = glm(A~., family = binomial, data = dat.glm)
        
        # Generate natural spline basis for the covariates in dat.est
        Z_ns_new <- predict(Z_ns, newx = dat.est[, cov.names])
        colnames(Z_ns_new) = paste("Z", colnames(Z_ns_new), sep = "")
        PS <- predict(model, newdata = as.data.frame(Z_ns_new), type = "response")
        
    }else if(model =="gbm"){
        if(is.null(ntree)){ ntree = 2000 }
        cov.names = c(cov.names, cov.names.binary)
        formula.A = formula(paste(A.name, "~", paste(cov.names, collapse = " + ")))
        psfit = gbm(formula.A, data = dat.fit, distribution = 'bernoulli', 
                    weights = weights, n.trees = ntree)
     
        PS <- predict(psfit, newdata = dat.est, n.trees = ntree, type = 'response')
        
    }
    # else if(model == "pns.logistic"){
    #     # logistic regression with natural splines basis
    #     if(length(cov.names)>1){
    #         stop("The logistic multivariate regression with natural splines is not implemented yet.")
    #     }
    #     
    #     Z_ns <- ns(dat.fit[, cov.names], df = df)
    #     XX = cbind(dat.fit[, cov.names], Z_ns)
    #     
    #     ### Continue here ....... check if glmnet handle logistic regression
    # 
    # }
    else{
        stop("This model for A is not implemented yet.")
    }
    
    PS = pmin(pmax(PS, trim), 1-trim)
    return(PS)
}


truth_PS <- function(dat.est, cov.names, gamma.A){
    
    linpred = as.matrix(cbind(1,dat.est[,cov.names])) %*% gamma.A
    PS = as.vector(exp(linpred)/(1+exp(linpred)))
    
    return(PS)
}








### function that output the conditional CDF matrices for F and G at given time points ----------------
#' @title Estimate the Conditional CDF of the Event Time given Covariates
#' @description Estimate the conditional cumulative distribution function (CDF) of the event time given covariates evaluated at given time points. The options implemented in this function are: Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen'.
#' @param dat.fit data frame that is used to fit the model for the full data conditional distribution of the event time given the covariates.
#' @param dat.est data frame that contains the subjects for which the estimated conditional CDF is computed.
#' @param time.eval vector of time points at which the conditional CDF is evaluated.
#' @param model method used to estimate the conditional CDF. The options available are "Cox" and "survPen", corresponding to Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen', respectively.
#' @param time.name name of the event time variable.
#' @param Q.name name of the left truncation time variable.
#' @param event.name name of the event indicator.
#' @param cov.names vector of the names of covariates. For using splines in model = "pCox", this is the names of continuous covariates. 
#' @param trim constant for bounding the estimated conditional CDF from 1.
#' @param formula.survPen the formula when applying the hazard model with penalized splines implemented in \code{survPen::survPen}.
#' @param cov.names.binary The names of binary covariates when using splines in the model = 'pCox' option.
#' @param probs vector os probabilities at which the quantiles are used as knots for generating spline basis
#' @return \code{F_est()} returns a matrix of the estimated conditional CDF for subjects in `\code{data.est}' evaluated at the time points in the vector `\code{time.eval}'. Each row corresponds to a subject and each column corresponds to a time point. The column names of the matrix are the times in `\code{time.eval}'.
#' @export
#' @import survival stats survPen
#' @seealso  \code{\link{G_est}}
#' @examples
#' data("simu")
#' u = c(1, 1.5, 2, 2.5, 3, 3.5, 4)
#' Fuz.mx = est_F(simu, simu[1:10,], u, "Cox", "time", "Q", "delta", c("Z1","Z2"))
est_F <- function(dat.fit, dat.est = dat.fit, time.eval, model,
                  time.name, Q.name = NULL, event.name, cov.names, trim = 0,
                  weights = NULL, 
                  mtry = NULL, ntree = NULL, 
                  formula.survPen = NULL, 
                  x.fit = NULL, x.est = NULL, cov.names.binary = NULL,
                  nfolds = 10, s = "lambda.min", alpha = 1, lambda = NULL, 
                  df = 5){
    
    # if(is.null(trim)){ trim = 1e-7 }
    if(is.null(nfolds)){ nfolds = 10 }
    if(is.null(s)){ s = "lambda.min" }
    if(is.null(alpha)){ alpha = 1 }
    if(is.null(df)){ df = 5 }
    
    if((!is.null(weights)) & model == "RF"){
        stop("The function LTRCforests::ltrcrrf() does not take in case weights.")
    }
    if((!is.null(weights)) & model == "survPen"){
        stop("The function survPen::survPen() does not take in case weights.")
    }
    
    
    u = time.eval
    
    if(model == "RF"){
        cov.names = c(cov.names, cov.names.binary)
        F_hat = F_hat.ltrcrrf(dat.fit, time.name, Q.name, event.name, cov.names,
                              mtry, ntree)
        Fuz.mx = F_hat(newdata = dat.est, time.eval = u)   #***

    }else if(model == "Cox"){
        cov.names = c(cov.names, cov.names.binary)
        if(is.null(Q.name)){
            formula.T = formula(paste("Surv(", time.name, ",",
                                      event.name, ") ~ ",
                                      paste(cov.names, collapse = " + "), collapse = ""))
            id.cox = 1:nrow(dat.fit)
            dat.cox = dat.fit[id.cox, ]
        }else{
            formula.T = formula(paste("Surv(", Q.name, ",", time.name, ",",
                                      event.name, ") ~ ",
                                      paste(cov.names, collapse = " + "), collapse = ""))
            id.cox = ((dat.fit[,time.name] - dat.fit[,Q.name]) >= 10^{-7})
            dat.cox = dat.fit[id.cox, ]
        }
        
        # fit Cox-PH model for T|Z
        if(is.null(weights)){
            fit.T = coxph(formula.T, data = dat.cox)
        }else{
            fit.T = coxph(formula.T, data = dat.cox, weights = weights[id.cox])
        }
        
        basehaz.T = basehaz(fit.T, centered = FALSE)
        beta.T = coef(fit.T)
        
        Z.T = as.matrix(dat.est[,cov.names])
        Fuz.mx = Fz(u, Z.T, basehaz.T, beta.T)
        
        # # Visualize the estimates --------------------------------------------
        # nplot = 8
        # plot(time.eval, Fuz.mx[1,], type = "l", col = 1)
        # for(i in 2:nplot){
        #     lines(time.eval, Fuz.mx[i,], col = i)
        # }
        # # End visualize ------------------------------------------------------
        
    }else if(model == "survPen"){
        
        cov.names = c(cov.names, cov.names.binary)
        names = c(cov.names, time.name, Q.name, event.name)
        if(sum("delta.T"==names)>0){
            stop("The variables name 'delta.T' is used in the intermediate steps and cannot be used as inputs. ")
        }
    
        # time = dat.fit[ ,time.name]
        Q = dat.fit[ ,Q.name]
        dat.fit$delta.T = dat.fit[ ,event.name]
        n.u = length(u)
        n.est = nrow(dat.est)
        
        mod.T = survPen(formula.survPen, data = dat.fit, t0 = Q, t1 = X, event = delta.T)
        
        cov.mx = dat.est[,cov.names]
        cov.mx.rep = matrix(rep(t(cov.mx), n.u), ncol = ncol(cov.mx), byrow = TRUE)
        TZ.mx = cbind(rep(u, each = n.est), cov.mx.rep)
        colnames(TZ.mx) <- c(time.name, cov.names)
        dat.est.u = as.data.frame(TZ.mx)
        
        Fuz.mx = 1 - matrix(predict(mod.T, dat.est.u)$surv, ncol = n.u, byrow = FALSE)
        colnames(Fuz.mx) = u
        
    }else if(model == "pCox"){
        
        if(is.null(cov.names)){
            if(is.null(x.fit)|is.null(x.est)){
                stop("Need to input 'cov.names' if 'x.fit' or 'x.est' is NULL. ")
            }
        }else{
            dat.cmb = rbind(dat.fit[, c(cov.names, cov.names.binary)], dat.est[, c(cov.names, cov.names.binary)])            
            XX = ns_mx(dat.cmb, cov.names, cov.names.binary, df)
            x.fit = XX[1:nrow(dat.fit), ]
            x.est = XX[-(1:nrow(dat.fit)), ]
        }
        
        if(is.null(Q.name)){
            yss = Surv(dat.fit[,time.name], dat.fit[,event.name])
        }else{
            yss = Surv(dat.fit[,Q.name], dat.fit[,time.name], dat.fit[,event.name])
        }
        
        if(is.null(weights)){
            cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, 
                                alpha = alpha, lambda = lambda)
            fit.pred = survival::survfit(cv.fit, s = s, x = x.fit, y = yss, newx = x.est,
                                         alpha = alpha, lambda = lambda)
        }else{
            cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, weights = weights, 
                                alpha = alpha, lambda = lambda)
            fit.pred = survival::survfit(cv.fit, s = s, x = x.fit, y = yss, newx = x.est,
                                         weights = weights, alpha = alpha, lambda = lambda)
        }
        
        Fuz.mx0 = 1 - t(fit.pred$surv)
        colnames(Fuz.mx0) <- fit.pred$time
        
        Fuz.mx = CDF_eval.mx(time.eval, Fuz.mx0)
        # # Visualize the estimates --------------------------------------------
        # nplot = 8
        # plot(time.eval, Fuz.mx[1,], type = "l", col = 1)
        # for(i in 2:nplot){
        #     lines(time.eval, Fuz.mx[i,], col = i)
        # }
        # # End visualize ------------------------------------------------------
        
    }else{
        stop("This T model is not implemented in this function!")
    }
    
    Fuz.mx = pmin(Fuz.mx, 1-trim)
    
    return(Fuz.mx)
    # return(list(Fuz.mx = Fuz.mx, beta.T = beta.T))
}


# computing the truth F when T^* is simulated from a Weibull distribution
truth_F <- function(dat.est, time.eval, 
                    cov.names, T.min, beta.T, shape.T){
    
    F.est <- function(t){
        cumpr = pweibull(rep(t-T.min, nrow(dat.est)), shape = shape.T, 
                        scale = exp(-1/shape.T * (cbind(1, as.matrix(dat.est[,cov.names]))%*%beta.T)),
                        lower = TRUE)
        return(cumpr)
    }
    
    Fuz.mx = sapply(time.eval, F.est)
    colnames(Fuz.mx) <- time.eval
    
    return(Fuz.mx)
}


# Compute the truth F when T^* is simulated from AFT model
truth_F.AFT <- function(dat.est, time.eval, 
                        cov.names, beta.T, 
                        err.dist = c("weibull", "lognormal"),
                        err.sigma = NULL, err.shape = NULL, err.scale = NULL){
    
    y_pred = cbind(1, as.matrix(dat.est[,cov.names]))%*%beta.T
    
    if(err.dist == "lognormal"){
        F.est <- function(t){
            cumpr = pnorm(log(t) - y_pred, mean = 0, sd = err.sigma,
                          lower = TRUE)
            return(cumpr)
        }
    }else if(err.dist == "weibull"){
        F.est <- function(t){
            cumpr = pweibull(log(t) - y_pred + err.scale*gamma(1+1/err.shape), 
                             shape = err.shape, scale = err.scale,
                             lower = TRUE)
            return(cumpr)
        }
    }else{
        stop("This distribution is not implemented in the function yet.")
    }
    
    
    Fuz.mx = sapply(time.eval, F.est)
    colnames(Fuz.mx) <- time.eval
    
    return(Fuz.mx)
    
}




## 2024-12-07: Mimicking the spline basis in Nie & Wager code
ns_mx <- function(data, cov.names, cov.names.binary, df){

    X_ns = do.call(cbind, lapply(cov.names, function(cov){splines::ns(data[, cov], df = df)}))
    X_ns_wbinary = cbind(X_ns, data[, cov.names.binary])  # combine with the binary covariates
    dim_ns = dim(X_ns)[2]
    X_ns_linear_interaction = stats::model.matrix(~.*.-1, data.frame(X_ns_wbinary)) # pairwise interaction (not including squared term for each column), including the first order terms
    X_ns_sq = do.call(cbind, lapply(1:dim_ns, function(col){matrix(X_ns[,col]^2)})) # squared term for each column
    X_ns = cbind(X_ns_linear_interaction, X_ns_sq)

    return(X_ns)
}


# ## Oracle basis for HTE4 simulation setup -- only for testing purposes
# ns_mx <- function(data, cov.names, cov.names.binary, df){
#     
#     X_ns_conti_binary = cbind(data[, cov.names], data[, cov.names.binary])  # combine with the binary covariates
#     X_ns_sin_sqrt = cbind(sin(pi*data$Z1), sin(pi*data$Z2), sqrt(abs(data$Z1)), sqrt(abs(data$Z2)))
#     X_ns_interaction = cbind(data$A*data$Z1, data$A*data$Z2, 
#                              data$A*sin(pi*data$Z1), data$A*sin(pi*data$Z2),
#                              data$A*sqrt(abs(data$Z1)), data$A*sqrt(abs(data$Z2)))
#     X_ns = as.matrix(cbind(X_ns_conti_binary, 
#                            X_ns_sin_sqrt,
#                            X_ns_interaction))
#     
#     return(X_ns)
# }




# # Generate natural cubic spline basis and their pairwise interactions - old function -- before 2024-12-07
# ns_mx <- function(data, cov.names, cov.names.binary, df){
# 
#     dat.cmb = data
# 
# 
#     # Create spline basis functions for continuous covariates
#     spline_basis_list_continuous <- lapply(cov.names, function(cov) {
#         return(ns(dat.cmb[,cov], df = df, intercept = TRUE))   ## probably should change to false
#     })
# 
#     # Create pairwise interactions for continuous covariates
#     interaction_matrices_continuous <- lapply(seq_along(spline_basis_list_continuous), function(i) {
#         do.call(cbind, lapply(seq_along(spline_basis_list_continuous), function(j) {
#             if (i < j) {
#                 comb.mx = expand.grid(i = 1:ncol(spline_basis_list_continuous[[i]]), j = 1:ncol(spline_basis_list_continuous[[j]]))
#                 apply(comb.mx, MARGIN = 1, FUN = function(ids){spline_basis_list_continuous[[i]][,ids[1]]*spline_basis_list_continuous[[j]][,ids[2]]})
#             } else {
#                 NULL
#             }
#         }))
#     })
# 
#     # Create pairwise interactions between continuous and binary covariates
#     interaction_matrices_continuous_binary <- lapply(seq_along(spline_basis_list_continuous), function(i) {
#         do.call(cbind, lapply(seq_along(cov.names.binary), function(j) {
#             n.replicate = ncol(spline_basis_list_continuous[[i]])
#             spline_basis_list_continuous[[i]] * replicate(n.replicate, dat.cmb[, cov.names.binary[j]] )
#         }))
#     })
# 
#     # Create pairwise interactions for binary covariates
#     interaction_matrices_binary <- lapply(seq_along(cov.names.binary), function(i) {
#         do.call(cbind, lapply(seq_along(cov.names.binary), function(j) {
#             if (i < j) {
#                 dat.cmb[, cov.names.binary[i]] * dat.cmb[, cov.names.binary[j]]
#             } else {
#                 NULL
#             }
#         }))
#     })
# 
#     matrices_binary = list(dat.cmb[, c(cov.names.binary)])
# 
#     # Combine all the interaction terms and covariates
#     XX <- do.call(cbind, c(
#         spline_basis_list_continuous,
#         interaction_matrices_continuous,
#         interaction_matrices_continuous_binary,
#         interaction_matrices_binary,
#         matrices_binary
#     ))
#     
#     # XX = cbind(XX, 
#     #            sin(pi*data$Z1),
#     #            sin(pi*data$Z2),
#     #            sqrt(abs(data$Z1)), 
#     #            sqrt(abs(data$Z2)),
#     #            data$A * sin(pi*data$Z1),
#     #            data$A * sin(pi*data$Z2),
#     #            data$A * sqrt(abs(data$Z1)), 
#     #            data$A * sqrt(abs(data$Z2)))
# 
#     return(XX)
# }





#' @title Estimate the Conditional CDF for the Left Truncation Time given Covariates
#' @description Estimate the conditional cumulative distribution function (CDF) of the left truncation time given covariates evaluated at given time points. The options implemented in this function are: Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen'.
#' @param dat.fit data frame that is used to fit the model for the full data conditional distribution of the event time given the covariates.
#' @param dat.est data frame that contains the subjects for which the estimated conditional CDF is computed.
#' @param time.eval vector of time points at which the conditional CDF is evaluated.
#' @param model method used to estimate the conditional CDF. The options available are "Cox" and "survPen", corresponding to Cox proportional hazards regression using function \code{coxph()} from  R package `survival', and the hazard model with penalized splines using function \code{survPen()} from R package `survPen', respectively.
#' @param time.name name of the event time variable.
#' @param Q.name name of the left truncation time variable.
#' @param event.name name of the event indicator.
#' @param cov.names vector of the names of covariates.
#' @param trim constant for bounding the estimated conditional CDF from 0.
#' @param weights vector of case weights.
#' @param trunc a logical indicator indicating whether there is truncation or not. If \code{trunc = TRUE}, then the reversed time technique is applied.
#' @param formula.survPen the formula when applying the hazard model with penalized splines implemented in \code{survPen::survPen}.
#' @return \code{G_est()} returns a matrix of the estimated conditional CDF for subjects in \code{data.est} evaluated at the time points in the vector \code{time.eval}. Each row corresponds to a subject and each column corresponds to a time point. The column names of the matrix are the times in `\code{time.eval}'.
#' @export
#' @import survival stats survPen
#' @seealso  \code{\link{F_est}}
#' @examples
#' data("simu")
#' v = c(0.5, 1, 1.5, 2, 2.5, 3)
#' Gvz.mx = est_G(simu, simu[1:10,], v, "Cox", "time", "Q", "delta", c("Z1","Z2"))
est_G <- function(dat.fit, dat.est = dat.fit, time.eval, model,
                  time.name, Q.name, event.name, cov.names, trim = 0, weights = NULL,
                  formula.survPen = NULL, tau = NULL, trunc = TRUE,
                  x.fit = NULL, x.est = NULL, cov.names.binary = NULL,
                  nfolds = 10, s = "lambda.min", alpha = 1, lambda = NULL, df = 5){
    
    if((!is.null(weights)) & model == "RF"){
        stop("The function LTRCforests::ltrcrrf() does not take in case weights.")
    }
    if((!is.null(weights)) & model == "survPen"){
        stop("The function survPen::survPen() does not take in case weights.")
    }
    
    if(mean(dat.fit[,event.name])<1){
        stop("The truncation time is always observed, so dat.fit[,event.name] should be a vector of one's.")
    }
    
    names = c(time.name, Q.name, event.name, cov.names)
    if(sum(names == "Q2")){
        stop("The names of the variables cannot be 'Q2'. It is used in the middle of the computation.")
    }
    if(sum(names == "T2")){
        stop("The names of the variables cannot be 'T2'. It is used in the middle of the computation.")
    }
    
    v = time.eval
    
    if(is.null(tau)){
        tau = max(c(dat.fit[,time.name], dat.est[,time.name])) + 1
    }
    
    dat.fit$Q2 = tau - dat.fit[ ,Q.name]
    dat.fit$T2 = tau - dat.fit[ ,time.name]
    dat.fit$delta.1 = rep(1, nrow(dat.fit))
    
    dat.est$Q2 = tau - dat.est[ ,Q.name]
    dat.est$T2 = tau - dat.est[ ,time.name]
    dat.est$delta.1 = rep(1, nrow(dat.est))
    
    # if(model == "RF"){
    # 
    #     G_hat = G_hat.ltrcrrf(dat.fit, time.name, Q.name, event.name, cov.names,
    #                           mtry, ntree, tau)
    #     Gvz.mx = pmax(G_hat(newdata = dat.est, time.eval = v), trim)   #***
    # 
    # }else
    if(model == "Cox"){
        
        cov.names = c(cov.names, cov.names.binary)
        id.cox = ((dat.fit[, time.name] - dat.fit[, Q.name]) >= 10^{-7})
        dat.cox = dat.fit[id.cox, ]
        Z.Q = as.matrix(dat.est[,cov.names])
        
        if(trunc){
            formula.Q2 = formula(paste("Surv(tau-", time.name, ", tau-", Q.name,",", event.name, ") ~ ",
                                       paste(cov.names, collapse = " + "), collapse = ""))
        }else{
            formula.Q2 = formula(paste("Surv(Q2,", event.name, ") ~ ",
                                       paste(cov.names, collapse = " + "), collapse = ""))
        }
        
        
        # fit Cox-PH model for (tau-Q)|Z
        if(is.null(weights)){
            fit.Q2 = coxph(formula.Q2, data = dat.cox)
        }else{
            fit.Q2 = coxph(formula.Q2, data = dat.cox, weights = weights[id.cox])
        }
        
        basehaz.Q2 = basehaz(fit.Q2, centered = FALSE)
        beta.Q2 = coef(fit.Q2)
        
        Gvz.mx = Gz(v, Z.Q, basehaz.Q2, beta.Q2, tau)
    
        
    }else if(model == "pCox"){
        
        if(is.null(cov.names)){
            if(is.null(x.fit)|is.null(x.est)){
                stop("Need to input 'cov.names' if 'x.fit' or 'x.est' is NULL. ")
            }
        }else{
            dat.cmb = rbind(dat.fit[, c(cov.names, cov.names.binary)], dat.est[, c(cov.names, cov.names.binary)])            
            XX = ns_mx(dat.cmb, cov.names, cov.names.binary, df)
            x.fit = XX[1:nrow(dat.fit), ]
            x.est = XX[-(1:nrow(dat.fit)), ]
        }
        
        if(trunc){
            yss = Surv(dat.fit$T2, dat.fit$Q2, dat.fit$delta.1)
        }else{
            yss = Surv(dat.fit$Q2, dat.fit$delta.1)
        }
        
        
        if(is.null(weights)){
            cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, alpha = alpha, lambda = lambda)
        }else{
            cv.fit <- cv.glmnet(x.fit, yss, family = "cox", nfolds = nfolds, weights = weights, alpha = alpha, lambda = lambda)
        }
        
        fit.pred = survival::survfit(cv.fit, s = s, x = x.fit, y = yss, newx = x.est, 
                                     weights = weights, alpha = alpha, lambda = lambda)
        Gvz.mx0 = t(fit.pred$surv)
        colnames(Gvz.mx0) <- tau - fit.pred$time
        id0 = order(tau-fit.pred$time)
        Gvz.mx0 = Gvz.mx0[,id0]
        
        Gvz.mx = CDF_eval.mx(time.eval, Gvz.mx0)
        # # Visualize the estimates --------------------------------------------
        # nplot = 8
        # plot(time.eval, Gvz.mx[1,], type = "l", col = 1)
        # for(i in 2:nplot){
        #     lines(time.eval, Gvz.mx[i,], col = 4+i)
        # }
        # # End visualize ------------------------------------------------------
        
    
    }else if(model == "survPen"){
        
        cov.names = c(cov.names, cov.names.binary)
        T2 = dat.fit$T2
        Q2 = dat.fit$Q2
        delta.Q = dat.fit$delta.1
        
        n.v = length(v)
        n.est = nrow(dat.est)
        
        mod.Q2 = survPen(formula.survPen, data = dat.fit, t0 = T2, t1 = Q2, event = delta.Q)
        
        cov.mx = dat.est[,cov.names]
        cov.mx.rep = matrix(rep(t(cov.mx), n.v), ncol = ncol(cov.mx), byrow = TRUE)
        Q2Z.mx = cbind(tau - rep(v, each = n.est), cov.mx.rep)
        colnames(Q2Z.mx) <- c("Q2", cov.names)
        dat.est.v = as.data.frame(Q2Z.mx)
        
        Gvz.mx = matrix(predict(mod.Q2, dat.est.v)$surv, ncol = n.v, byrow = FALSE)
        
        
    }else{
        stop("This Q model is not implemented in this function!")
    }
    
    colnames(Gvz.mx) = v
    Gvz.mx = pmax(Gvz.mx, trim)
    
    return(Gvz.mx)
    # return(list(Gvz.mx = Gvz.mx, b.Q = beta.Q2))
}



truth_G <- function(dat.est, time.eval, cov.names,
                    Q.min, Q.max, tau, beta.Q){
    
    Q2.min = tau-Q.max
    Q2.max = tau-Q.min
    
    ### continue here - debug
    G.est <- function(t){
        t2 = tau - t
        cumpr = (1-pmin(pmax((t2 - Q2.min)/(Q2.max - Q2.min),0),1))^(1/exp(-cbind(1,as.matrix(dat.est[, cov.names])) %*% beta.Q))
        return(cumpr)
    }
    
    Gvz.mx = sapply(time.eval, G.est)
    colnames(Gvz.mx) <- time.eval
    
    return(Gvz.mx)
}









### Functions needed when using Cox models ---------------------------------------
# function for computing the baseline CDF or survival function
baseCDF.single <- function(t, basehaz){
    if(t<min(basehaz[,'time'])){
        return(0)   # check 1 or 0
    }else{
        id = which((t - basehaz[,'time'])>=0)
        index = which.min((t - basehaz[,'time'])[id])
        
        return(1-exp(-basehaz[index,'hazard']))
    }
}

baseCDF <- function(t, basehaz){
    cdf = sapply(t, baseCDF.single, basehaz = basehaz)
    return(cdf)
}

baseS.single <- function(t, basehaz){
    if(t<min(basehaz[,'time'])){
        return(1)   # check 1 or 0
    }else{
        id = which((t - basehaz[,'time'])>=0)
        index = which.min((t - basehaz[,'time'])[id])
        
        return(exp(-basehaz[index,'hazard']))
    }
}

baseS <- function(t, basehaz){
    S = sapply(t, baseS.single, basehaz = basehaz)
    return(S)
}


## function for computing G(t|z)
Gz <- function(t, Z, basehaz.Q2, beta.Q2, tau){
    if(is.matrix(Z)){
        result = matrix(nrow = nrow(Z), ncol = length(t))
        for(j in 1:length(t)){
            if(tau-t[j]<min(basehaz.Q2$time)){
                result[,j] = 1
            }else{
                result[,j] = (baseS(tau-t[j], basehaz.Q2))^exp(Z %*% beta.Q2)
            }
        }
    }else{
        result = matrix(nrow = 1, ncol = length(t))
        for(j in 1:length(t)){
            if(tau-t[j]<min(basehaz.Q2$time)){
                result[,j] = 1
            }else{
                result[,j] = (baseS(tau-t[j], basehaz.Q2))^exp(sum(Z*beta.Q2))
            }
        }
    }
    colnames(result) = t
    
    return(result)
}

## function for computing F(t|z)
Fz <- function(t, Z, basehaz.T, beta.T){
    if(is.matrix(Z)){
        result = matrix(nrow = nrow(Z), ncol = length(t))
        for(j in 1:length(t)){
            result[,j] = 1-(baseS(t[j], basehaz.T))^exp(Z %*% beta.T)
        }
    }else{
        result = matrix(nrow = 1, ncol = length(t))
        for(j in 1:length(t)){
            result[,j] = 1-(baseS(t[j], basehaz.T))^exp(sum(Z*beta.T))
        }
    }
    colnames(result) = t
    
    return(result)
}














### Functions needed when using LTRCforests::ltrcrrf -------------------------------------

# returns a function F_hat that takes 'newdata' as the input, and outputs a matrix
# of the estimated survival probabilities for each subject in 'newdata' (row)
# at the times in time.eval (col).

F_hat.ltrcrrf <- function(dat, time.name, Q.name, event.name, cov.names,
                          mtry, ntree){

    # Z.names = c(A.names, cov.names)
    Z.names = cov.names
    formula.T = formula(paste(paste("Surv(", Q.name,",", time.name, ",", event.name,") ~ ", collapse = ""),
                              paste(Z.names, collapse = "+"), collapse = ""))

    T.obj <- ltrcrrf(formula = formula.T, data = dat, mtry = mtry, ntree = ntree)
    jumps.T = sort(dat[,time.name])

    F_hat <- function(newdata, time.eval = jumps.T){ # by default, evaluated at jumps.T
        newdata[, Q.name] = 0
        T.pred = predictProb(object = T.obj, newdata = newdata, time.eval = time.eval)  #**
        survProb.mx = t(T.pred$survival.probs)
        rownames(survProb.mx) = rownames(newdata)
        colnames(survProb.mx) = time.eval

        return(1-survProb.mx)
    }

    return(F_hat)

}


# returns a function F_hat that takes 'newdata' as the input, and outputs a matrix
# of the estimated survival probabilities for each subject in 'newdata' (row)
# at the times in time.eval (col).

G_hat.ltrcrrf <- function(dat, time.name, Q.name, event.name, cov.names,
                          mtry, ntree, tau = max(dat[,time.name])+1){

    if(tau < max(dat[,time.name])+1e-10){
        stop("Need a larger tau.")
    }

    names = c(time.name, Q.name, event.name, cov.names)
    if(is.na(match("Q2", names)) & is.na(match("T2", names))){
        dat$Q2 = tau - dat[,Q.name]
        dat$T2 = tau - dat[,time.name]
    }else{
        stop("The names of the variables cannot be 'T2' or 'Q2'.")
    }

    # Z.names = c(A.names, cov.names)
    Z.names = cov.names
    formula.Q2 = formula(paste(paste("Surv(T2, Q2,", event.name,") ~ ", collapse = ""),
                               paste(Z.names, collapse = "+"), collapse = ""))

    Q2.obj <- ltrcrrf(formula = formula.Q2, data = dat, mtry = mtry, ntree = ntree)
    jumps.Q = sort(dat[, Q.name])

    G_hat <- function(newdata, time.eval = jumps.Q){ # by default, evaluated at jumps.Q
        time.eval2 = tau-time.eval-1e-10
        if(min(time.eval2)<0){
            stop("Need a larger tau.")
        }

        newdata$Q2 = tau - newdata[,Q.name]
        newdata$T2 = 0

        odr = order(time.eval2)
        Q.pred = predictProb(object = Q2.obj, newdata = newdata,
                             time.eval = (time.eval2)[odr])   #**
        odr2 = match((time.eval2)[odr], time.eval2)
        prob.mx = t(Q.pred$survival.probs)[,odr2]
        rownames(prob.mx) = rownames(newdata)
        colnames(prob.mx) = time.eval

        return(prob.mx)
    }

    return(G_hat)
}


