### truncAC_AIPW


# DR estimator for triple biases: confounding, left truncation, right censoring
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
    
    # if(sum(names == "MBweights")){
    #     stop("The names of the variables cannot be 'MBweights'. It is used as the name for the weights involved in multiplier bootstrap")
    # }
    # 
    # if(mBootstrap){
    #     mbw = rexp(nrow(dat))
    #     dat$MBweights = nrow(dat)*mbw/sum(mbw)
    # }else{
    #     dat$MBweights = rep(1, nrow(dat))
    # }
    
    
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
    
    
    # Previous code that output teuncAC_AIPW, truncA_AIPW_IPCW, and IPW estimators
    # est_truncAC_AIPW = truncAC_AIPW(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, PS, mu_A1, mu_A0, 
    #                                 X.name, Q.name, event.name.T, A.name, trim)
    # 
    # return(est_truncAC_AIPW)
    
    
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
    Delta = dat[, event.name]
    GXZ_hat = CDF_eval(dat.est[,X.name], Gvz.mx)
    SdXZ_hat = 1 - CDF_eval(dat.est[,X.name] - dat.est[,Q.name], 1-Sdz.mx)
    
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
    
    # est = rbind(truncAC_AIPW = est_truncAC_AIPW$truncAC_AIPW,
    #             truncA_AIPW_IPCW = est_truncAC_AIPW$truncA_AIPW_IPCW,
    #             IPW = est_truncAC_AIPW$IPW)
    # est.a1 = est[,1]
    # est.a0 = est[,2]
    # 
    # est.a1a0 = c(est.a1, est.a0)
    # names(est.a1a0) = c(paste(rownames(est), ".a1", sep = ""),
    #                     paste(rownames(est), ".a0", sep = ""))
    # 
    # return(est.a1a0)
    
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




## after updating the dr_truncAC_AIPW() function, the `truncAC_AIPW()` function seems to be not used. --  Need to double check if it is used in `truncC_AIPW()`
truncAC_AIPW <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, PS, mu_A1, mu_A0, 
                         X.name, Q.name, event.name, A.name, trim = 1e-7){
    
    truncAC_EF = truncAC_AIPW_transMean_EF(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, PS, mu_A1, mu_A0,
                                           X.name, Q.name, event.name, A.name, trim)

    # truncAC_AIPW
    Num_truncAC_AIPW.a1 = truncAC_EF$Num_truncAC_AIPW.a1
    Num_truncAC_AIPW.a0 = truncAC_EF$Num_truncAC_AIPW.a0
    Den_truncAC_AIPW = truncAC_EF$Den_truncAC_AIPW

    est_truncAC_AIPW.a1 = mean(Num_truncAC_AIPW.a1)/mean(Den_truncAC_AIPW)
    est_truncAC_AIPW.a0 = mean(Num_truncAC_AIPW.a0)/mean(Den_truncAC_AIPW)

    est_truncAC_AIPW = c(est.a1 = est_truncAC_AIPW.a1, est.a0 = est_truncAC_AIPW.a0)


    # truncA_AIPW_IPCW
    Num_truncA_AIPW_IPCW.a1 = truncAC_EF$Num_truncA_AIPW_IPCW.a1
    Num_truncA_AIPW_IPCW.a0 = truncAC_EF$Num_truncA_AIPW_IPCW.a0
    Den_truncA_AIPW_IPCW = truncAC_EF$Den_truncA_AIPW_IPCW

    est_truncA_AIPW_IPCW = c(est.a1 = mean(Num_truncA_AIPW_IPCW.a1)/mean(Den_truncA_AIPW_IPCW),
                             est.a0 = mean(Num_truncA_AIPW_IPCW.a0)/mean(Den_truncA_AIPW_IPCW))


    # IPW
    Num_IPW.a1 = truncAC_EF$Num_IPW.a1
    Num_IPW.a0 = truncAC_EF$Num_IPW.a0
    Den_IPW.a1 = truncAC_EF$Den_IPW.a1
    Den_IPW.a0 = truncAC_EF$Den_IPW.a0

    est_IPW = c(est.a1 = mean(Num_IPW.a1)/mean(Den_IPW.a1),
                est.a0 = mean(Num_IPW.a0)/mean(Den_IPW.a0))

    
    
    # ### Debug - change to compute using truncC_AIPW -----------------------------------------
    # truncC_AIPW_result <- truncC_AIPW_transMean(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
    #                                             X.name, Q.name, event.name, trim)
    # 
    # # !!! save the temp results
    # dr_truncC_AIPW_result = truncC_AIPW_result
    # dr_truncC_details = list(dat = dat, nu = nu, 
    #                          Fuz.mx = Fuz.mx, Fuz_A1.mx = Fuz_A1.mx, Fuz_A0.mx = Fuz_A0.mx,
    #                          Gvz.mx = Gvz.mx, Sdz.mx = Sdz.mx, PS = PS, 
    #                          mu_A1 = mu_A1, mu_A0 = mu_A0) 
    # save(dr_truncC_AIPW_result, dr_truncC_details, 
    #      file = "debug_results/dr_truncC_AIPW_results.rda")
    # 
    # 
    # truncC_AIPW_1 = truncC_AIPW_result$truncC_const1
    # truncC_AIPW_nuT = truncC_AIPW_result$truncC_nu
    # 
    # # truncAC_AIPW
    # A = dat[, A.name]
    # PS_hat = PS
    # mu1_hat = mu_A1
    # mu0_hat = mu_A0
    # Num.a1 = A/pmax(PS_hat, trim) * truncC_AIPW_nuT - (A-PS_hat)/pmax(PS_hat, trim) * mu1_hat * truncC_AIPW_1
    # Num.a0 = (1-A)/pmax(1-PS_hat, trim) * truncC_AIPW_nuT + (A-PS_hat)/pmax(1-PS_hat, trim) * mu0_hat * truncC_AIPW_1
    # Den = truncC_AIPW_1
    # 
    # est_truncAC = c(est.a1 = mean(Num.a1)/mean(Den),
    #                 est.a0 = mean(Num.a0)/mean(Den))
    # est_truncAC
    # 
    # ## They give the same results
    # 
    # ### End debug ---------------------------------------------------------------------------
    # 
    
    
    return(list(truncAC_AIPW = est_truncAC_AIPW,
                truncA_AIPW_IPCW = est_truncA_AIPW_IPCW,
                IPW = est_IPW))
}


truncAC_AIPW_transMean_EF <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, PS, mu_A1, mu_A0, 
                                      X.name, Q.name, event.name, A.name, trim = 1e-7){
    
    X.res = dat[,X.name] - dat[,Q.name]
    A = dat[, A.name]
    
    efs = truncAIPW_transMean_EF(dat, nu, Fuz.mx, Gvz.mx, X.name, Q.name, trim)
    Num_AIPW = efs$Num_AIPW
    Den_AIPW = efs$Den_AIPW
    
    Num_IPW.Q = efs$Num_IPW.Q
    Den_IPW.Q = efs$Den_IPW.Q
    
    # compute the IPTW weights
    id_A1 = (A == 1)
    id_A0 = (A == 0)
    w_IPTW.1 = rep(0, nrow(dat))
    w_IPTW.0 = rep(0, nrow(dat))
    w_IPTW.1[id_A1] = 1/pmax(PS[id_A1], trim)      # possibly considering also bounding PS in the program? - added 2025-02-08
    w_IPTW.0[id_A0] = 1/pmax(1-PS[id_A0], trim)    # also for Num_truncAC_3.a1 and Num_truncAC_3.a0 - added long before 2025-02-08
    
    # compute the IPCW weights
    Sdyz.vec = diag(1 - CDF_eval.mx(X.res, 1-Sdz.mx))
    id1 = (dat[, event.name] == 1)
    w_IPCW = rep(0, nrow(dat))
    w_IPCW[id1] = 1/pmax(Sdyz.vec[id1], trim)
    
    # The terms for the numerator
    Num_truncAC_1.a1 = w_IPTW.1*w_IPCW*Num_AIPW
    Num_truncAC_1.a0 = w_IPTW.0*w_IPCW*Num_AIPW
    
    Num_truncAC_3.a1 = - (A-PS)/pmax(PS, trim) * mu_A1 * w_IPCW * Den_AIPW
    Num_truncAC_3.a0 = (A-PS)/pmax(1-PS, trim) * mu_A0 * w_IPCW * Den_AIPW
    
    result_Aug_QD = aug_QD(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, Sdyz.vec,
                           X.name, Q.name, event.name, trim)
    Aug_QD_nu = result_Aug_QD$Aug_QD_nu
    Aug_QD_const1 = result_Aug_QD$Aug_QD_const1
   
    
    # Num2 and Num4
    Num_truncAC_2.a1 = w_IPTW.1 * Aug_QD_nu
    Num_truncAC_2.a0 = w_IPTW.0 * Aug_QD_nu
    
    Num_truncAC_4.a1 = - (A-PS)/pmax(PS, trim) * mu_A1 * Aug_QD_nu
    Num_truncAC_4.a0 = (A-PS)/pmax(1-PS, trim) * mu_A0 * Aug_QD_nu
    
    # the terms for the denominator
    Den_truncAC_1 = w_IPCW * Den_AIPW
    Den_truncAC_2 = Aug_QD_const1
    
    ## compute the numerator and denominator of the estimators
    # truncAC_AIPW
    Num_truncAC_AIPW.a1 = Num_truncAC_1.a1 + Num_truncAC_2.a1 + Num_truncAC_3.a1 + Num_truncAC_4.a1
    Num_truncAC_AIPW.a0 = Num_truncAC_1.a0 + Num_truncAC_2.a0 + Num_truncAC_3.a0 + Num_truncAC_4.a0
    Den_truncAC_AIPW = Den_truncAC_1 + Den_truncAC_2
    
    # truncA_AIPW_IPCW
    Num_truncA_AIPW_IPCW.a1 = Num_truncAC_1.a1 + Num_truncAC_3.a1
    Num_truncA_AIPW_IPCW.a0 = Num_truncAC_1.a0 + Num_truncAC_3.a0
    Den_truncA_AIPW_IPCW = Den_truncAC_1
    
    
    # IPW
    Num_IPW.a1 = w_IPTW.1 * w_IPCW * Num_IPW.Q
    Num_IPW.a0 = w_IPTW.0 * w_IPCW * Num_IPW.Q
    Den_IPW.a1 = w_IPTW.1 * w_IPCW * Den_IPW.Q
    Den_IPW.a0 = w_IPTW.0 * w_IPCW * Den_IPW.Q
    
    return(list(Num_truncAC_AIPW.a1 = Num_truncAC_AIPW.a1, 
                Num_truncAC_AIPW.a0 = Num_truncAC_AIPW.a0,
                Den_truncAC_AIPW = Den_truncAC_AIPW,
                Num_truncA_AIPW_IPCW.a1 = Num_truncA_AIPW_IPCW.a1,
                Num_truncA_AIPW_IPCW.a0 = Num_truncA_AIPW_IPCW.a0,
                Den_truncA_AIPW_IPCW = Den_truncA_AIPW_IPCW,
                Num_IPW.a1 = Num_IPW.a1, Num_IPW.a0 = Num_IPW.a0,
                Den_IPW.a1 = Den_IPW.a1, Den_IPW.a0 = Den_IPW.a0))
}




aug_QD <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, Sdyz.vec = NA,
                   X.name, Q.name, event.name, trim = 1e-7){
    if(sum(is.na(Sdyz.vec))>0){
        X.res = dat[,X.name] - dat[,Q.name]
        Sdyz.vec = diag(1 - CDF_eval.mx(X.res, 1-Sdz.mx))
    }
    
    aug_result_1 = aug_QD_1(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, X.name, Q.name, event.name, trim)
    aug_result_2 = aug_QD_2(dat, nu, Fuz.mx, Gvz.mx, Sdyz.vec, Q.name, event.name, trim)
    aug_result_3 = aug_QD_3(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, X.name, Q.name, event.name, trim)
    
    Aug_QD_nu_1 = aug_result_1$aug_1
    Aug_QD_nu_2 = aug_result_2$aug_2
    Aug_QD_nu_3 = aug_result_3$aug_3
    Aug_QD_nu = Aug_QD_nu_1 + Aug_QD_nu_2 - Aug_QD_nu_3
    
    Aug_QD_const1_1 = aug_result_1$aug_1_const1
    Aug_QD_const1_2 = aug_result_2$aug_2_const1
    Aug_QD_const1_3 = aug_result_3$aug_3_const1
    Aug_QD_const1 = Aug_QD_const1_1 + Aug_QD_const1_2 - Aug_QD_const1_3
    
    # Aug_11 = aug_result_1$aug_11
    # Aug_21 = aug_result_1$aug_21
    # Aug_31 = aug_result_1$aug_31
    
    return(list(Aug_QD_nu = Aug_QD_nu,
                Aug_QD_const1 = Aug_QD_const1))
}



aug_QD_1 <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                     X.name, Q.name, event.name, trim = 1e-7){
    
    nn = nrow(dat)
    X = dat[,X.name]
    Q = dat[,Q.name]
    X.res = dat[,X.name] - dat[,Q.name]
    Delta = dat[, event.name]
    
    u = as.numeric(colnames(Fuz.mx))
    v = as.numeric(colnames(Gvz.mx))
    d = as.numeric(colnames(Sdz.mx))
    
    tau.Tmax = max(u) + 1e-5
    tau.Dmax = max(d) + 1e-5
    
    # compute the IPCW weights
    Sdyz.vec = diag(1 - CDF_eval.mx(X.res, 1-Sdz.mx))   # S_D(X-Q|Q,A,Z)
    id0 = (dat[, event.name] == 0)
    w_D = rep(0, nrow(dat))
    w_D[id0] = 1/pmax(Sdyz.vec[id0], trim)
    
    ## compute aug_1
    # compute aug_11
    Fxz.vec = diag(CDF_eval.mx(X, Fuz.mx))  # F(X|A,Z)
    Guz.mx = CDF_eval.mx(u, Gvz.mx)    # G(u|A,Z).mx
    
    
    nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)  # row - subject; col - u
    
    ind_Xu.mx = matrix(nrow = nn, ncol = length(u))     #  I(X <= u)
    for(i in 1:nn){
        ind_Xu.mx[i,] = as.numeric(X[i] <= u)
    }  
    fu1.mx = ind_Xu.mx * nuu.mx / pmax(Guz.mx, trim)
    aug_11 = w_D * as.vector(int_fmx_dF(tau.Tmax, fu1.mx, Fuz.mx)) / pmax(1-Fxz.vec, trim)   # this is faster than diag(int_infty_fmx_dF(X,...))
    
    fu1_const1.mx = ind_Xu.mx / pmax(Guz.mx, trim)
    aug_11_const1 = w_D * as.vector(int_fmx_dF(tau.Tmax, fu1_const1.mx, Fuz.mx)) / pmax(1-Fxz.vec, trim)   # this is faster than diag(int_infty_fmx_dF(X,...))
    
    
    # compute aug_12
    # start_time <- Sys.time()
    # system.time({
    wnu.mx = nuu.mx / pmax(Guz.mx, trim)
    int_wnu_dF = matrix(nrow = nn, ncol = length(d))  # int \ind(Q+d<=u)\nu(u)/G(u|A,Z) dF(u|A,Z)
    F_Qu_z.mx = matrix(nrow = nn, ncol = length(d)) # F(Q+u|A,Z)
    
    wnu_const1.mx = 1 / pmax(Guz.mx, trim)
    int_wnu_const1_dF = matrix(nrow = nn, ncol = length(d))  # int \ind(Q+d<=u) 1/G(u|A,Z) dF(u|A,Z)

    for(i in 1:nn){
        wnu.mx.i = matrix(wnu.mx[i,], nrow=1)
        wnu_const1.mx.i = matrix(wnu_const1.mx[i,], nrow=1)
        Fuz.mx.i = matrix(Fuz.mx[i,], nrow=1)
        colnames(Fuz.mx.i) = colnames(Fuz.mx)
        
        int_wnu_dF[i,] = int_infty_fmx_dF(Q[i]+d, wnu.mx.i, Fuz.mx.i)  
        int_wnu_const1_dF[i,] = int_infty_fmx_dF(Q[i]+d, wnu_const1.mx.i, Fuz.mx.i)  
        
        F_Qu_z.mx[i,] = CDF_eval.mx(Q[i]+d, Fuz.mx.i)
    }
    # })
    # now_time <- Sys.time()
    # now_time - start_time      # ~ 10.66324 secs
    
    FDdz.mx = 1-Sdz.mx
    ind_dXQ.mx = matrix(nrow = nn, ncol = length(d))
    for(i in 1:nn){
        ind_dXQ.mx[i,] = as.numeric(d <= X[i]-Q[i])
    }
    
    fd2.mx = ind_dXQ.mx * int_wnu_dF/ (pmax(1-F_Qu_z.mx, trim) * (pmax(Sdz.mx, trim))^2)
    aug_12 = as.vector(int_fmx_dF(tau.Dmax, fd2.mx, FDdz.mx))
    
    fd2_const1.mx = ind_dXQ.mx * int_wnu_const1_dF/ (pmax(1-F_Qu_z.mx, trim) * (pmax(Sdz.mx, trim))^2)
    aug_12_const1 = as.vector(int_fmx_dF(tau.Dmax, fd2_const1.mx, FDdz.mx))
    
    aug_1 = aug_11 - aug_12
    
    aug_1_const1 = aug_11_const1 - aug_12_const1
    
    
    return(list(aug_1 = aug_1, # aug_11 = aug_11, aug_12 = aug_12,
                aug_1_const1 = aug_1_const1))
}




aug_QD_2 <- function(dat, nu, Fuz.mx, Gvz.mx, Sdyz.vec,
                     Q.name, event.name, trim = 1e-7){
    
    nn = nrow(dat)
    Q = dat[,Q.name]
    Delta = dat[, event.name]
    
    u = as.numeric(colnames(Fuz.mx))
    nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)  # row - subject; col - u
    
    mqz = diag(int_fmx_dF(Q, nuu.mx, Fuz.mx))
    Gqz = CDF_eval(Q, Gvz.mx)
    Fqz = CDF_eval(Q, Fuz.mx)
    
    # aug_21 = (1-Delta)/pmax(Sdyz.vec, trim) * mqz / (pmax(1-Fqz,trim) * pmax(Gqz,trim))  # for testing purposes
    aug_2 = (1 - Delta/pmax(Sdyz.vec, trim)) * mqz / (pmax(1-Fqz,trim) * pmax(Gqz,trim))  
    
    aug_2_const1 = (1 - Delta/pmax(Sdyz.vec, trim)) * Fqz / (pmax(1-Fqz,trim) * pmax(Gqz,trim))
    
    
    return(list(aug_2 = aug_2, 
                aug_2_const1 = aug_2_const1))
}




aug_QD_3 <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                     X.name, Q.name, event.name, trim = 1e-7){
    
    nn = nrow(dat)
    X = dat[,X.name]
    Q = dat[,Q.name]
    X.res = dat[,X.name] - dat[,Q.name]
    Delta = dat[, event.name]
    
    u = as.numeric(colnames(Fuz.mx))
    v = as.numeric(colnames(Gvz.mx))
    d = as.numeric(colnames(Sdz.mx))
    
    nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)  # row - subject; col - u
    
    tau.Tmax = max(u) + 1e-5
    tau.Qmax = max(v) + 1e-5
    tau.Dmax = max(d) + 1e-5
    
    # compute the IPCW weights
    Sdyz.vec = diag(1 - CDF_eval.mx(X.res, 1-Sdz.mx))   # S_D(X-Q|Q,A,Z)
    id0 = (dat[, event.name] == 0)
    w_D = rep(0, nrow(dat))
    w_D[id0] = 1/pmax(Sdyz.vec[id0], trim)
    
    # compute aug_31
    mvz.mx = int_fmx_dF(v, nuu.mx, Fuz.mx)
    Fvz.mx = CDF_eval.mx(v, Fuz.mx)
    
    F_Xv_z.mx = matrix(nrow = nn, ncol = length(v))  # used in the computation of aug_31
    ind_Qv.mx = matrix(nrow = nn, ncol = length(v))  # used in both aug_31 and aug_32
    for(i in 1:nn){
        Fuz.mx.i = matrix(Fuz.mx[i,], nrow=1)
        colnames(Fuz.mx.i) = colnames(Fuz.mx)
        
        F_Xv_z.mx[i,] = CDF_eval.mx(pmin(v, X[i]), Fuz.mx.i)
        
        ind_Qv.mx[i,] = (Q[i] <= v)
    }
    
    fv31.mx = mvz.mx * ind_Qv.mx / ( pmax(1-F_Xv_z.mx, trim) * (pmax(Gvz.mx, trim))^2 )
    aug_31 = w_D * int_fmx_dF(tau.Qmax, fv31.mx, Gvz.mx)
    
    fv31_const1.mx = Fvz.mx * ind_Qv.mx / ( pmax(1-F_Xv_z.mx, trim) * (pmax(Gvz.mx, trim))^2 )
    aug_31_const1 = w_D * int_fmx_dF(tau.Qmax, fv31_const1.mx, Gvz.mx)
    
    
    ## Compute aug_32
    # start_time <- Sys.time()
    # system.time({
    aug32_int_dG = matrix(nrow = nn, ncol = length(d))   # int ... dG(v|A,Z) in aug_32
    aug32_int_dG_const1 = matrix(nrow = nn, ncol = length(d))   # int F(v|A,Z)/... dG(v|A,Z) in aug_32_const1
    v.vec = rep(v, length(d))
    d.vec = rep(d, each = length(v))
    for(i in 1:nn){
        Fuz.mx.i = matrix(Fuz.mx[i,], nrow=1)
        colnames(Fuz.mx.i) = colnames(Fuz.mx)
        
        F_Qdv_z.mx = matrix(CDF_eval.mx(pmin(Q[i]+d.vec, v.vec), Fuz.mx.i), ncol = length(v), byrow = TRUE)
        
        fv32_i.mx = t(replicate(length(d), mvz.mx[i,]*ind_Qv.mx[i,])) / ( pmax(1-F_Qdv_z.mx, trim) * t( replicate(length(d), (pmax(Gvz.mx[i,], trim))^2) ) )
        aug32_int_dG[i,] = int_fmx_dF(tau.Qmax, fv32_i.mx, t(replicate(length(d), Gvz.mx[i,])) )
        
        fv32_i_const1.mx = t(replicate(length(d), Fvz.mx[i,]*ind_Qv.mx[i,])) / ( pmax(1-F_Qdv_z.mx, trim) * t( replicate(length(d), (pmax(Gvz.mx[i,], trim))^2) ) )
        aug32_int_dG_const1[i,] = int_fmx_dF(tau.Qmax, fv32_i_const1.mx, t(replicate(length(d), Gvz.mx[i,])) )
    }
    # })
    # now_time <- Sys.time()
    # now_time - start_time   # ~ 3min
    
    atrisk_XQ.mx = matrix(nrow = nn, ncol = length(d))
    for(i in 1:nn){
        atrisk_XQ.mx[i,] = (X[i]-Q[i] >= d)
    }
    
    FDdz.mx = 1-Sdz.mx
    fd32.mx = aug32_int_dG * atrisk_XQ.mx / ((pmax(Sdz.mx, trim))^2)
    aug_32 = int_fmx_dF(tau.Dmax, fd32.mx, FDdz.mx)
    
    fd32_const1.mx = aug32_int_dG_const1 * atrisk_XQ.mx / ((pmax(Sdz.mx, trim))^2)
    aug_32_const1 = int_fmx_dF(tau.Dmax, fd32_const1.mx, FDdz.mx)
    
    
    aug_3 = aug_31 - aug_32
    
    aug_3_const1 = aug_31_const1 - aug_32_const1
    
    return(list(aug_3 = aug_3, # aug_31 = aug_31, aug_32 = aug_32,
                aug_3_const1 = aug_3_const1))
}




