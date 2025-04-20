

# Return applying truncC_AIPW to the function \nu(T) and the constant function 1
#### ! Need testing
#### Double check plus or minus augmentation term - should be correct, but be aware if the results are not as expected.
truncC_AIPW_transMean <- function(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx,
                                  X.name, Q.name, event.name, trim = 1e-7){
    
    X.res = dat[,X.name] - dat[,Q.name]
    
    # apply truncAIPW
    efs = truncAIPW_transMean_EF(dat, nu, Fuz.mx, Gvz.mx, X.name, Q.name, trim)
    Num_AIPW = efs$Num_AIPW
    Den_AIPW = efs$Den_AIPW
    
    Num_IPW.Q = efs$Num_IPW.Q
    Den_IPW.Q = efs$Den_IPW.Q
    
    # compute the IPCW weights
    Sdyz.vec = diag(1 - CDF_eval.mx(X.res, 1-Sdz.mx))
    id1 = (dat[, event.name] == 1)
    w_IPCW = rep(0, nrow(dat))
    w_IPCW[id1] = 1/pmax(Sdyz.vec[id1], trim)
    
    result_Aug_QD = aug_QD(dat, nu, Fuz.mx, Gvz.mx, Sdz.mx, Sdyz.vec,
                           X.name, Q.name, event.name, trim)
    Aug_QD_nu = result_Aug_QD$Aug_QD_nu
    Aug_QD_const1 = result_Aug_QD$Aug_QD_const1
    
    truncC_nu = w_IPCW * Num_AIPW + Aug_QD_nu 
    truncC_const1 = w_IPCW * Den_AIPW + Aug_QD_const1
    
    
    return(list(truncC_nu = truncC_nu, 
                truncC_const1 = truncC_const1))
}