
## function that compute different parts involved in the estimating function (EF)
# Fuz.mx: the matrix representing the CDF of T|(A,Z) for each subject in dat
# Gvz.mx: the matrix representing the CDF of Q|(A,Z) for each subject in dat
truncAIPW_transMean_EF <- function(dat, nu, Fuz.mx, Gvz.mx,
                                   T.name, Q.name, trim = 0){
    
    if(nrow(dat) != nrow(Fuz.mx) | nrow(dat) != nrow(Gvz.mx)){
        stop("The number of rows of dat, Fuz.mx and Gvz.mx are not the same. ")
    }
    
    time = dat[, T.name]
    Q = dat[, Q.name]
    
    u = as.numeric(colnames(Fuz.mx))  # jumps.T
    v = as.numeric(colnames(Gvz.mx))  # jumps.Q
    
    tau2 = max(v)+1e-10
    
    Gtz = CDF_eval(time, Gvz.mx)
    Gqz = CDF_eval(Q, Gvz.mx)
    Fqz = CDF_eval(Q, Fuz.mx)
    
    DDen1 = 1/pmax(Gtz, trim)
    DDen2 = Fqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))
    
    Fvz.mx = CDF_eval.mx(v, Fuz.mx)
    
    nn = nrow(dat)
    atrisk.mx = matrix(nrow = nn, ncol = length(v))
    for(i in 1:nn){
        atrisk.mx[i,] = (Q[i] <= v & v < time[i])    # at risk indicator for subject i at times in jumps.Q
    }
    f.mx = atrisk.mx*Fvz.mx/(pmax(Gvz.mx, trim)^2*pmax(1-Fvz.mx, trim))
    DDen3 = as.vector(int_fmx_dF(tau2, f.mx, Gvz.mx))
    
    
    NNum1 = nu(time)/pmax(Gtz, trim)
    
    nuu.mx = matrix(rep(nu(u), nrow(dat)), nrow = nrow(dat), byrow = TRUE)
    mqz = diag(int_fmx_dF(Q, nuu.mx, Fuz.mx))
    NNum2 = mqz/(pmax(Gqz, trim)*pmax(1-Fqz, trim))
    
    mvz.mx = int_fmx_dF(v, nuu.mx, Fuz.mx)
    fnu.mx = atrisk.mx*mvz.mx/(pmax(Gvz.mx, trim)^2*pmax(1-Fvz.mx, trim))
    NNum3 = as.vector(int_fmx_dF(tau2, fnu.mx, Gvz.mx))
    
    Num_AIPW = NNum1 + NNum2 - NNum3
    Den_AIPW = DDen1 + DDen2 - DDen3
    
    
    # For the estimating equation of the IPW and Reg.T1, Reg.T2 estimators
    DDen4 = 1/pmax(1-Fqz, trim)
    NNum4 = nu(time) + mqz/pmax(1-Fqz, trim)
    
    tau.Tmax = max(u) + 1
    Enutz = as.vector(int_fmx_dF(tau.Tmax, nuu.mx, Fuz.mx))
    NNum5 = Enutz/pmax(1-Fqz, trim)
    
    
    
    return(list(Num_AIPW = Num_AIPW, Den_AIPW = Den_AIPW,
                Num_IPW.Q = NNum1, Den_IPW.Q = DDen1,
                Num_Reg.T1 = NNum4, Num_Reg.T2 = NNum5, Den_Reg = DDen4))
    
}


