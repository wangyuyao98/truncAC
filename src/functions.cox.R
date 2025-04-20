## functions used for the code of the DR estimator with coxph() to estimate the nuisance parameters

# # function for computing cumulative hazard
# # haz: the second column is the jump times and the first column is the cumulative hazard
# cumhaz <- function(t, basehaz){
#     id = which((t - basehaz[,'time'])>=0)
#     index = which.min((t - basehaz[,'time'])[id])
#     
#     return(basehaz[index,'hazard'])
# }

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
    cdf = sapply(t, baseCDF.single, basehaz = basehaz.T)
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
            result[,j] = (baseS(tau-t[j], basehaz.Q2))^exp(Z %*% beta.Q2)
        }
    }else{
        result = matrix(nrow = 1, ncol = length(t))
        for(j in 1:length(t)){
            result[,j] = (baseS(tau-t[j], basehaz.Q2))^exp(sum(Z*beta.Q2))
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


## computing integral with respect to F - m(v,Z;F)= \ing_0^v nu(t) dF(t)
# nu: the given function in the parameter of interest
# time: the jump times for F
# Ftz: the corresponding values of F(t|z) at the jump times of T
# return a matrix: each row corresponds to Z_i, each column corresponds to an element in v
integral_F <- function(v, nu, Z, basehaz.T, beta.T){
    n = nrow(Z)
    time = basehaz.T[,'time']
    Ftz = Fz(time, Z, basehaz.T, beta.T)  # each row corresponds to Z_i, col corresponds to jump times of CDF of T
    dF = Ftz - cbind(0, Ftz[,-ncol(Ftz)])  # each row correspond to a subject
    nut = matrix(rep(nu(time),n), nrow = n, byrow = TRUE)
    vals = cbind(0,nut*dF)
    m = length(v)
    
    inte = matrix(nrow = n, ncol = m)
    for(j in 1:m){
        id = findInterval(v[j], c(0, time, Inf))
        inte[,j] = rowSums(cbind(0,vals[,(1:id)]))
    }
    
    return(inte)
}


## computing the integral \int F(v|Z)/((G(v|Z))^2*(1-F(v|Z))) J(v) dG(v|Z)
# Q_orignal, time_original: the Q, T in the sample
integral_JG_F <- function(Q_original, time_original, Z.Q, Z.T, basehaz.Q2, basehaz.T, 
                          beta.Q2, beta.T, tau){
    Q2 = basehaz.Q2[,'time']   # increasing order
    Q = tau - Q2   # decreasing order
    Q = sort(Q)    # change to increasing order
    Fqz = Fz(Q, Z.T, basehaz.T, beta.T)
    Gqz = Gz(Q, Z.Q, basehaz.Q2, beta.Q2, tau)
    dG = Gqz - cbind(0, Gqz[,-ncol(Gqz)])
    n = length(Q_original)
    # m = length(Q)
    
    result = Fqz*dG/((Gqz)^2*(1-Fqz))
    inte = rep(NA, n)
    for(i in 1:n){
        atrisk = (Q_original[i] <= Q & Q <= time_original[i])  # at risk indicator for subject i
        inte[i] = sum(result[i,atrisk])
    }
    
    return(inte)
}


## computing the integral \int m(v,Z;F)/((G(v|Z))^2*(1-F(v|Z))) J(v) dG(v|Z)
# Q_orignal, time_original: the Q, T in the sample
integral_JG_m <- function(Q_original, time_original, Z.Q, Z.T, basehaz.Q2, basehaz.T,beta.Q2, beta.T, tau, nu){
    Q2 = basehaz.Q2[,'time']   # increasing order
    Q = tau - Q2   # decreasing order
    Q = sort(Q)    # change to increasing order
    mqz = integral_F(Q, nu, Z.T, basehaz.T, beta.T)
    Fqz = Fz(Q, Z.T, basehaz.T, beta.T)
    Gqz = Gz(Q, Z.Q, basehaz.Q2, beta.Q2, tau)
    dG = Gqz - cbind(0, Gqz[,-ncol(Gqz)])
    n = length(Q_original)
    
    result = mqz*dG/((Gqz)^2*(1-Fqz))
    inte = rep(NA, n)
    for(i in 1:n){
        atrisk = (Q_original[i] <= Q & Q <= time_original[i])  # at risk indicator for subject i
        inte[i] = sum(result[i,atrisk])
    }
    
    return(inte)
}


# 
# # function for computing the summary table for each column
# summary.mx <- function(X){
#     return(cbind(apply(X, 2, min),
#                  apply(X, 2, quantile, probs = 0.25),
#                  apply(X, 2, quantile, probs = 0.5),
#                  apply(X, 2, mean),
#                  apply(X, 2, quantile, probs = 0.75),
#                  apply(X,2, max)))
# }