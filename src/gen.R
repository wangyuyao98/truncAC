gen <- function(n, multi = 100, gamma, T.min, Q.min, D.min, Q.max, tau, 
                beta.Q, beta.T, beta.C, shape.T, shape.C){
    
    Z1 = runif(multi*n, -1, 1)
    linpred = cbind(1,Z1) %*% gamma
    prob.A = exp(linpred)/(1+exp(linpred))
    A = rbinom(multi*n, size = 1, prob = prob.A)
    
    Q2.max = tau - Q.min
    Q2.min = tau - Q.max
    U2 = runif(multi*n, min = 0, max = 1)
    Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-cbind(1,A,Z1) %*% beta.Q))))
    QQ = tau - Q2
    
    TT1 = T.min + rweibull(multi*n, shape = shape.T, scale = exp(-1/shape.T * (cbind(1, 1, Z1)%*%beta.T)))
    TT0 = T.min + rweibull(multi*n, shape = shape.T, scale = exp(-1/shape.T * (cbind(1, 0, Z1)%*%beta.T)))
    TT = TT0
    TT[A==1] = TT1[A==1]
    
    D = D.min + rweibull(multi*n, shape = shape.C, scale = exp(-1/shape.C * (cbind(1, QQ, A, Z1)%*%beta.C)))
    # D = D.min + rweibull(multi*n, shape = shape.C, scale = exp(-1/shape.C * (cbind(1, A, Z1)%*%beta.C)))
    C = QQ + D
    
    X = pmin(TT, C)
    delta = as.numeric(TT<C)
    
    
    dat.full = data.frame(time.1 = TT1, time.0 = TT0, time = TT, C = C, 
                          Q = QQ, X = X, delta = delta, A = A, Z1 = Z1)
    dat.all = dat.full
    obs.id = which(dat.full$Q < dat.full$time)
    dat.obs = dat.full[obs.id,]
    if(length(obs.id)<n){
        stop('Truncation rate is high. Need to increase the mutiplier.')
    }
    dat = dat.obs[1:n, -(1:4)]
    dat.full = dat.full[(1:obs.id[n]), ]
    
    return(list(dat = dat, dat.full = dat.full, dat.all = dat.all))
}