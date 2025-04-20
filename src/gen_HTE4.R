## Simulate T^*(a) from AFT model with weibull error - so T^*|Z follows a Cox model
## Complex tau()

## simulate data with two covariates
# d number of covariates generated
gen <- function(n, multi = 100, gamma.A, T.min, Q.min, D.min, Q.max, tau, 
                beta.Q, beta.C, shape.C,
                beta.T, err.dist = c("weibull", "lognormal"),
                err.sigma = NULL, err.bdd = Inf,
                err.shape = NULL, err.scale = NULL,
                d = 2){
    
    if(d<2){
        stop("This function only implement the case with more than 2 covariates.")
    }
    
    Z1 = runif(multi*n, -1, 1)
    Z2 = runif(multi*n, -1, 1)
    
    Z_rest = matrix(runif(multi*n*(d-2), -1, 1), nrow = multi*n, ncol = d-2)
    if(d>2){
        colnames(Z_rest) = paste("Z", (3:d), sep = "")
    }
    
    linpred = cbind(1,Z1,Z2) %*% gamma.A
    prob.A = exp(linpred)/(1+exp(linpred))
    A = rbinom(multi*n, size = 1, prob = prob.A)
    
    Q2.max = tau - Q.min
    Q2.min = tau - Q.max
    U2 = runif(multi*n, min = 0, max = 1)
    Q2 = Q2.min + as.vector((Q2.max-Q2.min)*(1-U2^(exp(-cbind(1, A, Z1, Z2, A*Z1, A*Z2) %*% beta.Q))))
    QQ = tau - Q2
    
    # Simulate T from AFT model - complex tau()
    if(err.dist == "lognormal"){
        eps1 = pmax(pmin(rnorm(multi*n, mean = 0, sd = err.sigma), err.bdd), -err.bdd)
        eps0 = pmax(pmin(rnorm(multi*n, mean = 0, sd = err.sigma), err.bdd), -err.bdd)
    }else if(err.dist == "weibull"){
        eps1 = pmax(pmin(rweibull(multi*n, shape = err.shape, scale = err.scale) - err.scale*gamma(1+1/err.shape),
                         err.bdd), -err.bdd)
        eps0 = pmax(pmin(rweibull(multi*n, shape = err.shape, scale = err.scale) - err.scale*gamma(1+1/err.shape), 
                         err.bdd), -err.bdd)
    }else{
        stop("This setting is not impletmented yet")
    }
    
    TT1 = exp(cbind(1,1,sin(pi*Z1), sqrt(abs(Z2))) %*% beta.T + eps1)
    TT0 = exp(cbind(1,0,0,rep(0, length(Z1))) %*% beta.T + eps0)
    TT = TT0
    TT[A==1] = TT1[A==1]
    # exp(0.2*(cbind(1,A,A*Z1,sqrt(abs(Z2))) %*% beta.T + eps))
    
    
    D = D.min + rweibull(multi*n, shape = shape.C, scale = exp(-1/shape.C * (cbind(1, QQ, A, Z1, Z2, A*Z1, A*Z2)%*%beta.C)))
    # D = D.min + rweibull(multi*n, shape = shape.C, scale = exp(-1/shape.C * (cbind(1, A, Z1)%*%beta.C)))
    C = QQ + D
    
    X = pmin(TT, C)
    delta = as.numeric(TT<C)
    
    
    dat.full = data.frame(time.1 = TT1, time.0 = TT0, time = TT, C = C, 
                          Q = QQ, X = X, delta = delta, A = A, Z1 = Z1, Z2 = Z2, Z_rest)
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



