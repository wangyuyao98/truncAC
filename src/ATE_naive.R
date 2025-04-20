## naive estimator that ignore all three biases
ATE_naive_boot <- function(dat, id, nu, X.name, A.name){
    
    dat <- dat[id,]
    
    est.a1 = mean(nu(dat[, X.name])[dat[,A.name] == 1]) 
    est.a0 = mean(nu(dat[, X.name])[dat[,A.name] == 0]) 
    
    return(c(naive.a1 = est.a1, naive.a0 = est.a0))
}

