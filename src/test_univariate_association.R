## functions for testing the univariate association 

get_km_test <- function(cov.name, cov.name.cat, data, Q.name, T.name, status.name,
                        isplot = FALSE){
    
    data = data[,c(cov.name, cov.name.cat, Q.name, T.name, status.name)]
    data = data[complete.cases(data), ]
    # print(paste("n =", nrow(data)))
    
    formula = paste("Surv(", Q.name, ",", T.name, ",", status.name, ")~", cov.name) %>% as.formula()
    cox_fit = coxph(formula, data = data)
    p = summary(cox_fit)$sctest[3]   # p-value from the score test
    
    formula_cat = paste("Surv(", Q.name, ",", T.name, ",", status.name, ")~", cov.name.cat) %>% as.formula()
    cox_fit.cat = coxph(formula_cat, data = data)
    p.cat = summary(cox_fit.cat)$sctest[3]   # p-value from the score test
    
    
    fit = surv_fit(formula_cat, data = data)  # surv_fit() [in survminer package], which is an extension to the R base function survfit() with more functionalities -  allowing decompose the formula in survfit()  #dat[,status.name]
    
    if(isplot){
        xlim = c(min(data[,T.name][data[,status.name]==1])-0.1, max(data[,T.name]))
        ggsurvplot(fit, conf.int = F, data = data,
                   pval = round(p, 4), # show p-value of log-rank test
                   pval.coord = c(xlim[1]+2, 0.1), # p-value location
                   pval.size = 5, # adjust p-value size
                   xlab = "Age", xlim = xlim, break.x.by = 2) %>% print()
    }
    
    return(c(p = p, p.cat = p.cat))
}


get_km_test_Q <- function(cov.name, cov.name.cat, data, Q.name, T.name,
                          isplot = FALSE){
    
    data = data[,c(cov.name, cov.name.cat, Q.name, T.name)]
    data = data[complete.cases(data), ]
    # print(paste("n =", nrow(data)))
    
    tau = max(data[,T.name])+0.1
    data$Q2 = tau - data[,Q.name]
    data$T2 = tau - data[,T.name]
    data$delta.1 = rep(1, nrow(data))
    
    formula = paste("Surv(T2, Q2, delta.1)~", cov.name) %>% as.formula()
    cox_fit = coxph(formula, data = data)
    p = summary(cox_fit)$sctest[3]   # p-value from the score test
    
    formula_cat = paste("Surv(T2, Q2, delta.1)~", cov.name.cat) %>% as.formula()
    cox_fit.cat = coxph(formula_cat, data = data)
    p.cat = summary(cox_fit.cat)$sctest[3]   # p-value from the score test
    
    fit = surv_fit(formula_cat, data = data)  # surv_fit() [in survminer package], which is an extension to the R base function survfit() with more functionalities -  allowing decompose the formula in survfit()  #dat[,status.name]
    
    
    if(isplot){
        xlim = c(min(fit$time)-0.1, max(fit$time))
        ggsurvplot(fit, conf.int = F, data = data,
                   pval = round(p, 4), # show p-value of log-rank test
                   pval.coord = c(xlim[1]+2, 0.1), # p-value location
                   pval.size = 5, # adjust p-value size
                   xlab = paste("Reversed Age (", round(tau,1), " - Age)", sep = ""), 
                   xlim = xlim, break.x.by = 2) %>% print()
    }
    
    return(c(p = p, p.cat = p.cat))
    
}


get_km_test_D <- function(cov.name, cov.name.cat, data, Q.name, T.name, status.name,
                          isplot = FALSE){
    
    data = data[,c(cov.name, cov.name.cat, Q.name, T.name, status.name)]
    data = data[complete.cases(data), ]
    # print(paste("n =", nrow(data)))
    
    if(sum("Y" == c(Q.name, T.name, status.name, cov.name, cov.name.cat))>0|
       sum("delta.D" == c(Q.name, T.name, status.name, cov.name, cov.name.cat))>0){
        stop("The variable names 'Y' and 'delta.D' are used in the intemediate computation, so they cannot be used as the inputs. ")
    }
    
    data$Y = data[,T.name] - data[,Q.name] 
    data$delta.D = 1-data[,status.name] 
    formula = paste("Surv(Y,delta.D)~", cov.name) %>% as.formula()
    cox_fit = coxph(formula, data = data)
    p = summary(cox_fit)$sctest[3]   # p-value from the score test
    
    formula_cat = paste("Surv(Y,delta.D)~", cov.name.cat) %>% as.formula()
    cox_fit.cat = coxph(formula_cat, data = data)
    p.cat = summary(cox_fit.cat)$sctest[3]   # p-value from the score test
    
    
    fit = surv_fit(formula_cat, data = data)  # surv_fit() [in survminer package], which is an extension to the R base function survfit() with more functionalities -  allowing decompose the formula in survfit() 
    
    if(isplot){
        xlim = c(min(data$Y)-0.1, max(data$Y))
        # xlim = c(0, max(data[,T.name]))
        ggsurvplot(fit, conf.int = F, data = data,
                   pval = round(p,4), # show p-value of log-rank test
                   pval.coord = c(xlim[1]+0.5, 0.1), # p-value location
                   pval.size = 5, # adjust p-value size
                   xlab = "Residual age", xlim = xlim, break.x.by = 2) %>% print()
        
        # print(paste("p-value from the score test of fitting a Cox model:", p))
    }
    
    return(c(p = p, p.cat = p.cat))
}