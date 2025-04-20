## main function for simulating data from HTE3 setup and estimate CATE
#- Simulate T^* from AFT model
rm(list = ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

start_time <- Sys.time()
system.time({
    
# library(boot)
library(cvTools)  # used for creating folds for cross-fitting
library(survival)
library(randomForestSRC)
library(LTRCforests)
library(glmnet)
library(Matrix)
library(splines)
library(survPen)
library(gbm)
library(xgboost)

source("est_nuisance.R")
source("misc.R")
source("truncAIPW.R")
source("truncAC_AIPW.R")
source("truncC_AIPW.R")
# source("CATE.R")
# source("development.R")
source("gen_HTE3.R")  # simulate T(a) from AFT model - the true tau() with nu = log() is a linear function in Z1
source("trunclearner.R")
source("cvboost.R")
source("cvboost_wsq.R")

load("seeds_input.rda")
load("simu_setting2.rda")  #inputs/simu_setting2.rda


rm(gamma.A, T.min, beta.T, beta.C, beta.Q)
gamma.A = c(0, 1, -1)
beta.T = 0.2 * c(4, 1, -1, 1)
beta.Q = c(-0.6,         0,    0.4,    0,     0,   0.5)
beta.C = c(-1,    0,   0.3,    0.1,  0.2,     0,     0)

# Distribution for the error term
err.sigma = 0.2*0.1
err.bdd = 0.2*3

err.shape = 2
err.scale = 0.04

# S_D
err.scale * sqrt(gamma(1+2/err.shape) - (gamma(1+1/err.shape))^2)
err.dist  = "weibull"
D.min = 0

# For tuning of xgboost
num_search_rounds = 50
ntrees_max = 1000

# Sample size
n = 500

est_approach_G = "IPCW"   #  "IPCW" or "truncIPW.F"
est_approach_PS = "truncIPCW"   # "truncIPCW" or "truncIPW.F"

############## prepare the inputs ###################
K = 5   # number of folds for cross-fitting for the cf-estimator
trim = 0.05   
trim.est = 0


t00 = 6
nu <- function(t,t0=t00){
    # indicator function
    # result = as.numeric(t>t0)
    
    # result = pmin(t, t0)
    
    result = log(t)
    
    return(result)
}


### Inputs used for testing the truncrlearner() function step by step
cov.names.CATE = c("Z1","Z2")


################ simulate data ##################
set.seed(seed)
genresult = gen(n, multi, gamma.A, T.min, Q.min, D.min, Q.max, tau,
                beta.Q, beta.C, shape.C,
                beta.T, err.dist = err.dist,
                err.bdd = err.bdd,
                err.shape = err.shape, err.scale = err.scale)     # for HTE3 setup
dat.full = genresult$dat.full

dat = genresult$dat
dat$delta.1 = rep(1, nrow(dat))
dat$AZ1 = dat$A * dat$Z1
dat$AZ2 = dat$A * dat$Z2
# dat$Z1sqrt = sqrt(abs(dat$Z1))
dat$Z2sqrt = sqrt(abs(dat$Z2))

CATE.truth = cbind(1, dat$Z1) %*% c(beta.T[2], beta.T[3])

############### oracle-truncR-learner #################

X.name = "X"
Q.name = "Q"
A.name = "A"
event.name = "delta"
event.name.T = "delta"

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

jumps.X = sort(dat[,X.name])
jumps.Q = sort(dat[,Q.name])
jumps.Y = sort(dat[,X.name] - dat[,Q.name])

tau1 = min(jumps.X)
tau2 = max(jumps.Q)
tau3 = max(jumps.Y)

tau.max = max(jumps.X)



## Compute the true nuisance parameters

cov.names.T.oracle = c("A","AZ1","Z2sqrt")
u = c(min(jumps.X)-1e-10, jumps.X, max(jumps.X)+1e-10)
Fuz.mx = truth_F.AFT(dat, u, 
                     cov.names = cov.names.T.oracle, beta.T, 
                     err.dist = err.dist, err.sigma = err.sigma, 
                     err.shape = err.shape, err.scale = err.scale)
Fuz_A1.mx = truth_F.AFT(dat_A1, u, 
                        cov.names = cov.names.T.oracle, beta.T, 
                        err.dist = err.dist, err.sigma = err.sigma, 
                        err.shape = err.shape, err.scale = err.scale)
Fuz_A0.mx = truth_F.AFT(dat_A0, u, 
                        cov.names = cov.names.T.oracle, beta.T, 
                        err.dist = err.dist, err.sigma = err.sigma, 
                        err.shape = err.shape, err.scale = err.scale)

d = c(min(jumps.Y)-1e-10, jumps.Y, max(jumps.Y)+1e-10)
cov.names.C.oracle = c("Q","A","Z1","Z2","AZ1","AZ2")
Sdz.mx = truth_Sd(dat, d, 
                  cov.names = cov.names.C.oracle, D.min, beta.C, shape.C)

v = c(tau1-1e-10, jumps.Q[jumps.Q>=tau1], max(jumps.Q)+1e-10)
cov.names.Q.oracle = c("A","Z1","Z2","AZ1","AZ2")
Gvz.mx = truth_G(dat, v, 
                 cov.names = cov.names.Q.oracle, Q.min, Q.max, tau, beta.Q)

cov.names.A.oracle = c("Z1","Z2")
PS = truth_PS(dat, cov.names.A.oracle, gamma.A)


trunclearner_result_oracle = trunclearner(dat, nu, cov.names.CATE, cov.names.CATE.binary = NULL, K,
                                          X.name, Q.name, event.name, A.name, trim,
                                          Fuz.mx = Fuz.mx, Fuz_A1.mx = Fuz_A1.mx, Fuz_A0.mx = Fuz_A0.mx,
                                          Gvz.mx = Gvz.mx, Sdz.mx = Sdz.mx, PS = PS,
                                          num_search_rounds=num_search_rounds, ntrees_max=ntrees_max)




########### truncR-learner with estimated nuisance parameters ############

CATE.alg = "xgb"

# models - note that "pCox" needs sample size ~400 for the algorithm to converge with the current setup
model.T = "pCox"  # "pCox", "survPen", "RF"
model.Q = "pCox"
model.A = "gbm"
model.D = "pCox"  # "pCox", "survPen", "RF"

X.name = "X"
Q.name = "Q"
A.name = "A"
event.name = "delta"
event.name.T = "delta"

cov.names.T = c("Z1","Z2")
cov.names.Q = c("Z1","Z2")
cov.names.A = c("Z1","Z2")
cov.names.D = c("Q","Z1","Z2")
cov.names.binary.T = c("A")
cov.names.binary.Q = c("A")
cov.names.binary.A = NULL
cov.names.binary.D = c("A")
cov.names.CATE = c("Z1","Z2")



formula_survPen.T = ~ smf(X, df = 5) + smf(Z1, df = 5) + smf(Z2, df = 5) + A # + tint(Z1, by=A, df = 5)
formula_survPen.D = ~ smf(Y, df = 5) + smf(Z1, df = 5) + smf(Z2, df = 5) + A

options.F = list(trim = trim.est,
                 ntree = 500, mtry = 2,
                 df = 7, nfolds = 10, s = "lambda.1se", alpha = 0.5,
                 formula.survPen = formula_survPen.T)
options.Sd = list(trim = trim.est,
                  ntree = 500, mtry = 2,
                  df = 7, nfolds = 10, s = "lambda.1se", alpha = 0.5,
                  formula.survPen = formula_survPen.D,
                  nfolds.OOF = 5)
options.G = list(trim = trim.est,
                 df = 7, nfolds = 10, s = "lambda.1se", alpha = 0.5,
                 trunc.weighting = TRUE)
options.PS = list(trim = trim.est,
                  df = 7,
                  ntree = 500)

options.CATE.xgb = list(booster = "gbtree",
                        nrounds = 100)


trunclearner_result = trunclearner(dat, nu, cov.names.CATE, cov.names.CATE.binary = NULL, K,
                                   X.name, Q.name, event.name, A.name, trim = trim,
                                   model.T, model.Q, model.A, model.D, 
                                   cov.names.T, cov.names.Q, cov.names.A, cov.names.D, 
                                   cov.names.binary.T, cov.names.binary.Q, cov.names.binary.A, cov.names.binary.D,
                                   options.F, options.G, 
                                   options.PS, options.Sd,
                                   est_approach_G = est_approach_G, est_approach_PS,
                                   num_search_rounds=num_search_rounds, ntrees_max=ntrees_max)



})
now_time <- Sys.time()
time_elapsed = now_time - start_time  
time_elapsed



################# Save the results ##################

save(trunclearner_result, trunclearner_result_oracle,
     CATE.truth,
     dat, beta.T,
     model.T, model.Q, model.D, model.A, 
     options.F, options.G, options.PS, options.Sd,
     trim, trim.est, 
     time_elapsed,
     file = "result.rda")



