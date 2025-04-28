# main function for generating data and compute the estimator
rm(list = ls())
# setwd("~/Documents/research/LeftTruncation/github/left_trunc_causal_C")

library(boot)
library(cvTools)  # used for creating folds for cross-fitting
library(survival)
# library(randomForestSRC)
# library(LTRCforests)
library(glmnet)
library(Matrix)
library(splines)
library(survPen)
library(gbm)

# Code under the `src/` folder
source("gen2.R")
source("est_nuisance.R")
source("misc.R")
source("truncAC_AIPW.R")  # Contains both the dr and cf estimators
# Load C++ implementation
library("Rcpp")
library("RcppArmadillo")
Rcpp::sourceCpp("fast_integrals.cpp")


# load("seeds_input.rda") # load seeds
seed = 123
seed.b = 213

## Input parameters under the `inputs` folder.
load("simu_setting2.rda")

n = 500
n.boot = 1  # number of bootstrap

## prepare the inputs 
K = 5   # number of folds for cross-fitting
trim = 0.05   
trim.est = 0
df = 7

# models
model.T = "pCox"  # "Cox", "pCox", "survPen", "RF"
model.Q = "pCox"  # "Cox", "pCox"
model.A = "gbm" # "logistic", "gbm"
model.D = "pCox"  # "Cox", "pCox", "survPen", "RF"

est_approach_G = "IPCW"   # "truncIPW.F", "IPCW"
est_approach_PS = "truncIPCW"   # "truncIPW.F", "truncIPCW"

X.name = "X"
Q.name = "Q"
A.name = "A"
event.name = "delta"

cov.names.T = c("Z1","Z2")
cov.names.Q = c("Z1","Z2")
cov.names.A = c("Z1","Z2")
cov.names.D = c("Q","Z1","Z2")
cov.names.binary.T = c("A")
cov.names.binary.Q = c("A")
cov.names.binary.A = NULL
cov.names.binary.D = c("A")



t00 = 3
nu <- function(t,t0=t00){
    # indicator function
    result = as.numeric(t>t0)
    
    return(result)
}



formula_survPen.T = ~ smf(X, df = 5) + smf(Z1, df = 5) + smf(Z2, df = 5) + A # + tint(Z1, by=A, df = 5)
formula_survPen.D = ~ smf(Y, df = 5) + smf(Z1, df = 5) + smf(Z2, df = 5) + A

options.F = list(trim = trim.est, 
                 ntree = 500, mtry = 2,
                 df = df, nfolds = 10, s = "lambda.min", alpha = 0,
                 formula.survPen = formula_survPen.T)
options.Sd = list(trim = trim.est, 
                  ntree = 500, mtry = 2,
                  df = df, nfolds = 10, s = "lambda.min", alpha = 0,
                  formula.survPen = formula_survPen.D,
                  nfolds.OOF = 5)
options.G = list(trim = trim.est, 
                 df = df, nfolds = 10, s = "lambda.min", alpha = 0,
                 trunc.weighting = FALSE)
options.PS = list(trim = trim.est,
                  df = df,
                  ntree = 500)




### simulate data
# load("/Users/wangyuyao/Documents/research/LeftTruncation/github/left_trunc_causal_C/seeds.rda")
# seed = seeds[1]
# seed.b = seeds[2]

set.seed(seed)
genresult = gen(n, multi, gamma.A, T.min, Q.min, C.min, Q.max, tau, 
                beta.Q, beta.T, beta.C, shape.T, shape.C)
dat = genresult$dat


dat$delta.1 = rep(1, nrow(dat))
dat$AZ1 = dat$A * dat$Z1
dat$Z1sq = (dat$Z1)^2
dat$Z2sq = (dat$Z2)^2


# ### Test line-by-line for the cf_truncAC_AIPW() function ----------------------------
# est_approach_G = "IPCW"
# est_approach_PS = "truncIPCW"
# Fuz.mx = NULL
# Fuz_A1.mx = NULL
# Fuz_A0.mx = NULL
# Gvz.mx = NULL
# Sdz.mx = NULL
# PS = NULL
# simulation = TRUE
# ### End test ------------------------------------------------------------------------




# estimation and bootstrap
set.seed(seed.b)

start_time <- Sys.time()
system.time({

bootresult = boot(dat, cf_truncAC_AIPW_boot, R = n.boot,
                  K = K, nu = nu, X.name = X.name, Q.name = Q.name, event.name = event.name, A.name = A.name,
                  model.T = model.T, model.Q = model.Q, model.A = model.A, model.D = model.D, 
                  cov.names.T = cov.names.T, cov.names.Q = cov.names.Q, 
                  cov.names.A = cov.names.A, cov.names.D = cov.names.D,
                  trim = trim,
                  cov.names.binary.T = cov.names.binary.T, cov.names.binary.Q = cov.names.binary.Q, 
                  cov.names.binary.A = cov.names.binary.A, cov.names.binary.D = cov.names.binary.D,
                  options.F = options.F, options.G = options.G, 
                  options.PS = options.PS, options.Sd = options.Sd,
                  est_approach_G = est_approach_G, est_approach_PS = est_approach_PS,
                  simulation = TRUE)


})
now_time <- Sys.time()
now_time - start_time

# bootresult$t0
# bootresult$t



save(bootresult, file = "result.rda")






