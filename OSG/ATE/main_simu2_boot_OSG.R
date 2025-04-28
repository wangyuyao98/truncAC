# main function for generating data and compute the estimator
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

library(survival)
library(boot)

# Load code under the `src/` folder
source("gen2.R")
source("est_nuisance.R")
source("misc.R")
source("truncAC_AIPW.R")
# Load C++ implementation
library("Rcpp")
library("RcppArmadillo")
Rcpp::sourceCpp("fast_integrals.cpp")


# load("seeds_input.rda") # load seeds
seed = 123
seed.b = 213

## Input parameters under the `inputs` folder.
load("simu_setting2.rda")



## prepare the inputs
n = 500 
n.boot = 1

# models
model.T = "Cox"
model.Q = "Cox"
model.A = "logistic"
model.D = "Cox"

X.name = "X"
Q.name = "Q"
A.name = "A"
event.name.T = "delta"
# event.names.Q = "delta.1"

cov.names.T = c("A","Z1","Z2")
# cov.names.T = c("AZ1", "Z2sq")
cov.names.Q = c("A","Z1","Z2")
# cov.names.Q = c("AZ1", "Z2sq")
cov.names.A = c("Z1","Z2")
# cov.names.A = c("Z1sq","Z2sq")
cov.names.D = c("A","Q","Z1","Z2")
# cov.names.D = c("Q", "AZ1", "Z2sq")

trim = 0.05
trim.est = trim

t00 = 3
nu <- function(t,t0=t00){
    # indicator function
    result = as.numeric(t>t0)
    
    return(result)
}



# simulate data
set.seed(seed)
genresult = gen(n, multi, gamma.A, T.min, Q.min, C.min, Q.max, tau, 
                beta.Q, beta.T, beta.C, shape.T, shape.C)


dat = genresult$dat

dat$delta.1 = rep(1, nrow(dat))
dat$AZ1 = dat$A * dat$Z1
dat$Z1sq = (dat$Z1)^2
dat$Z2sq = (dat$Z2)^2


# # Compute the estimators
# start_time <- Sys.time()
# system.time({
#     
#     result = dr_truncAC_AIPW(dat, nu, X.name, Q.name, event.name.T, A.name, 
#                              model.T, model.Q, model.D, 
#                              cov.names.T, cov.names.Q, cov.names.A, cov.names.D,
#                              trim)
#     
# })
# now_time <- Sys.time()
# now_time - start_time



# estimation and bootstrap
set.seed(seed.b)

start_time <- Sys.time()
system.time({

bootresult = boot(dat, dr_truncAC_AIPW_boot, R = n.boot,
                  nu = nu, X.name = X.name, Q.name = Q.name, event.name = event.name.T, A.name = A.name,
                  model.T = model.T, model.Q = model.Q, model.A = model.A, model.D = model.D,
                  cov.names.T = cov.names.T, cov.names.Q = cov.names.Q, cov.names.A = cov.names.A, cov.names.D = cov.names.D,
                  trim = trim, trim.est = trim.est, simulation = TRUE)

})
now_time <- Sys.time()
now_time - start_time

save(bootresult, trim, n.boot, n, 
     file = "result.rda")










