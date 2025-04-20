# ### organize OSG results
# 
# ### Run on OSG -------------------------------------------------------------------
# rm(list = ls())
# 
# n.boot = 200  # total number of bootstrap runs
# nb.each = 5     # nb.each = 20 for dr-Param;  nb.each = 5 for cf-nonparam
# ns = ceiling(n.boot/nb.each)
# itern = 500   # number of simulation runs
# njobs = itern*ns
# 
# est.list = list()
# bootresult.list = list()
# 
# missing_index <- NULL
# missing_sample_index <- NULL
# for(i in 1:itern){
#     bresult = matrix(nrow = nb.each*ns, ncol = 10)
#     # bresult = matrix(nrow = nb.each*ns, ncol = 8)  # old results that does not contain SE for ATE
# 
#     for(j in ns:1){
#         id = (i-1)*ns + j
#         missing_ij = tryCatch(load(paste("result_", id, ".rda", sep = "")),
#                               error = function(cond){NA})
#         if(is.na(missing_ij)){
#             bresult[((j-1)*nb.each+1):(j*nb.each),] = NA
#             missing_index = cbind(missing_index, id)
#             missing_sample_index = cbind(missing_sample_index, c(i,j))
#         }else{
#             load(paste("result_", id, ".rda", sep = ""))
#             bresult[((j-1)*nb.each+1):(j*nb.each),] = bootresult$t
#         }
#     }
#     bootresult.list[[i]] = bresult
#     est.list[[i]] = tryCatch(bootresult$t0, error = function(cond){NA})
# }
# 
# missing_index
# sort(unique(missing_sample_index[1,]))
# 
# 
# # save the combined results
# save(est.list,
#      bootresult.list,
#      missing_index,
#      missing_sample_index,
#      itern, n.boot, nb.each,
#      file = paste("cmb_bootresults_b", nb.each*ns, ".rda", sep = ""))
# 
# nb.each*ns
# 
# 
# # Save the missing indeces
# write(missing_index,
#       file = paste("indices_missing.txt", sep = ""),
#       sep = "\n")








### Run on my local computer -----------------------------------------------------------

rm(list = ls())
setwd("~/Documents/research/LeftTruncation/github/left_trunc_causal_C/")

# load("OSG_results/ATE/dr/simu2/n1000_FCox2_GCox1_PSlgs2_SDCox1/cmb_bootresults_b200.rda")
# load("OSG_results/ATE/dr/simu2/n1000_FCox1_GCox1_PSlgs1_SDCox1_trimest0/cmb_bootresults_b200.rda")
# load("OSG_results/ATE/cf/simu2/n1000_pCox_pCox_gbm_pCox_df7/cmb_bootresults_b200.rda")
# load("OSG_results/ATE/cf/simu2/n1000_Cox1_Cox1_lgs1_Cox1_trim05_trimest0/cmb_bootresults_b200.rda")
# load("OSG_results/ATE/cf/simu2/n1000_pCox_pCox_gbm_pCox_df7_trimest0/cmb_bootresults_b200.rda")
# load("OSG_results/ATE/dr/simu2/n500_FCox2_GCox1_PSlgs1_SDCox2_trimest0/cmb_bootresults_b200.rda")
load("OSG_results/ATE/cf/simu2/df7_alpha0_lambdamin_trim1/n500_pCox_pCox_gbm_pCox_trimest0/cmb_bootresults_b200.rda")


narm = T
itern = 500 
alpha = 0.05
qz = qnorm(alpha/2, lower = FALSE)

# est.list[[1]]

est_a1 = matrix(nrow = itern, ncol = 2)
est_a0 = matrix(nrow = itern, ncol = 2)
est_ATE = matrix(nrow = itern, ncol = 2)

SE_a1 = matrix(nrow = itern, ncol = 2)
SE_a0 = matrix(nrow = itern, ncol = 2)
SE_ATE = matrix(nrow = itern, ncol = 2)

bootSE_a1 = matrix(nrow = itern, ncol = 2)
bootSE_a0 = matrix(nrow = itern, ncol = 2)
bootSE_ATE = matrix(nrow = itern, ncol = 2)

for(i in 1:itern){
    est_a1[i,] <- est.list[[i]][1:2]
    est_a0[i,] <- est.list[[i]][3:4]
    est_ATE[i,] <- est.list[[i]][1:2] - est.list[[i]][3:4]
    
    SE_a1[i,] <- est.list[[i]][5:6]
    SE_a0[i,] <- est.list[[i]][7:8]
    SE_ATE[i,] <- est.list[[i]][9:10]     # update after updating the code and recollect results from OSG
    
    bootSE_a1[i,] <- apply(bootresult.list[[i]][,1:2], 2, sd)
    bootSE_a0[i,] <- apply(bootresult.list[[i]][,3:4], 2, sd)
    bootSE_ATE[i,] <- apply(bootresult.list[[i]][,1:2] - bootresult.list[[i]][,3:4], 2, sd)
}



# load the truth
load("theta_truth_simu2.rda")
ATE.truth = theta.truth[1] - theta.truth[2]

bias_a1 = apply(est_a1, 2, mean) - theta.truth[1]
bias_a0 = apply(est_a0, 2, mean) - theta.truth[2]
bias_ATE = apply(est_ATE, 2, mean) - ATE.truth

SD_a1 = apply(est_a1, 2, sd)
SD_a0 = apply(est_a0, 2, sd)
SD_ATE = apply(est_ATE, 2, sd)

meanSE_a1 = apply(SE_a1, 2, mean)
meanSE_a0 = apply(SE_a0, 2, mean)
meanSE_ATE = apply(SE_ATE, 2, mean)
meanbootSE_a1 = apply(bootSE_a1, 2, mean, na.rm = narm)
meanbootSE_a0 = apply(bootSE_a0, 2, mean, na.rm = narm)
meanbootSE_ATE = apply(bootSE_ATE, 2, mean, na.rm = narm)

# # median
# meanSE_a1 = apply(SE_a1, 2, median)
# meanSE_a0 = apply(SE_a0, 2, median)
# meanSE_ATE = apply(SE_ATE, 2, median)
# meanbootSE_a1 = apply(bootSE_a1, 2, median, na.rm = narm)
# meanbootSE_a0 = apply(bootSE_a0, 2, median, na.rm = narm)
# meanbootSE_ATE = apply(bootSE_ATE, 2, median, na.rm = narm)

lowerCI_a1 = est_a1 - qz*SE_a1
upperCI_a1 = est_a1 + qz*SE_a1
lowerCI_a0 = est_a0 - qz*SE_a0
upperCI_a0 = est_a0 + qz*SE_a0
lowerCI_ATE = est_ATE - qz*SE_ATE
upperCI_ATE = est_ATE + qz*SE_ATE

lowerbootCI_a1 = est_a1 - qz*bootSE_a1
upperbootCI_a1 = est_a1 + qz*bootSE_a1
lowerbootCI_a0 = est_a0 - qz*bootSE_a0
upperbootCI_a0 = est_a0 + qz*bootSE_a0
lowerbootCI_ATE = est_ATE - qz*bootSE_ATE
upperbootCI_ATE = est_ATE + qz*bootSE_ATE

CP_a1 = apply(lowerCI_a1<theta.truth[1] & theta.truth[1]<upperCI_a1,  2, mean, na.rm = narm)
CP_a0 = apply(lowerCI_a0<theta.truth[2] & theta.truth[2]<upperCI_a0,  2, mean, na.rm = narm)
CP_ATE = apply(lowerCI_ATE<ATE.truth & ATE.truth<upperCI_ATE,  2, mean, na.rm = narm)

bootCP_a1 = apply(lowerbootCI_a1<theta.truth[1] & theta.truth[1]<upperbootCI_a1,  2, mean, na.rm = narm)
bootCP_a0 = apply(lowerbootCI_a0<theta.truth[2] & theta.truth[2]<upperbootCI_a0,  2, mean, na.rm = narm)
bootCP_ATE = apply(lowerbootCI_ATE<ATE.truth & ATE.truth<upperbootCI_ATE,  2, mean, na.rm = narm)



# round of each column of X according to the numbers in digits
format_round <- function(X){
    X0 = cbind(bias = sprintf("%.4f", X[,1]),
               SD = sprintf("%.4f", X[,2]),
               SEbootSE = sprintf("%.4f/%.4f", X[,3], X[,4]),
               CPbootCP = sprintf("%.3f/%.3f", X[,5], X[,6]))
    colnames(X0) <- c("bias", "SD", "SE/bootSE", "CP/bootCP")
    rownames(X0) <- rownames(X)
    
    return(X0)
}




tab_a1 = cbind(bias = bias_a1,
               SD = SD_a1,
               SE = meanSE_a1,
               bootSE = meanbootSE_a1,
               CP = CP_a1,
               bootCP = bootCP_a1)
rownames(tab_a1) <- c("truncAC_AIPW", "IPW")



tab_a0 = cbind(bias = bias_a0,
               SD = SD_a0,
               SE = meanSE_a0,
               bootSE = meanbootSE_a0,
               CP = CP_a0,
               bootCP = bootCP_a0)
rownames(tab_a0) <- c("truncAC_AIPW", "IPW")



tab_ATE = cbind(bias = bias_ATE,
                SD = SD_ATE,
                SE = meanSE_ATE,
                bootSE = meanbootSE_ATE,
                CP = CP_ATE,
                bootCP = bootCP_ATE)
rownames(tab_ATE) <- c("truncAC_AIPW", "IPW")



tab_a1_organized = format_round(tab_a1)
tab_a0_organized = format_round(tab_a0)
tab_ATE_organized = format_round(tab_ATE)

method.names = c("Param", "IPW")
# method.names = c("cf", "IPW")

library(xtable)
print(xtable(cbind(method.names, tab_a1_organized),
             type = "latex", 
             align = rep("c",6),
             file = "tab.tex"),
      include.rownames = FALSE)


print(xtable(cbind(method.names, tab_a0_organized),
             type = "latex", 
             align = rep("c",6),
             file = "tab.tex"),
      include.rownames = FALSE)


print(xtable(cbind(method.names, tab_ATE_organized),
             type = "latex", 
             align = rep("c",6),
             file = "tab.tex"),
      include.rownames = FALSE)


