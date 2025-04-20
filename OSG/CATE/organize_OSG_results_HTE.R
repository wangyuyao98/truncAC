### organize OSG results

# Run the following on OSG

# Function to check missing files
check_missing_files <- function(directory) {
  missing_files <- integer()

  # Loop through integers 1 to 500
  for (i in 1:500) {
    file_name <- paste0(directory, "/result.", i, ".rda")
    if (!file.exists(file_name)) {
      missing_files <- c(missing_files, i)
    }
  }

  # Print out the missing files
  if (length(missing_files) > 0) {
    cat("Missing files for the following integers:\n")
    print(missing_files)
  } else {
    cat("No files are missing.\n")
  }

  return(missing_files)
}

# Check for missing files in the current working directory
missing_files = check_missing_files(getwd())
write(missing_files,
      file = paste("indices_missing.txt", sep = ""),
      sep = "\n")




rm(list = ls())
itern = 500

bdd = 1

methods = c("R","DR","S")
method.names = c(methods, paste(methods, "-o", sep = ""))
MSE = matrix(nrow = itern, ncol = length(method.names))
MSE_bdd = matrix(nrow = itern, ncol = length(method.names))
MAPE = matrix(nrow = itern, ncol = length(method.names))
colnames(MSE) = method.names
colnames(MAPE) = method.names

missing_index = NULL
tau_hat.list = vector(mode = "list", length = itern)
CATE.truth.list = vector(mode = "list", length = itern)
for(i in 1:itern){
  missing_i = tryCatch(load(paste("result.", i, ".rda", sep = "")),
                       error = function(cond){NA})
  if(is.na(missing_i)){
    MSE[i,] = NA
    MAPE[i,] = NA
    missing_index = cbind(missing_index, i)
  }else{
    load(paste("result.", i, ".rda", sep = ""))
    CATE.truth.list[[i]] = CATE.truth
    id.bdd = (dat$Z1>-bdd & dat$Z1<bdd & dat$Z2>-bdd & dat$Z2<bdd)

    tau_hat.mx = cbind(R = trunclearner_result$est_R,
                       DR = trunclearner_result$est_DR,
                       S = trunclearner_result$est_S,
                       R_oracle = trunclearner_result_oracle$est_R,
                       DR_oracle = trunclearner_result_oracle$est_DR,
                       S_oracle = trunclearner_result_oracle$est_S)

    tau_hat.list[[i]] = tau_hat.mx

    CATE.truth.mx = replicate(ncol(tau_hat.mx), as.vector(CATE.truth))
    MSE[i,] = colMeans((tau_hat.mx - CATE.truth.mx)^2)
    MAPE[i,] = colMeans(abs(tau_hat.mx - CATE.truth.mx)/CATE.truth.mx)
    
    MSE_bdd[i,] = colMeans((tau_hat.mx[id.bdd,] - CATE.truth.mx[id.bdd,])^2)
  }
}

RMSE = sqrt(MSE)

colMeans(MSE, na.rm = T)
colMeans(MAPE, na.rm = T)

apply(RMSE, 2, median, na.rm = T)

apply(MSE_bdd, 2, median, na.rm = T)


save(tau_hat.list, CATE.truth.list,
     MSE, MAPE,
     MSE_bdd,
     missing_index,
     file = paste("cmb_result_", bdd, ".rda", sep = ""))

