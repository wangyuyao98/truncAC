
## functions for computing the CDF at given times ----------------------------------

# Return a vector of CDF(time.eval_i|Z_i)
CDF_eval <- function(time.eval, CDF.mx){
    if(length(time.eval) != nrow(CDF.mx)){
        stop("The number of time points does not equal the number of subjects!")
    }
    
    jumps = as.numeric(colnames(CDF.mx))
    CDF.mx = cbind(0, CDF.mx)
    
    probs = rep(NA, length(time.eval))
    for(i in 1:length(time.eval)){
        id = findInterval(time.eval[i], c(0, jumps, Inf))
        probs[i] = CDF.mx[i,id]
    }
    
    return(probs)
}



# Return a matrix of CDF(time.eval_j|Z_i) - (i,j)-th element
CDF_eval.mx <- function(time.eval, CDF.mx){
    jumps = as.numeric(colnames(CDF.mx))
    CDF.mx = cbind(0, CDF.mx)
    id = findInterval(time.eval, c(0, jumps, Inf))
    probs = CDF.mx[,id]
    
    if(length(time.eval)>1 & is.null(dim(probs))){
        names(probs) = time.eval
    }else if(length(time.eval)>1 & (!is.null(dim(probs)))){
        colnames(probs) = time.eval
    }
    
    return(probs)
}




# functions for computing L-S integrals ---------------------------------------------

int_fmx_dF <- function(v, f.mx, F.mx){
    
    if(mean(dim(f.mx) == dim(F.mx))<1){
        stop("The dimensions of f.mx and F.mx are not the same!")
    }
    
    jumps = as.numeric(colnames(F.mx))
    dF.mx = F.mx - cbind(0, F.mx[,-ncol(F.mx)]) 
    # dF.mx = F.mx - cbind(F.mx[,1], F.mx[,-ncol(F.mx)]) 
    id = findInterval(v, c(0, jumps, Inf))
    vals = cbind(0, f.mx*dF.mx)
    
    inte = int_sum(id, vals)
    
    rownames(inte) = rownames(F.mx)
    colnames(inte) = v
    
    return(inte)
}


# # compute the integral for one subject at a single time point, i.e., f.mx and F.mx are just vectors, and v is a given time
# int_fmx_dF.single <- function(v, f.vec, F.vec, jumps){
#     if(length(f.vec) != length(F.vec) | length(f.vec) != length(jumps)){
#         stop("The lengths of f.mx, F.mx, and jumps are not all the same!")
#     }
#     
#     dF.vec = F.vec - c(0, F.vec[-length(F.vec)])
#     id = findInterval(v, c(0, jumps, Inf))
#     vals = c(0, f.vec*dF.vec)
#     
#     inte = sum(vals[1:id])
#     
#     return(inte)
# }



# int_fmx_dS <- function(v, f.mx, F.mx){
#     
#     if(mean(dim(f.mx) == dim(F.mx))<1){
#         stop("The dimensions of f.mx and F.mx are not the same!")
#     }
#     
#     jumps = as.numeric(colnames(F.mx))
#     dF.mx = F.mx - cbind(1, F.mx[,-ncol(F.mx)]) 
#     # dF.mx = F.mx - cbind(F.mx[,1], F.mx[,-ncol(F.mx)]) 
#     id = findInterval(v, c(0, jumps, Inf))
#     vals = cbind(0, f.mx*dF.mx)
#     
#     inte = int_sum(id, vals)
#     
#     rownames(inte) = rownames(F.mx)
#     colnames(inte) = v
#     
#     return(inte)
# }






# function for returning a matric for \int_v^\infty fmx dF
int_infty_fmx_dF <- function(v, f.mx, F.mx){   
    if(mean(dim(f.mx) == dim(F.mx))<1){
        stop("The dimensions of f.mx and F.mx are not the same!")
    }
    jumps = as.numeric(colnames(F.mx))
    if(nrow(f.mx)==1){
        dF.mx = F.mx - cbind(0, matrix(F.mx[,-ncol(F.mx)], nrow=1) ) 
    }else{
        dF.mx = F.mx - cbind(0, F.mx[,-ncol(F.mx)])
    }
    
    id = findInterval(v, c(0, jumps, Inf))
    vals = cbind(0, f.mx*dF.mx)
    
    inte = int_sum_infty(id, vals)
    
    # rownames(inte) = rownames(F.mx)
    # colnames(inte) = v
    
    return(inte)
}



# id: a single index denoting the 
int_sum.single <- function(id, vals){
    if(id==1){
        inte = vals[,1]
    }else{
        inte = rowSums(vals[,(1:id)])
    }
    return(inte)
}

int_sum <- function(id, vals){
    inte = sapply(id, int_sum.single, vals = vals)
    return(inte)
}


# id: a single index denoting the 
int_sum_infty.single <- function(id, vals){
    if(id == ncol(vals)){
        inte = rep(0, nrow(vals))
    }else if(id == ncol(vals)-1){
        inte = vals[,ncol(vals)]
    }else{
        if(nrow(vals) == 1){
            inte = sum(vals[,((id+1):ncol(vals))])
        }else{
            inte = rowSums(vals[,((id+1):ncol(vals))])
        }
    }
    return(inte)
}

int_sum_infty <- function(id, vals){
    inte = sapply(id, int_sum_infty.single, vals = vals)
    return(inte)
}



