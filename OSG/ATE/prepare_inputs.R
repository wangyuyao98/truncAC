# code for preparing the inputs for OSG
rm(list = ls())
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# ### code for writing indices into a txt file
# itern = 500
# write(1:itern, file = paste("indices_1-", itern, ".txt", sep = ""), sep = "\n")

# itern = 2
# write(1:itern, file = paste("indices_1-", itern, ".txt", sep = ""), sep = "\n")


### create folders that contains the seeds

# load the seeds for simulating datasets
load("seeds.rda")

# simulate the seeds for boostrap
itern = 500
n.boot = 100    # the total number of bootstrap runs
nb.each = 20      # the number of bootstrap in each job
ns = ceiling(n.boot/nb.each)

set.seed(213)
seeds.boot = matrix(sample(100000, size = ns*itern), nrow = itern, ncol = ns)

dir.create("input_seeds")
for(i in 1:itern){
    for(j in 1:ns){
        index = ns*(i-1) + j
        folder = paste("input_seeds/input_", index, sep = "")
        dir.create(folder)
        seed = seeds[i]
        seed.b = seeds.boot[i,j]
        save(seed, seed.b, file = paste(folder, "/seeds_input.rda", sep = ""))
    }
}

write(1:index, file = paste("indices", itern, ".txt", sep = ""), sep = "\n")





# ### Previous code that need to submit the jobs multiple times
# ### code for writing indices into a txt file
# itern = 500
# write(1:itern, file = paste("indices_1-", itern, ".txt", sep = ""), sep = "\n")
# 
# # itern = 2
# # write(1:itern, file = paste("indices_1-", itern, ".txt", sep = ""), sep = "\n")
# 
# 
# ### create folders that contains the seeds
# 
# # load the seeds for simulating datasets
# load("seeds.rda")
# 
# # simulate the seeds for boostrap
# n.boot = 100
# nb.each = 100
# ns = ceiling(n.boot/nb.each)
# set.seed(123)
# seeds.boot = matrix(sample(100000, size = ns*itern), nrow = itern, ncol = ns)
# 
# for(j in 1:ns){
#     dir.create(paste("input_b", j, sep = ""))
#     for(i in 1:itern){
#         folder = paste("input_b", j, "/data", i, sep = "")
#         dir.create(folder)
#         seed = seeds[i]
#         seed.b = seeds.boot[i,j]
#         save(seed, seed.b, file = paste(folder, "/seeds_input.rda", sep = ""))
#     }
# }






