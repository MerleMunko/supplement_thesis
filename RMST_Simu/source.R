### Simulation for RMST ANOVA and Multiple Tests ###
library(multcomp)
library(parallel)

# Set the parameters for the simulation
Nsim  <- 5000
Nres  <- 1999

# Set the options for the parameters
# vector of group sizes
n_vec_list       <- list(balanced_large    = 4*c(15,15,15,15),
                         unbalanced_large  = 4*c(10,20,10,20),
                         balanced_medium   = 2*c(15,15,15,15),
                         unbalanced_medium = 2*c(10,20,10,20),
                         balanced_small    = c(15,15,15,15),
                         unbalanced_small  = c(10,20,10,20))
# The survival distribution setting
#distr_vec       <- c("pwExp", "exp_early", "exp_late", "exp_prop", "Weib_pwExp", 
#                     "Weib_scale", "Weib_prop", "Weib_shape", "Weib_late", "Weib_early", "logn")  
distr_vec       <- c("pwExp", "exp_early", "exp_late", "exp_prop", "Weib_pwExp", 
                     "Weib_scale", "Weib_prop", "Weib_shape", "Weib_late", "Weib_early", "logn")[setting[1]]  

# The censoring distribution setting
#cens_dis_vec    <- c("weib_eq", "weib_uneq", "unif_eq", "unif_uneq", "exp_eq", "exp_uneq")
#cens_dis_vec    <- c("weib_eq", "weib_uneq", "unif_eq", "unif_uneq", "exp_eq", "exp_uneq")[setting[2]]
cens_dis_vec    <- c("weib_eq", "weib_uneq+low", "weib_uneq+high", "weib_eq+high")[setting[2]]

# group(s) following the alternative, 0 = all groups under the null hypothesis
alt_ind     <- 4             

# difference of the RMSTs for alt_ind \neq 0 # for delta = 0, we choose the data of group(s) alt_ind from the alternative distr. with the same RMST.
#delta_vec       <- c(0,0.5,1,1.5)                
delta_vec       <- c(0,1.5)                


Settings <- expand.grid(distr_vec, cens_dis_vec, 1:length(n_vec_list), delta_vec)
n_names <- names(n_vec_list)
rownames(Settings) <- set_names <- apply(Settings, 1, function(x){ 
  paste(x[4], n_names[as.numeric(x[3])], x[1], x[2], sep = "_")
})

# source the functions
source("functions_multiple.R")
source("resampling_multiple.R")

# define the methods which should be considered
#methods <- c("asymptotic_global", "permutation", "pooled", "wild, Rademacher", "wild, Gaussian", "groupwise", "parametric")
methods <- c("asymptotic_global", "permutation", "asymptotic", "wild, Rademacher", "wild, Gaussian", "groupwise", "asymptotic_bonf", "permutation_bonf")
m <- length(methods)

# which method for calculating the critical value?
# equicoordinate means calculating the quantiles of the maximum statistic; only possible for multiple contrast matrices with the same rank
# inequicoordinate means calculating the quantiles for each coordinate seperately; possible for all types of multiple contrast matrices
crit.value.method <- c("equi", "inequi")

for(set.ind in 1:nrow(Settings)){
  print(paste("Setting", set.ind, "of", nrow(Settings)))
  
  # Set the parameters for the data generation
  n_vec       <- n_vec_list[[Settings[set.ind,3]]]   # vector of group sizes
  distr       <- Settings[set.ind,1]                 # The survival distribution setting
  cens_dis    <- Settings[set.ind,2]                 # The censoring distribution setting
  delta       <- Settings[set.ind,4]                 # difference of the RMSTs
  #alt_ind     <- as.numeric(delta != 0)              # group(s) following the alternative, 0 = null hypothesis
  
  source("Simulation_multiple.R")
  
  # censoring rates
  save(censp_s, file = paste(path,"/results/censrate/cens_",set_names[set.ind],".Rdata", sep=""))
  # results for the global error rates
  save(results, file = paste(path, "/results/errorrate/",set_names[set.ind],".Rdata", sep=""))
  # results for the local error rates
  save(results_multiple, file = paste(path, "/results/multiple_errorrate/multiple_",set_names[set.ind],".Rdata", sep=""))
  # decision list for all decisions
  save(decision_arr, file = paste(path, "/results/decision_list/decision_",set_names[set.ind],".Rdata", sep=""))
  
}


