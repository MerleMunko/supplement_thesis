### Simulation for RMST ANOVA and Multiple Tests ###
library(multcomp)
library(parallel)
#library(survival)

# Set the parameters for the simulation
Nsim  <- 5000
#Nres <- 1999
Nres  <- 19999

# Set the options for the parameters
# vector of group sizes
n_vec_list       <- list(data = c(450,481,654,649))

# group(s) following the alternative, 0 = all groups under the null hypothesis
alt_ind     <- c(0,1)          


Settings <- expand.grid(1:length(n_vec_list), alt_ind)
n_names <- names(n_vec_list)
rownames(Settings) <- set_names <- apply(Settings, 1, function(x){ 
  as.character(x[2])
})

# source the functions
source("functions_multiple.R")
source("resampling_multiple.R")

# define the methods which should be considered
#methods <- c("asymptotic_global", "permutation", "pooled", "wild, Rademacher", "wild, Gaussian", "groupwise", "parametric")
methods <- c("asymptotic_global", "permutation", "asymptotic", "wild, Rademacher", 
             "wild, Gaussian", "groupwise", "asymptotic_bonf", "permutation_bonf")
m <- length(methods)

# which method for calculating the critical value?
# equicoordinate means calculating the quantiles of the maximum statistic; only possible for multiple contrast matrices with the same rank
# inequicoordinate means calculating the quantiles for each coordinate seperately; possible for all types of multiple contrast matrices
crit.value.method <- c("equi", "inequi")


for(set.ind in 1:nrow(Settings)){
  print(paste("Setting", set.ind, "of", nrow(Settings)))
  
  # Set the parameters for the data generation
  n_vec       <- n_vec_list[[Settings[set.ind,1]]]   # vector of group sizes
  #distr       <- Settings[set.ind,1]                 # The survival distribution setting
  #cens_dis    <- Settings[set.ind,2]                 # The censoring distribution setting
  alt_ind       <- Settings[set.ind,2]                
  
  
  source("Simulation_multiple_2x2.R")
  
  # censoring rates
  save(censp_s, file = paste(path,"/results_data/censrate/cens_",set_names[set.ind],"_unbalanced_large_Kaplan_Meier_Kaplan_Meier.Rdata", sep=""))
  # results for the global error rates
  save(results, file = paste(path, "/results_data/errorrate/",set_names[set.ind],"_unbalanced_large_Kaplan_Meier_Kaplan_Meier.Rdata", sep=""))
  # results for the local error rates
  save(results_multiple, file = paste(path, "/results_data/multiple_errorrate/multiple_",set_names[set.ind],"_unbalanced_large_Kaplan_Meier_Kaplan_Meier.Rdata", sep=""))
  # decision list for all decisions
  save(decision_arr, file = paste(path, "/results_data/decision_list/decision_",set_names[set.ind],"_unbalanced_large_Kaplan_Meier_Kaplan_Meier.Rdata", sep=""))
  
}


