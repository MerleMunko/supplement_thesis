
# load packages
library(parallel)
library(GFDmcv)

# source functions
source("functions.R")

# Set the parameters for the simulation
Nsim  <- 2000
Nres  <- 1000

# Set the options for the parameters
# vector of group sizes
n_vec_list       <- list(balanced_large    = 25*c(60,60,60,60),
                         unbalanced_large  = 25*c(128,44,52,16),
                         balanced_medium   = 5*c(60,60,60,60),
                         unbalanced_medium = 5*c(128,44,52,16),
                         balanced_small    = c(60,60,60,60),
                         unbalanced_small  = c(128,44,52,16)
)

# The survival distribution setting
distr_vec       <- c("pwExp", "exp_early", "exp_late", "exp_prop", "Weib_scale",
                     "Weib_prop", "Weib_shape", "Weib_late", "logn")[setting[1]]

# The censoring distribution setting
cens_dis_vec    <- c("weib_eq", "weib_uneq+low", "weib_uneq+high")

# discreet or continuous case?
discreet_vec <- c("continuous"=F,"discreet"=T)

# difference of the RMSTs for alt_ind \neq 0 
# for delta = 0, we choose the data of group(s) alt_ind from the alternative 
# distr. with the same RMST.
if(setting[1] %in% c(3,4,8)){
  delta_vec       <- 1.5
}else{ delta_vec       <- c(0,1.5)   }


Settings <- expand.grid(distr_vec, cens_dis_vec, 1:length(n_vec_list), 
                        delta_vec, discreet_vec)
n_names <- names(n_vec_list)
d_names <- names(discreet_vec)
rownames(Settings) <- set_names <- apply(Settings, 1, function(x){
  paste(x[4], n_names[as.numeric(x[3])], x[1], x[2], 
        d_names[as.logical(x[5])+1], sep = "_")
})


# define the methods which should be considered
methods <- c("asymptotic", "asymptotic_bonf", 
             #"pooledBS", "pooledBS_bonf", "wild_poisson", "wild_gaussian", 
             #"wild_rademacher", "wild_mammen","groupwiseBS", 
             "perm_bonf") # "rand_bonf"
m <- length(methods)

# which method for calculating the critical value?
crit.value.method <- c("equi", "inequi")
# group(s) following the alternative, 0 = all groups under the null hypothesis
alt_ind     <- 4
# number of groups
k <- 4
# number of risks
M <- 3
kM <- k*M
# probabilities for the different risks
D_probs <- rbind(c(0.33,0.25,0.42),
                 c(0.33,0.25,0.42),
                 c(0.33,0.25,0.42),
                 c(0.33,0.25,0.42))
# construct hypothesis matrices and vectors
I <- diag(M)
c_mat_list <- list("2by2" = rows2list(rbind(cbind(I,-I,I,-I),cbind(I,I,-I,-I),
                                            cbind(I,-I,-I,I))),
                   "Dunnett" = rows2list(kronecker(contr_mat(4, "Dunnett"),
                                                   diag(M))),
                   "Tukey" = rows2list(kronecker(contr_mat(4, "Tukey"),diag(M))))
c_vec_list <- list("2by2" = lapply(1:(3*M),function(x) 0),
                   "Dunnett" = lapply(1:(M*3), function(x) 0),
                   "Tukey" = lapply(1:(M*6), function(x) 0))
C <- length(c_mat_list)
namen <- names(c_mat_list)



for(set.ind in 1:nrow(Settings)){
  print(paste("Setting", set.ind, "of", nrow(Settings)))
  
  # Set the parameters for the data generation
  n_vec       <- n_vec_list[[Settings[set.ind,3]]]   # vector of group sizes
  distr       <- Settings[set.ind,1]                 # The survival distribution setting
  cens_dis    <- Settings[set.ind,2]                 # The censoring distribution setting
  delta       <- Settings[set.ind,4]                 # difference of the RMSTs
  discreet    <- Settings[set.ind,5]                 # discreet or continuous distributions?
  
  # generate and save data
  source("data.R")
  set.seed(1)
  #  if (file.exists(paste(path,"/data/data_",set_names[set.ind],".Rdata", sep=""))) {
  #    load(paste(path,"/data/data_",set_names[set.ind],".Rdata", sep=""))
  #  }else{
  all_data <- replicate(Nsim, data_gen(n_vec, sum(n_vec)), simplify = FALSE)
  #save(all_data, file = paste(path,"/data/data_",set_names[set.ind],".Rdata", 
  #                            sep=""))
  #  }
  
  #if (!file.exists(paste(path, "/results/errorrate/",
  #                       set_names[set.ind],".Rdata", sep=""))) {
    source("Simulation.R")
    # results for the global error rates
    save(results, file = paste(path, "/results/errorrate/",set_names[set.ind],
                               ".Rdata", sep=""))
    # results for the local error rates
    save(results_multiple, file = paste(path, 
                                        "/results/multiple_errorrate/multiple_",
                                        set_names[set.ind],".Rdata", sep=""))
  #}
  
}


