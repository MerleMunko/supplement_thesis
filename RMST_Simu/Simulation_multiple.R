
### The simulation for one setting

# define the contrast matrices as partitionized matrices in a list
k          <- length(n_vec)
c_mat_list <- list(Dunnett   = rows2list(contrMat(rep(1,k), type="Dunnett")),
                   Tukey     = rows2list(contrMat(rep(1,k), type="Tukey")),
                   GrandMean = rows2list(contrMat(rep(1,k), type="GrandMean")))
C          <- length(c_mat_list)
namen      <- names(c_mat_list)


# check whether the single contrast matrices have the same rank
if(any(sapply(c_mat_list, function(c_mat) length(unique(sapply(c_mat, function(x) length(x))))!=1)) & "equi" %in% crit.value.method)  stop("Equicoordinate quantiles are not possible for different ranks.")

# check whether the single contrast matrices are vectors
if(any(sapply(c_mat_list, function(c_mat) any(sapply(c_mat, function(x) nrow(x))!=1))) & "asymptotic" %in% methods) stop("Method \"asymptotic\" only works for contrast vectors.")


# source the data generating function
source("data.R")

# estimate the censoring rate
set.seed(1)
d <- data_gen( rep(10^6, k), k*10^6)
censp_s <- numeric(k)
for(i in 1:k){
  censp_s[i] <- mean(1 - d$status[d$group == i])
}

# start the simulation
RNGkind("L'Ecuyer-CMRG") # Needed that the parallel kernels do not start with the same seed (set also mc.set.seed = TRUE)
set.seed(1)

# measure the time
start <- Sys.time()

decision_arr <- simplify2array(mclapply(1:Nsim, function(huhu){
  # information, where we are in our simulation
  if(huhu %% 500 == 0) print(paste("Iteration", huhu, "of", Nsim))
  
  # generate and sort data
  my_data <- data_gen(n_vec, n)
  
  # compute the test statistics and the variances of the rmst
  test_stat_list <- wrap_multiple_test_stat(my_data, c_mat_list, tau)
  
  # compute the global matrices and global teststatistics
  if("asymptotic_global" %in% methods | "permutation" %in% methods){
    global_c_mat_list <- lapply(c_mat_list, function(X) global_mat(X, k))
    global_teststats  <- sapply(global_c_mat_list, function(H) test_stat(my_data, tau, H))
  }
  
  # empty decision matrices
  decision_vec <- list()
  #names(decision_vec) <- list(methods, namen, crit.value.method)
  
  # Global Hypothesis: Asymptotics
  if("asymptotic_global" %in% methods)  decision_vec["asymptotic_global"] <- list(wrap_asymptotics_global(global_teststats, global_c_mat_list, crit.value.method))
  
  # Global Hypothesis: Permutation
  if("permutation" %in% methods)  decision_vec["permutation"] <- list(wrap_perm(my_data, tau, global_c_mat_list, Nres, global_teststats, crit.value.method))
  
  ## Multiple Testing
  # Asymptotics
  if("asymptotic" %in% methods)  decision_vec["asymptotic"] <- list(wrap_asymptotics(test_stat_list, c_mat_list, k, crit.value.method))
  
  # pooled Bootstrap test
  if("pooled" %in% methods)  decision_vec["pooled"] <- list(wrap_pooledBS(my_data, tau, c_mat_list, Nres, test_stat_list, crit.value.method))
  
  # wild Bootstrap test, Rademacher multipliers
  if("wild, Rademacher" %in% methods)  decision_vec["wild, Rademacher"] <- list(wrap_wildBS(my_data, tau, c_mat_list, Nres, test_stat_list, n_vec, multiplier = "Rademacher", crit.value.method))
  
  # wild Bootstrap test, Gaussian multipliers
  if("wild, Gaussian" %in% methods)  decision_vec["wild, Gaussian"] <- list(wrap_wildBS(my_data, tau, c_mat_list, Nres, test_stat_list, n_vec, multiplier = "Gaussian", crit.value.method))
  
  # groupwise Bootstrap test
  if("groupwise" %in% methods)  decision_vec["groupwise"] <- list(wrap_groupwiseBS(my_data, tau, c_mat_list, Nres, test_stat_list, n_vec, crit.value.method))
  
  # parametric Bootstrap test
  if("parametric" %in% methods)  decision_vec["parametric"] <- list(wrap_parametricBS(my_data, c_mat_list, Nres, test_stat_list, crit.value.method))
  
  ## Bonferroni approaches
  # Global Hypothesis: Asymptotics
  if("asymptotic_bonf" %in% methods)  decision_vec["asymptotic_bonf"] <- list(wrap_asymptotics_bonf(test_stat_list$teststats, c_mat_list, crit.value.method))
  
  # Global Hypothesis: Permutation
  if("permutation_bonf" %in% methods)  decision_vec["permutation_bonf"] <- list(wrap_perm_bonf(my_data, tau, c_mat_list, Nres, test_stat_list$teststats, crit.value.method))
  
  
  # empty output list
  decision_out <- lapply(1:C, function(c_ind){
    out_arr           <- array(NA, dim = c(m, length(c_mat_list[[c_ind]]), length(crit.value.method)))
    dimnames(out_arr) <- list(methods, names(c_mat_list[[c_ind]]), crit.value.method) 
    return(out_arr)
  })
  names(decision_out) <- namen
  if(C>1){
    for(c_ind in 1:C){
      liste <- lapply(decision_vec, function(x) x[[c_ind]])
      for(method in 1:m){
        decision_out[[c_ind]][method,,] <- liste[[method]]
      }}
  }else{
    for(method in 1:m){
      if(is.list(decision_vec[[method]])){
        decision_out[[1]][method,,] <- decision_vec[[method]][[1]]
      }else{
        decision_out[[1]][method,,] <- decision_vec[[method]]
      } 
      
    }
  }
  
  # print the time until now
  if(huhu %% 500 == 0) print(difftime(Sys.time(), start))
  
  return(decision_out)
  
}, mc.cores = cores, mc.set.seed = TRUE))

# results for the local and global error rates

results_multiple <- lapply(1:C, function(c_ind) rowMeans(simplify2array(decision_arr[c_ind,]),dims=3))
names(results_multiple) <- namen
results          <- simplify2array(lapply(1:C, function(c_ind){
  rowMeans(apply(simplify2array(decision_arr[c_ind,]),c(1,3,4),function(x) sum(x) > 0), dims=2)
}))
dimnames(results)[[3]] <- namen
results <- aperm(results, c(1,3,2))



