

# function for one simulation run with data huhu
simu <- function(huhu){
  if(huhu %% 500 == 0) print(paste("Simulation run", huhu))
  my_data <- all_data[[huhu]]
  
  test_stat_list <- wrap_multiple_test_stat(my_data, tau, c_mat_list, c_vec_list, M, k)
  times <- sort(unique(my_data$X[my_data$D > 0]))
  # calculate sqrt_Sigma_hat and est_list
  sqrt_Sigma_hat <- matrix(0,ncol=kM,nrow=kM)
  est_list <- list()
  for(i in 1:k) {
    values2 <- my_data[my_data$group == i,]
    
    #calculate for each group
    temp <- RMTL(values2$X, values2$D, tau, M = M)
    est_list[[i]] <- estimators(values2$X, values2$D, tau, M, times = times)
    
    # eigen value decomposition of var
    eigen_temp <- eigen(temp[["var_rmtl"]])
    Mitte <- matrix(0,nrow=M,ncol=M)
    diag(Mitte)[eigen_temp$values > 0] <- ((eigen_temp$values)[eigen_temp$values > 0])^(1/2)
    # determine the inverse root of var
    sqrt_Sigma_hat[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- eigen_temp$vectors%*%Mitte%*%t(eigen_temp$vectors)
  }
  ginvs <- lapply(c_mat_list, function(x) lapply(x, function(c_mat){
    MASS::ginv(c_mat %*% test_stat_list$var %*% t(c_mat))
  }))
  
  # empty decision matrices
  decision_vec <- list()
  
  if("asymptotic" %in% methods)  decision_vec["asymptotic"] <- list(wrap_asymptotics(test_stat_list, c_mat_list, kM, crit.value.method, Nres, ginvs))
  if("asymptotic_bonf" %in% methods)  decision_vec["asymptotic_bonf"] <- list(wrap_asymptotics_bonf(test_stat_list$teststats, c_mat_list, kM, crit.value.method))
  
  if("pooledBS" %in% methods | "pooledBS_bonf" %in% methods)  decision_vec[c("pooledBS", "pooledBS_bonf")] <- wrap_pooledBS(my_data, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, sqrt_Sigma_hat, ginvs)
  
  if("wild_gaussian" %in% methods)  decision_vec["wild_gaussian"] <- list(wrap_wildBS(my_data, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, est_list, bs_function = rnorm)) # gaussian
  if("wild_rademacher" %in% methods)  decision_vec["wild_rademacher"] <- list(wrap_wildBS(my_data, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, est_list, bs_function = function(n) 2*rbinom(n,1,0.5)-1)) # rademacher
  if("wild_mammen" %in% methods) decision_vec["wild_mammen"] <- list(wrap_wildBS(my_data, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, est_list, bs_function = function(n) sample(mammen, n, replace=TRUE, prob = probs_mammen))) # mammen
  if("wild_poisson" %in% methods)  decision_vec["wild_poisson"] <- list(wrap_wildBS(my_data, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, est_list, bs_function = function(n) rpois(n,1) - 1)) # poisson
  
  if("groupwiseBS" %in% methods)  decision_vec["groupwiseBS"] <- list(wrap_groupwiseBS(my_data, tau, c_mat_list, Nres, test_stat_list, n_vec, crit.value.method, M, kM))
  
  if("perm_bonf" %in% methods)  decision_vec["perm_bonf"] <- list(wrap_perm_bonf(my_data, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM))
  if("perm" %in% methods)  decision_vec["perm"] <- list(wrap_perm(my_data, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM))
  
  if("rand_bonf" %in% methods)  decision_vec["rand_bonf"] <- list(wrap_rand(my_data, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, ginvs))
  
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
      for(method in methods){
        decision_out[[c_ind]][method,,] <- liste[[method]]
      }}
  }else{
    for(method in methods){
      if(is.list(decision_vec[[method]])){
        decision_out[[1]][method,,] <- decision_vec[[method]][[1]]
      }else{
        decision_out[[1]][method,,] <- decision_vec[[method]]
      } 
    }
  }
  
  
  return(decision_out)
}

# start the simulation
RNGkind("L'Ecuyer-CMRG") # Needed that the parallel kernels do not start with the same seed (set also mc.set.seed = TRUE)
set.seed(1)

# measure the time
start <- Sys.time()

# simulate...
decision_list <- mclapply(1:Nsim, simu, mc.cores = cores)

# stop the time
(time <- difftime(start,Sys.time()))

# put the output into a nice array
decision_arr <- simplify2array(decision_list)
results_multiple <- lapply(1:C, function(c_ind) rowMeans(simplify2array(decision_arr[c_ind,]),dims=3))
names(results_multiple) <- namen
results          <- simplify2array(lapply(1:C, function(c_ind){
  rowMeans(apply(simplify2array(decision_arr[c_ind,]),c(1,3,4),function(x) sum(x) > 0), dims=2)
}))

#save(decision_arr, file = paste(path, "/results/decision_list/dec_",set_names[set.ind],".Rdata", sep=""))
