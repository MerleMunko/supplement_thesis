
####### Pooled Bootstrap #######

# function for calculating the test statistic for the pooled Bootstrap
# Sigma_hat: diagonal matrix of the variance estimators of the original data of the RMSTs
# mu_hat: pooled RMST estimator
test_stat_pooledBS <- function(bs_values, tau, c_mat, Sigma_hat, mu_hat, n_group){
  group <- bs_values$group
  rmst <- numeric(n_group)
  var <- numeric(n_group)
  
  for(i in 1:n_group) {
    values2 <- bs_values[group == i,]
    #calculate for each group
    temp <- RMST(values2$time, values2$status, tau)
    rmst[i] <- temp["rmst"]
    var[i] <- temp["var_rmst"]
  }
  
  #calculate some vectors and matrices needed more than once
  inv_sqrt_var_bs <- numeric(n_group)
  inv_sqrt_var_bs[var != 0] <- 1/sqrt(var[var != 0]) #if the BS variance is 0, we use the g-inverse
  inv_sqrt_Sigma_hat_bs <- diag(inv_sqrt_var_bs)
  vec <- c_mat %*% sqrt(Sigma_hat) %*% inv_sqrt_Sigma_hat_bs %*% (rmst- mu_hat)
  
  #calculate the pooled BS test statistic
  return( t(vec) %*% MASS::ginv(c_mat %*% Sigma_hat %*% t(c_mat)) %*% vec )  
  # the factor n is eliminated by the other factors
  # appearing in the formula for the test statistic and
  # the factor n_i in the variance formula
}

# wrapper for the pooled bootstrap
wrap_pooledBS <- function(values, tau, c_mat_list, Nres, test_stat_list, crit.value.method){
  
  # calculate Sigma_hat and the test statistic
  group <- values$group
  n_group <- max(group)
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)
  
  Sigma_hat <- diag(test_stat_list$var)
  teststats <- test_stat_list$teststats
  
  #calculate mu_hat for the pooled sample
  temp <- RMST(values$time, values$status, tau)
  mu_hat <- temp["rmst"]
  
  erg <- replicate(Nres , expr = {
    # generate a pooled Bootstrap sample
    bs_values <- values[sample(1:n,n,replace = TRUE),c("time", "status")]
    bs_values$group <- group
    bs_values <- sort_data(bs_values)
    
    # calculate the pooled BS test statistics for all matrices in c_mat_list
    erg_int <- list()
    for( hyp in 1:C){
      erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
        test_stat_pooledBS(bs_values=bs_values, tau=tau, c_mat=c_mat, 
                           Sigma_hat=Sigma_hat, mu_hat=mu_hat, n_group = n_group)
      })
    }
    names(erg_int) <- namen
    erg_int
  }, simplify = TRUE)
  dim(erg) <- c(C,Nres)
  
  # save the test statistics in arrays in a list
  ts_pooledBS <- list()
  for( hyp in 1:C){
    dim <- length(erg[hyp,][[1]])
    ts_pooledBS[[hyp]] <- matrix(unlist(erg[hyp,]), ncol = dim, byrow=TRUE)
  }
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out <- matrix(numeric(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- reject_multiple(teststats[[ind]], crit_values2(t(ts_pooledBS[[ind]])))
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- reject_multiple(teststats[[ind]], quantile(apply(ts_pooledBS[[ind]], 1, max), probs = 0.95))
    
    out
  })
  
  return(decision)
}


####### Wild Bootstrap ########

# function for calculating the test statistic for the wild Bootstrap
# bs_weights = G_ij are the bootstrap multipliers
test_stat_wildBS <- function(survtab_list, bs_weights, tau, c_mat, n_group, weight_ind){

  wild_rmst <- numeric(n_group)
  wild_var <- numeric(n_group)
  
  for(i in 1:n_group) {
    
    area <- survtab_list[[i]]$area
    d <- survtab_list[[i]]$d
    Y <- survtab_list[[i]]$Y
    
    weight <- bs_weights[weight_ind[[i]]]
    
    G <- G2 <- numeric(length(d))
    l <- 1
    for(j in 1:length(d)){
      if(d[j] > 0){
      G[j] <- sum(weight[l:(l+d[j]-1)])
      G2[j] <- sum((weight[l:(l+d[j]-1)])^2)
      
      l <- l+d[j]}
    }
    
    wild_rmst[i] <- sum(G * ( 1/sqrt(Y)/sqrt(Y-d)) * area, na.rm = TRUE) # \hat{mu}_i^G
    # ACHTUNG hier ist irgendwo ein Fehler
    
    wild_var[i] <- sum( G2 * ( 1/Y/(Y-d)) * area^2  , na.rm = TRUE)
    
  }
  
  #calculate some vectors and matrices needed more than once
  vec <- c_mat %*% wild_rmst
  
  #calculate the pooled BS test statistic
  return( t(vec) %*% MASS::ginv(c_mat %*% diag(wild_var) %*% t(c_mat)) %*% vec )  
  # the factor n is eliminated by the other factors
  # appearing in the formula for the test statistic and
  # the factor n_i in the variance formula
}

# wrapper for the wild bootstrap
wrap_wildBS <- function(values, tau, c_mat_list, Nres, test_stat_list, n_vec,
                          multiplier = c("Rademacher", "Gaussian"), crit.value.method){
  
  # calculate the test statistics
  group <- values$group
  n_group <- max(group)
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)
  
  teststats <- test_stat_list$teststats
  
  # generate wild bootstrap multipliers
  multiplier <- match.arg(multiplier)
  if(multiplier == "Rademacher") G_mat <- matrix(2*(rbinom(n*Nres, size = 1, prob = 0.5) - 0.5),ncol = n)
  if(multiplier == "Gaussian") G_mat <- matrix(rnorm(n*Nres),ncol = n)
  
  # save d, Y and area for all groups i
  survtab_list <- weight_ind <- list()
  for(i in 1:n_group){
    values2 <- values[values$group == i,]
    survtab <- simple_surv(values2$time, values2$status)
    
    t_i <- survtab[,1] <= tau # identify time points t_i <= tau
    t_sel <- survtab[,1][t_i] # select relevent time points
    S_km <- survtab[,2][t_i] # calculate Kaplan-Meier estimator
    
    
    #w.factor <- diff(c(0, t_sel,tau)) # width of the area under S(t) from t_0=0
    
    survtab_list[[i]] <- list(area = rev(cumsum(rev(diff(c(t_sel,tau)) * S_km))), # calculate areas under S(t) from t_i to tau for
                              d = survtab[,3][t_i], # determine number of events
                              Y = survtab[,4][t_i]) # determine number of individuals under risk
    
    if(i>1){
      weight_ind[[i]] <- sum(n_vec[1:(i-1)])+(1:(n_vec[i]))
    }else{ weight_ind[[i]] <- 1:n_vec[i]}
  }
  
  erg <- sapply(1:Nres , function(ind){
    
    
    # calculate the wild BS test statistics for all matrices in c_mat_list
    erg_int <- list()
    
    for( hyp in 1:C){
      erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
        test_stat_wildBS(survtab_list, bs_weights = G_mat[ind,], tau=tau, 
                         c_mat=c_mat, n_group = n_group, weight_ind = weight_ind)
      })
    }
    names(erg_int) <- namen
    erg_int
  }, simplify = TRUE)
  dim(erg) <- c(C,Nres)
  
  # save the test statistics in arrays in a list
  ts_wildBS <- list()
  for( hyp in 1:C){
    dim <- length(erg[hyp,][[1]])
    ts_wildBS[[hyp]] <- matrix(unlist(erg[hyp,]), ncol = dim, byrow=TRUE)
  }
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out <- matrix(numeric(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- reject_multiple(teststats[[ind]], crit_values2(t(ts_wildBS[[ind]])))
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- reject_multiple(teststats[[ind]], quantile(apply(ts_wildBS[[ind]], 1, max), probs = 0.95))
    
    out
  })
  
  return(decision)
}


####### groupwise Bootstrap #######

#the test statistic is calculated in the wrapper

# wrapper for the groupwise bootstrap
wrap_groupwiseBS <- function(values, tau, c_mat_list, Nres, test_stat_list, n_vec, crit.value.method){
  
  # calculate the test statistics
  group <- values$group
  n_group <- max(group)
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)
  
  teststats <- test_stat_list$teststats
  rmst <- test_stat_list$rmst
  
  erg <- replicate(Nres , expr = {
    
    # generate the groupwiese bootstrap sample
    #bs_values <- data.frame()
    #for(i in 1:n_group){
    #  bs_values <- rbind(bs_values, values[group == i,])
    #}
    #bs_values <- sort_data(bs_values)
    
    # calculate the rmst and variance estimator for the bootstrap sample
    bs_rmst <- numeric(n_group)
    bs_var <- numeric(n_group)
        
    for(i in 1:n_group) {
      # generate a bootstrap sample in group i
          values2 <- values[group == i,][sample(1:n_vec[i],n_vec[i],replace = TRUE),]
          values2 <- sort_data(values2)
          #calculate for each group
          temp <- RMST(values2$time, values2$status, tau)
          bs_rmst[i] <- temp["rmst"]
          bs_var[i] <- temp["var_rmst"]
    }
    
    # calculate the groupwise BS test statistics for all matrices in c_mat_list
    erg_int <- list()
    for( hyp in 1:C){
      erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
        vec <- c_mat %*% (bs_rmst- rmst)
        t(vec) %*% MASS::ginv(c_mat %*% diag(bs_var) %*% t(c_mat)) %*% vec  
      })
    }
    names(erg_int) <- namen
    erg_int
  }, simplify = TRUE)
  dim(erg) <- c(C,Nres)
  
  # save the test statistics in arrays in a list
  ts_groupwiseBS <- list()
  for( hyp in 1:C){
    dim <- length(erg[hyp,][[1]])
    ts_groupwiseBS[[hyp]] <- matrix(unlist(erg[hyp,]), ncol = dim, byrow=TRUE)
  }
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out <- matrix(numeric(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- reject_multiple(teststats[[ind]], crit_values2(t(ts_groupwiseBS[[ind]])))
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- reject_multiple(teststats[[ind]], quantile(apply(ts_groupwiseBS[[ind]], 1, max), probs = 0.95))
    
    out
  })
  
  return(decision)
}


####### parametric Bootstrap #######

# function for calculating the test statistic for the parametric Bootstrap
# Sigma_hat: diagonal matrix of the variance estimators of the original data of the RMSTs
test_stat_parametricBS <- function(Z_values, c_mat){

  n <- nrow(Z_values)
  
  #calculate some vectors and matrices needed more than once
  vec <- c_mat %*% colMeans(Z_values)
  
  #calculate the parametric BS test statistic
  return( n * t(vec) %*% MASS::ginv(c_mat %*% diag(apply(Z_values, 2, var)) %*% t(c_mat)) %*% vec )  
}

# wrapper for the parametric bootstrap
wrap_parametricBS <- function(values, c_mat_list, Nres, test_stat_list, crit.value.method){
  
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)
  
  teststats <- test_stat_list$teststats
  Sigma_hat <- 1/n * diag(test_stat_list$var) # here, we use the factor 1/n 
  
  # generate the Z values
  Z_list <- lapply(1:Nres, function(hallo) rmvnorm(n, sigma = Sigma_hat))
  
  erg <- sapply(1:Nres , function(ind){

    erg_int <- list()
    for( hyp in 1:C){
      erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
        test_stat_parametricBS(Z_values = Z_list[[ind]], c_mat=c_mat)
      })
    }
    names(erg_int) <- namen
    erg_int
  }, simplify = TRUE)
  dim(erg) <- c(C,Nres)
  
  # save the test statistics in arrays in a list
  ts_parametricBS <- list()
  for( hyp in 1:C){
    dim <- length(erg[hyp,][[1]])
    ts_parametricBS[[hyp]] <- matrix(unlist(erg[hyp,]), ncol = dim, byrow=TRUE)
  }
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out <- matrix(numeric(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- reject_multiple(teststats[[ind]], crit_values2(t(ts_parametricBS[[ind]])))
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- reject_multiple(teststats[[ind]], quantile(apply(ts_parametricBS[[ind]], 1, max), probs = 0.95))
    
    out
  })
  
  return(decision)
}


####### Permutation test for the global Hypothesis #######

#the test statistic is calculated as in test_stat

# wrapper for the studentized permutation approach
wrap_perm <- function(values, tau, global_c_mat_list, Nres, global_teststats, crit.value.method, alpha = 0.05){
  
  # calculate the test statistics
  group <- values$group
  C <- length(global_c_mat_list)
  
  perm_values <- values[,c("time", "status")]
  
  erg <- replicate(Nres , expr = {
    
    # generate the permuted sample
    perm_values$group <- sample(group)
    
    # calculate the perm. test statistics for all matrices in global_c_mat_list
    erg_int <- sapply(global_c_mat_list, function(hyp) test_stat(perm_values, tau, hyp))

    erg_int
  }, simplify = TRUE)
  dim(erg) <- c(C,Nres)
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out_temp <- reject(global_teststats[ind], crit_values2(t(erg[ind,]), alpha = alpha))
    out <- rep(out_temp,length(crit.value.method))
    
    return(out)
  })
  
  return(t(decision))
}

# wrapper for the permutation test for multiple hypotheses with bonferroni correction
wrap_perm_bonf <- function(values, tau, c_mat_list, Nres, teststats, crit.value.method, alpha = 0.05){
  
  l <- length(crit.value.method)
  out <- lapply(1:C, function(c_ind){
    bonf_ergs <- wrap_perm(values, tau, c_mat_list[[c_ind]], Nres, teststats[[c_ind]], crit.value.method, alpha = alpha/length(c_mat_list[[c_ind]]))
    dim(bonf_ergs) <- c(length(c_mat_list[[c_ind]]), l)
    colnames(bonf_ergs) <- crit.value.method
    rownames(bonf_ergs) <- names(c_mat_list[[c_ind]])
    bonf_ergs
  })
  
  out
}
