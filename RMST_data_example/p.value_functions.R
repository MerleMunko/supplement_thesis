# p.values
p.value2 <- function(data_mat, alpha = 0.05, teststat){
  
  n <- length(data_mat[1,])
  dimension <- length(data_mat[,1])
  
  # First forget the connection of the coordinates, and sort the values per coordinate
  data_order <- t( apply(data_mat, 1, sort) )
  
  # worst case1
  #for each point outside the box only one coordinate leads to a rejection. 
  # Thus, for each coordinate (alpha/dim) * number of obs are outside the box  (Bonferoni correction)
  j_low <- ceiling( (alpha/dimension) * n ) - 1
  A <- data_mat / data_order[,n - j_low]
  A[is.na(A)] <- 1
  alpha_low <- mean( apply(A, 2, max) > 1 ) # count the points being outside the box (where the box borders are given by the critical values)
  
  # worst case1
  # something like totally linear dependence (in the dimension two)
  j_high <- ceiling( alpha * n ) 
  A <- data_mat / data_order[,n - j_low]
  A[is.na(A)] <- 1
  alpha_high <- mean( apply(A, 2, max) > 1 ) # count the points being outside the box (where the box borders are given by the critical values)
  
  # now we search for values j_low and j_high = j_low + 1 , such that alpha_low <= alpha and alpha_high > alpha
  while( j_high - j_low > 1){
    j_mid <- ceiling( j_low + (j_high - j_low)/2) # approx. middle between j_low and j_high
    A <- data_mat / data_order[,n - j_mid]
    A[is.na(A)] <- 1
    alpha_sim <- mean( apply(A, 2, max) > 1 )
    ifelse(alpha_sim <= alpha, j_low <- j_mid, j_high <- j_mid)
  }
  beta <- 1-(n-j_low)/n # beta_tilde
  
  A <- data_mat / teststat
  A[is.na(A)] <- 1
  
  pvalues <- (rowMeans( A > 1 ))
  
  return(list(pvalues = pvalues, beta = beta))
  
}

global.p.value <- function(data_mat, teststat){
  beta_min <- min(rowMeans(data_mat > teststat))
  quants <- apply(data_mat, 1, quantile, probs = 1-beta_min , type=1)
  mean(apply(data_mat > quants, 2, any))
}

# p.value <- function(data_mat, teststat){
#   data_order <- t( apply(data_mat, 1, sort) )
#   pvalues <- numeric(nrow(data_mat))
#   
#   for(l in 1:length(teststat)){
#     
#     if(any(data_order[l,] == teststat[l])){
#       quants <- data_order[,which(data_order[l,] == teststat[l])[1]]
#     }else{if(any(data_order[l,] >= teststat[l])){
#       high <- which(data_order[l,] >= teststat[l])[1]
#       if(high > 1){
#         quants <- (data_order[,high]+data_order[,high-1])/2 
#       }else{
#         quants <- rep(-Inf, nrow(data_mat))
#       }
#     }else{
#       quants <- rep(Inf, nrow(data_mat))
#       }
#     }
#   pvalues[l] <- mean(apply(data_mat >= quants, 2, any))
#   }
#   pvalues
# }


p.value <- function(data_mat, teststat){
  data_order <- t( apply(data_mat, 1, sort) )
  pvalues <- numeric(nrow(data_mat))
  B <- ncol(data_mat)
  
  for(l in 1:length(teststat)){
    if(teststat[l] < data_order[l,1]){
      pvalues[l] <- 1
    }else{
      beta <- mean(data_mat[l,] > teststat[l])
      if(beta < 1){
        x <- round((1-beta)*B)
        quants <- data_order[,x]
      }else{
        quants <- -Inf
      }
      pvalues[l] <- mean(apply(data_mat > quants, 2, any))
    }
  }
  pvalues
}





wrap_asymptotics_global <- function(global_teststats, global_c_mat_list, crit.value.method){
  
  # p-value calculation
  l <- length(crit.value.method)
  C <- length(global_c_mat_list)
  out <- lapply(1:C, function(c_ind) t(matrix(rep(pchisq(global_teststats[[c_ind]], df=qr(global_c_mat_list[[c_ind]])$rank, lower.tail = FALSE), l),nrow=l,byrow=TRUE)))
  names(out) <- names(global_c_mat_list)
  out
}

wrap_asymptotics_bonf <- function(teststats, c_mat_list, crit.value.method){
  
  l <- length(crit.value.method)
  out <- lapply(1:C, function(c_ind){
    bonf_ergs <- wrap_asymptotics_global(teststats[[c_ind]], c_mat_list[[c_ind]], crit.value.method)
    bonf_mat <- matrix(unlist(bonf_ergs), byrow = TRUE, ncol = l)
    bonf_mat <- apply(bonf_mat, 2, function(p) p.adjust(p, method = "bonferroni"))
    bonf_mat <- as.matrix(bonf_mat)
    colnames(bonf_mat) <- crit.value.method
    rownames(bonf_mat) <- names(c_mat_list[[c_ind]])
    bonf_mat
  })
  
  out
}

wrap_perm <- function(values, tau, global_c_mat_list, Nres, global_teststats, crit.value.method){
  
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
  
  # p values
  lapply(1:C, function(c_ind){
      out_temp <- mean(erg[c_ind,] > global_teststats[[c_ind]][1])
      out <- rep(out_temp,length(crit.value.method))
      
      return(out)
    })
}

wrap_perm_bonf <- function(values, tau, c_mat_list, Nres, teststats, crit.value.method){
  
  l <- length(crit.value.method)
  out <- lapply(1:C, function(c_ind){
    bonf_ergs <- wrap_perm(values, tau, c_mat_list[[c_ind]], Nres, teststats[[c_ind]], crit.value.method)
    bonf_mat <- matrix(unlist(bonf_ergs), byrow = TRUE, ncol = l)
    bonf_mat <- apply(bonf_mat, 2, function(p) p.adjust(p, method = "bonferroni"))
    bonf_mat <- as.matrix(bonf_mat)
    colnames(bonf_mat) <- crit.value.method
    rownames(bonf_mat) <- names(c_mat_list[[c_ind]])
    bonf_mat
  })
  
  out
}

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
  
  # p values
  pvalues <- lapply(1:C, function(ind){
    out <- list()
    
    if("inequi" %in% crit.value.method) out[["inequi"]] <- p.value((t(ts_wildBS[[ind]])) , teststat = teststats[[ind]])
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[["equi"]] <- sapply(teststats[[ind]], function(x) mean(apply(ts_wildBS[[ind]], 1, max) > x))
    
    out
  })
  names(pvalues) <- paste("Matrix", 1:length(pvalues))
  
  return(pvalues)
}


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
  
  # p values
  pvalues <- lapply(1:C, function(ind){
    out <- list()
    
    if("inequi" %in% crit.value.method) out[["inequi"]] <- p.value((t(ts_groupwiseBS[[ind]])) , teststat = teststats[[ind]])
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[["equi"]] <- sapply(teststats[[ind]], function(x) mean(apply(ts_groupwiseBS[[ind]], 1, max) > x))
    
    out
  })
  names(pvalues) <- paste("Matrix", 1:length(pvalues))
  
  return(pvalues)
}


