
### General functions

# function for DA_hat, DN_hat and Y
estimators <- function(X, D, tau, M = NULL, times = NULL, var = TRUE){
  if(is.null(M)) M <- max(D)
  if(is.null(times)) times <- sort(unique(X[D > 0]))
  if(length(times) == 0){
    return(list(Yt = numeric(0), DNt = matrix(ncol=M, nrow=0), 
                DNjt = lapply(1:M,function(x) matrix(ncol=0, nrow=0)), 
                DA_hat = matrix(ncol=M, nrow=0), F_hat = matrix(ncol=M, nrow=0), 
                fm = matrix(ncol=M, nrow=0), gm = matrix(ncol=M, nrow=0), 
                sqrt_fac = numeric(0), ind = numeric(0), diffs = numeric(0)))
  }else{
  Yt <- sapply(times, function(t) sum(X >= t))
  #Nt <- matrix(sapply(1:M, function(m) sapply(times, function(t) sum(X[D == m] <= t))), ncol = M) # M cols, times rows
  DNt <- matrix(0, nrow = length(times), ncol = M)
  # Loop through each m and t combination
  for (m in 1:M) {
    # Find indices where D == m
    indices <- (D == m)
    # Count occurrences of each unique value of X within indices
    counts <- table(X[indices])
    # Update the result matrix with the counts
    for (t_index in 1:length(times)) {
      t <- times[t_index]
      DNt[t_index, m] <- ifelse(t %in% names(counts), counts[[as.character(t)]], 0)
    }
  }
  if(var){ DNjt <- lapply(1:M, function(m){
    if(any(D==m)){ input <- sapply(times, function(t) (X[D == m] == t)) }else{ input <- NA }
    matrix(input,nrow = length(times),ncol=sum(D==m),byrow = TRUE)
  })
  }
  A_hat  <- matrix(apply(ifelse(is.na(DNt / Yt),0,DNt / Yt),2,cumsum),ncol=M)
  DA_hat <- A_hat - rbind(0,A_hat[-length(times),,drop=FALSE])
  
  S_hat  <- cumprod(1 - rowSums(DA_hat))
  S_hat_minus  <- c(1,S_hat[-length(times)])
  F_hat  <- matrix(apply(DA_hat*S_hat_minus,2,cumsum), ncol = M)
  
  ind <- (times <= tau)
  if(!any(ind)){
    return(list(Yt = Yt, DNt = DNt, DNjt = DNjt, DA_hat = DA_hat, F_hat = F_hat, fm = matrix(ncol=M, nrow=0), gm = matrix(ncol=M, nrow=0), 
                sqrt_fac = numeric(0), ind = ind, diffs = numeric(0)))
  }else{
  timestau <- times[ind]
  diffs <- diff(c(timestau,tau))
  
  if(!var){ 
    return(list(Yt = Yt, DNt = DNt, DA_hat = DA_hat, F_hat = F_hat, ind = ind, diffs = diffs))
  }else{
  
  F_hat_diffs <- F_hat[ind,,drop=FALSE]*diffs
  IntF_hat <- matrix(apply(F_hat_diffs,2,function(x) rev(cumsum(rev(x)))),ncol=M)
  #IntF_hat <- t(sapply(timestau, 
  #                     function(t){ 
  #                       ind2 <- (timestau >= t)
  #                       colSums(F_hat_diffs[ind2,,drop=FALSE])
  #                     }))
  fm <- (tau-timestau)*matrix(sapply(1:M, function(m) 1 - rowSums(F_hat[ind,-m,drop=FALSE])),ncol=M) - IntF_hat
  gm <- (tau-timestau)*F_hat[ind,,drop=FALSE] - IntF_hat
  
  sqrt_fac <- numeric(sum(ind))
  nenner <- (1-rowSums(DA_hat[ind,,drop=FALSE]))
  sqrt_fac[nenner > 0]  <- 1/nenner[nenner > 0]
  
  return(list(Yt = Yt, DNt = DNt, DNjt = DNjt, DA_hat = DA_hat, F_hat = F_hat, fm = fm, gm = gm, 
              sqrt_fac = sqrt_fac, ind = ind, diffs = diffs))
  }}}
}

# function for the RMTL
RMTL <- function(X, D, tau, M = NULL, var = TRUE){
  
  # number of competing risks
  if(is.null(M)) M <- max(D)
  # number of observations
  #n <- length(X)
  
  # values of Y, N,.. at t_1,...,t_m
  temp <- estimators(X, D, tau, M, var = var)
  F_hat <- temp$F_hat
  ind <- temp$ind
  diffs <- temp$diffs

  #DF_hat <- F_hat - rbind(0,F_hat[-length(times),])
  mu_hat <- colSums(F_hat[ind,,drop=FALSE] * diffs) ## the RMTL estimator
  
  if(var){
    Yt <- temp$Yt
    DNt <- temp$DNt
    DA_hat <- temp$DA_hat
    
    fm <- temp$fm
    gm <- temp$gm
    fac <- (temp$sqrt_fac)^2
    
  # now: compute the covariance
  Sigma_hat <- matrix(NA, ncol = M, nrow = M)
  sigma_hat <- array(dim = c(sum(ind),M,M))
  
  for(m in 1:M){ # without factor n
    sigma_hat[,m,m] <- cumsum(DA_hat[ind,m] * (1 - DA_hat[ind,m]) / Yt[ind])
    if(m < M){for(m2 in (m+1):M){
      sigma_hat[,m,m2] <- sigma_hat[,m2,m] <- - cumsum(DA_hat[ind,m2] * DA_hat[ind,m] / Yt[ind])
    }
    }
  }
  
  Dsigma_hat <- sigma_hat
  Dsigma_hat[-1,,] <- sigma_hat[-1,,] - sigma_hat[-sum(ind),,]
  
  
  for(m in 1:M){
    Sigma_hat[m,m] <- sum(fac * (Dsigma_hat[,m,m] *  fm[,m]^2 +
      2*rowSums(Dsigma_hat[,m,-m,drop=FALSE]) * gm[,m] * fm[,m] +
      rowSums(Dsigma_hat[,-m,-m,drop=FALSE]) * gm[,m]^2), na.rm = TRUE)
    
    if(m < M){
      for(m2 in (m+1):M){
        Sigma_hat[m,m2] <- 
          Sigma_hat[m2,m] <- sum(fac * (Dsigma_hat[,m,m2] * fm[,m] * fm[,m2] +
          rowSums(Dsigma_hat[,m,-m2,drop=FALSE]) * gm[,m2] * fm[,m] +
          rowSums(Dsigma_hat[,m2,-m,drop=FALSE]) * gm[,m] * fm[,m2] +
          rowSums(Dsigma_hat[,-m,-m2,drop=FALSE]) * gm[,m] * gm[,m2]), na.rm = TRUE)
          
      }}
  }
  return(list("rmtl" = mu_hat, "var_rmtl" = Sigma_hat))
  }else{return(list("rmtl" = mu_hat))}
}

# wrapper for the test statistic for multiple sets of matrices, multiple groups, etc
wrap_multiple_test_stat <- function(values, tau, c_mat_list, c_vec_list, M = NULL, n_group = NULL){
  
  group <- values$group
  if(is.null(M)) M <- max(values$D)
  if(is.null(n_group)) n_group <- max(group)
  kM <- n_group*M
  rmtl <- numeric(kM)
  var <- matrix(0, ncol=kM,nrow=kM)
  
  for(i in 1:n_group) {
    values2 <- values[group == i,]
    #calculate for each group
    temp <- RMTL(values2$X, values2$D, tau, M)
    rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
    var[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
  }
  
  erg_int <- list()
  
  for( hyp in 1:length(c_mat_list)){
    C <- length(c_mat_list[[hyp]])
    erg_int[[hyp]] <- sapply(1:C, function(c_ind){ 
      c_mat <- c_mat_list[[hyp]][[c_ind]]
      vec <- c_mat %*% rmtl - c_vec_list[[hyp]][[c_ind]]
      return(t(vec) %*% MASS::ginv(c_mat %*% var %*% t(c_mat)) %*% vec)
    })
    
  }
  names(erg_int) <- names(c_mat_list)
  return(list(teststats = erg_int, var = var, rmtl=rmtl))
}

# function for getting the rows of a matrix as list
rows2list <- function(mat) apply(mat, 1, function(x)t(x), simplify = FALSE)

# function for multiple critical values
crit_values <- function(data_mat, alpha = 0.05){
  n <- ncol(data_mat)
  dimension <- nrow(data_mat)
  
  # First forget the connection of the coordinates, and sort the values per coordinate
  data_order <- t( apply(data_mat, 1, sort) )
  
  j <- min(c(floor(n*alpha/dimension), n-1))
  j_upper <- n
  #nonDupl <- which(!duplicated(t(data_order))) # instead j+1, nonDupl[nonDupl > j][1] can be used, because there is no change in between (faster for discrete data)
  if(j != j_upper - 1){
  while(mean(apply(data_mat > data_order[,n - j - 1],2,any)) <= alpha){ # check if j is the maximum 
    j <- j+1 # increase j
    if(j == j_upper-1) break # check if j is the maximum
  }}
  
  data_order[,n-j] # critical values
}

# function checks which of the hypotheses is rejected
reject_multiple <- function(teststat, calculated_quantiles) (teststat > calculated_quantiles)

# function for getting the global hypothesis matrix from a list of single hypothesis matrices
global_mat <- function(c_mat, kM) matrix(unlist(c_mat), ncol = kM, byrow=TRUE)


### Functions for asymptotic approaches
# wrapper for the asymptotic approach
wrap_asymptotics <- function(test_stat_list, c_mat_list, kM, crit.value.method, Nres, ginvs, alpha = 0.05){
  # extract Sigma_hat
  Sigma_hat <- test_stat_list$var
  C <- length(c_mat_list)
  
  # create output object
  # save the test decision
  l <- length(crit.value.method)
  out <- vector(mode = "list", C)
  names(out) <- names(c_mat_list)
  
  randomNumbers_needed <- TRUE
  # go through the global hypotheses
  for(c_ind in 1:C){
    # if all matrices are row vectors
    if(all(sapply(c_mat_list[[c_ind]], function(x) nrow(x)) == 1) ){
      # calculate the scaling matrix
      D <- diag(lapply(c_mat_list[[c_ind]], function(vec){
        MASS::ginv(sqrt(vec %*% Sigma_hat %*% t(vec))) #using the g-inverse to avoid problems with zero variances
      }))
      H <- global_mat(c_mat_list[[c_ind]], kM)
      # equicoordinate normal quantiles
      sigma <- (D%*%H%*%Sigma_hat%*%t(H)%*%D)
      # since errors occur when any of diag(sigma) = 0, we do not consider these components
      my_index <- (diag(sigma) != 0)
      if(any(my_index)){
        quant <- (mvtnorm::qmvnorm(1 - alpha , tail="both.tails", sigma=sigma[my_index, my_index])$quantile)^2
      }else{ quant <- 0 }
      
      
      out_temp <- t(matrix(rep(test_stat_list$teststats[[c_ind]] > quant,l), nrow=l, byrow=TRUE))
      colnames(out_temp) <- crit.value.method
      rownames(out_temp) <- names(test_stat_list$teststats[[c_ind]])
      
      out[[c_ind]] <- out_temp
      
    }else{
      
      # generate Nres random numbers if we do not have them yet
      if(randomNumbers_needed){
        random_numbers <- mvtnorm::rmvnorm(Nres, sigma = Sigma_hat)
        randomNumbers_needed <- FALSE
      } 
      
      
      hyp_mat <- c_mat_list[[c_ind]]
      # determine values for approximating the limiting distribution
      random_values <- t(sapply(1:length(hyp_mat), function(mat_ind) apply(random_numbers, 1, function(z){
        vec <- hyp_mat[[mat_ind]]%*%z
        t(vec)%*%ginvs[[c_ind]][[mat_ind]]%*%vec
      } )))
      
      
      out[[c_ind]] <- matrix(logical(l*length(test_stat_list$teststats[[c_ind]])),ncol=l)
      colnames(out[[c_ind]]) <- crit.value.method
      
      if("inequi" %in% crit.value.method) out[[c_ind]][,"inequi"] <- reject_multiple(test_stat_list$teststats[[c_ind]], crit_values(random_values, alpha = alpha))
      #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
      if("equi" %in% crit.value.method) out[[c_ind]][,"equi"] <- reject_multiple(test_stat_list$teststats[[c_ind]], quantile(apply(random_values, 2, max), probs = 1 - alpha))
      
    }
  }
  return(out)
}  

# wrapper for the global asymptotic approach (unadjusted)
wrap_asymptotics_global <- function(global_teststats, global_c_mat_list, crit.value.method, alpha = 0.05){
  
  C2 <- length(global_c_mat_list)
  # calculate the quantiles
  quant <- sapply(global_c_mat_list, function(H) qchisq(1 - alpha,df=qr(H)$rank))
  # test decision
  l <- length(crit.value.method)
  out <- lapply(1:C2, function(c_ind) t(matrix(rep(any(global_teststats[c_ind] > quant[c_ind]), l),nrow=l,byrow=TRUE)))
  names(out) <- names(global_c_mat_list)
  out
}

# wrapper for the asymptotic approach for multiple hypotheses with bonferroni correction
wrap_asymptotics_bonf <- function(teststats, c_mat_list, kM, crit.value.method, alpha = 0.05){
  C <- length(c_mat_list)
  l <- length(crit.value.method)
  out <- lapply(1:C, function(c_ind){
    bonf_ergs <- wrap_asymptotics_global(teststats[[c_ind]], c_mat_list[[c_ind]], crit.value.method, alpha = alpha/length(c_mat_list[[c_ind]]))
    bonf_mat <- matrix(unlist(bonf_ergs), byrow = TRUE, ncol = l)
    colnames(bonf_mat) <- crit.value.method
    rownames(bonf_mat) <- names(c_mat_list[[c_ind]])
    bonf_mat
  })
  
  out
}


### Functions for pooled bootstrap
# function for calculating the test statistic for the pooled Bootstrap
test_stat_pooledBS <- function(bs_values, tau, sqrt_Sigma_hat, mu_hat, n_group, M, kM, noinv = FALSE){
  group <- bs_values$group
  rmtl <- numeric(kM)
  var <- matrix(0, ncol=kM,nrow=kM)
  sqrt_var_inv <- matrix(0, ncol=kM,nrow=kM)
  
  for(i in 1:n_group) {
    values2 <- bs_values[group == i,]
    #calculate for each group
    temp <- RMTL(values2$X, values2$D, tau, M = M)
    rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
    var[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
    if(!noinv){
    # eigen value decomposition of var
    eigen_temp <- eigen(temp[["var_rmtl"]])
    Mitte <- matrix(0,nrow=M,ncol=M)
    diag(Mitte)[eigen_temp$values > 0] <- ((eigen_temp$values)[eigen_temp$values > 0])^(-1/2)
    # determine the inverse root of var
    sqrt_var_inv[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- eigen_temp$vectors%*%Mitte%*%t(eigen_temp$vectors)
  }}
  
  if(!noinv){return(list(vec1 = (rmtl- mu_hat), mat1 = sqrt_Sigma_hat %*% sqrt_var_inv, mat2 = var ))}
  if(noinv){return(list(vec1 = (rmtl- mu_hat), mat1 = var ))}
  #calculate the pooled BS test statistic
  #vec <- c_mat %*% sqrt_Sigma_hat %*% sqrt_var_inv %*% (rmtl- mu_hat)
  #return( t(vec) %*% MASS::ginv(c_mat %*% Sigma_hat %*% t(c_mat)) %*% vec )  
  # the factor n is eliminated by the other factors
}

# wrapper for both pooled bootstraps
wrap_pooledBS <- function(values, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, sqrt_Sigma_hat, ginvs, alpha = 0.05){
  
  # calculate Sigma_hat and the test statistic
  group <- values$group
  n_group <- max(group)
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)
  
  Sigma_hat <- test_stat_list$var
  teststats <- test_stat_list$teststats
  
  #calculate mu_hat for the pooled sample
  temp <- RMTL(values$X, values$D, tau, M = M, var = FALSE)
  mu_hat <- temp[["rmtl"]]
  
  # generate pooled Bootstrap samples
  bs_ergs <- replicate(Nres,{
    bs_values <- values[sample(1:n,n,replace = TRUE),c("X", "D")]
    bs_values$group <- group
    test_stat_pooledBS(bs_values=bs_values, tau=tau, 
                       sqrt_Sigma_hat = sqrt_Sigma_hat,
                       mu_hat=mu_hat, n_group = n_group, M = M, kM = kM, noinv = FALSE)
  }, simplify = FALSE)
  
  ## first test statistic
  erg <- sapply(bs_ergs, function(temp1){
    # calculate the two pooled BS test statistics for all matrices in c_mat_list
    erg_int <- list()
    for(hyp in 1:C){
      erg_int[[hyp]] <- sapply(1:length(c_mat_list[[hyp]]), function(c_mat_ind){ 
        vec <- c_mat_list[[hyp]][[c_mat_ind]] %*% temp1$mat1 %*% temp1$vec1
        return( t(vec) %*% ginvs[[hyp]][[c_mat_ind]] %*% vec )  
      })
    }
    names(erg_int) <- namen
    erg_int
  })
  dim(erg) <- c(C,Nres)
  
  # save the test statistics in arrays in a list
  ts_pooledBS <- list()
  for( hyp in 1:C){
    dim <- length(erg[hyp,][[1]])
    ts_pooledBS[[hyp]] <- matrix(unlist(erg[hyp,]), ncol = dim, byrow=TRUE)
  }
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out <- matrix(logical(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- reject_multiple(teststats[[ind]], crit_values(t(ts_pooledBS[[ind]]), alpha = alpha))
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- reject_multiple(teststats[[ind]], quantile(apply(ts_pooledBS[[ind]], 1, max), probs = 1 - alpha))
    
    out
  })
  
  ## second test statistic
  erg <- sapply(bs_ergs, function(temp1){
    # calculate the two pooled BS test statistics for all matrices in c_mat_list
    erg_int <- list()
    for(hyp in 1:C){
      erg_int[[hyp]] <- sapply(1:length(c_mat_list[[hyp]]), function(c_mat_ind){ 
        vec <- c_mat_list[[hyp]][[c_mat_ind]] %*% temp1$vec1
        return( t(vec) %*% MASS::ginv(c_mat_list[[hyp]][[c_mat_ind]] %*% temp1$mat2 %*%t(c_mat_list[[hyp]][[c_mat_ind]])) %*% vec )  
      })
    }
    names(erg_int) <- namen
    erg_int
  })
  dim(erg) <- c(C,Nres)
  
  # save the test statistics in arrays in a list
  ts_pooledBS <- list()
  for( hyp in 1:C){
    dim <- length(erg[hyp,][[1]])
    ts_pooledBS[[hyp]] <- matrix(unlist(erg[hyp,]), ncol = dim, byrow=TRUE)
  }
  
  # test decisions
  decision2 <- sapply(1:C, function(ind){
    out <- matrix(logical(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    out_temp <- reject_multiple(teststats[[ind]], apply(ts_pooledBS[[ind]], 2, function(x) quantile(x, prob = 1 - alpha/length(c_mat_list[[ind]]), type = 1)))
    if("inequi" %in% crit.value.method) out[,"inequi"] <- out_temp
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- out_temp
    
    out
  })
  return(list(decision,decision2))
}

### Functions for wild bootstrap
mammen <- 0.5*(1 + c(-1,1)*sqrt(5)); probs_mammen <- 0.5*(1 + c(1,-1)/sqrt(5))
# function for calculating the test statistic for the wild Bootstrap
test_stat_wildBS <- function(est_list, bs_weights, tau, n_group, M, kM){
  
  wild_rmtl <- numeric(kM)
  wild_var <- matrix(0,ncol=kM,nrow=kM)
  
  for(i in 1:n_group) {
    Yt <- est_list[[i]]$Yt
    DNt <- est_list[[i]]$DNt
    DNjt <- est_list[[i]]$DNjt
    DA_hat <- est_list[[i]]$DA_hat
    F_hat <- est_list[[i]]$F_hat
    fm <- est_list[[i]]$fm
    gm <- est_list[[i]]$gm
    sqrt_fac <- (est_list[[i]]$sqrt_fac)
    ind <- est_list[[i]]$ind
    nind <- sum(ind)
    diffs <- est_list[[i]]$diffs
    invsqrt2 <- 1/sqrt(2)
    fac <- sqrt_fac^2
    my_fac <- sqrt(1-rowSums(DA_hat[ind,,drop=FALSE]))
    my_fac <- ifelse(is.na(my_fac),0,my_fac)
    sqrt_DA_hat <- sqrt(DA_hat[ind,,drop=FALSE])
               
    weight <- bs_weights[[i]]
    
    DW_hat <- matrix(NA,ncol=M,nrow=nind)
    bs_mults <- bs_mults2 <- array(dim = c(M,M,nind))
    for(m in 1:M){
      bs_mults[m,,] <- t(sapply(1:M, function(l) DNjt[[l]][ind,,drop=FALSE] %*% weight[[m]][[l]])) 
      bs_mults2[m,,] <- t(sapply(1:M, function(l) DNjt[[l]][ind,,drop=FALSE] %*% (weight[[m]][[l]])^2 )) 
    }
    
    for(m in 1:M){
     DW_hat[,m] <- (bs_mults[m,m,] * my_fac + 
       invsqrt2 * rowSums(sapply(1:M, function(l) sign(l - m) * 
                                   (bs_mults[m,l,] * sqrt_DA_hat[,m] + 
                                       bs_mults[l,m,] * sqrt_DA_hat[,l]))))/Yt[ind]
    }
    DW_hat <- ifelse(is.na(DW_hat),0,DW_hat)
    
      
    wild_rmtl[((i-1)*M +1):(i*M)] <- colSums((fm * DW_hat + gm * 
                               sapply(1:M, function(m) rowSums(DW_hat[,-m,drop=FALSE]))) * 
                              sqrt_fac, na.rm = TRUE)
    
    Sigma_hat <- matrix(NA, ncol = M, nrow = M)
    Dsigma_hat <- array(dim = c(nind,M,M))
    
    for(m in 1:M){ # without factor n
      Dsigma_hat[,m,m] <- (bs_mults2[m,m,] * my_fac + ### ab hier weiter # DNjt ist noch in keinem schönen Format...
                             0.5 * rowSums(matrix(sapply(1:M, function(l) sign(abs(l-m)) * (bs_mults2[m,l,] * sqrt_DA_hat[,m] + 
                                                            bs_mults2[l,m,] * sqrt_DA_hat[,l])), nrow = nind)))/(Yt[ind])^2
      if(m < M){for(m2 in (m+1):M){
        Dsigma_hat[,m,m2] <- Dsigma_hat[,m2,m] <- -0.5* (bs_mults2[m,m2,] * sqrt_DA_hat[,m] + 
                                                          bs_mults2[m2,m,] * sqrt_DA_hat[,m2]) / (Yt[ind])^2
      }
      }
    }
    
    
    for(m in 1:M){
      Sigma_hat[m,m] <- sum(fac * (Dsigma_hat[,m,m] *  fm[,m]^2 +
                                     2*rowSums(Dsigma_hat[,m,-m,drop=FALSE]) * gm[,m] * fm[,m] +
                                     rowSums(Dsigma_hat[,-m,-m,drop=FALSE]) * gm[,m]^2), na.rm = TRUE)
      
      if(m < M){
        for(m2 in (m+1):M){
          Sigma_hat[m,m2] <- 
            Sigma_hat[m2,m] <- sum(fac * (Dsigma_hat[,m,m2] * fm[,m] * fm[,m2] +
                                            rowSums(Dsigma_hat[,m,-m2,drop=FALSE]) * gm[,m2] * fm[,m] +
                                            rowSums(Dsigma_hat[,m2,-m,drop=FALSE]) * gm[,m] * fm[,m2] +
                                            rowSums(Dsigma_hat[,-m,-m2,drop=FALSE]) * gm[,m] * gm[,m2]), na.rm = TRUE)
          
        }}
    }
    
    wild_var[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- Sigma_hat
    
  }
  
  #calculate some vectors and matrices needed more than once
  #vec <- c_mat %*% wild_rmtl
  
  #calculate the pooled BS test statistic
  #return( t(vec) %*% MASS::ginv(c_mat %*% (wild_var) %*% t(c_mat)) %*% vec )  
  # the factor n is eliminated by the other factors
  # appearing in the formula for the test statistic and
  # the factor n_i in the variance formula
  return(list(wild_rmtl = wild_rmtl, wild_var = wild_var ))
}

# wrapper for the wild bootstrap
wrap_wildBS <- function(values, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, est_list, bs_function, alpha = 0.05){
  
  # calculate the test statistic
  group <- values$group
  n_group <- max(group)
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)
  
  teststats <- test_stat_list$teststats
  
  erg <- sapply(1:Nres, function(huhu){
    # generate wild bootstrap weights
    bs_weights <- lapply(1:n_group, function(x) 
      lapply(1:M, function(y) 
        sapply(1:M, function(m) bs_function(ncol((est_list[[x]])$DNjt[[m]])), 
               simplify = FALSE)))
    
    # calculate the rmtl and var for wild bs
    RmtlAndVar <- test_stat_wildBS(est_list = est_list, bs_weights = bs_weights, 
                                   tau = tau, n_group = n_group, M = M, kM = kM)
    erg_int <- list()
    for(hyp in 1:C){
      erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
        vec <- c_mat %*% RmtlAndVar$wild_rmtl
        
        #calculate the pooled BS test statistic
        return( t(vec) %*% MASS::ginv(c_mat %*% (RmtlAndVar$wild_var) %*% t(c_mat)) %*% vec )  
      })
    }
    names(erg_int) <- namen
    erg_int
  })
  dim(erg) <- c(C,Nres)
  
  # save the test statistics in arrays in a list
  ts_wildBS <- list()
  for( hyp in 1:C){
    dim <- length(erg[hyp,][[1]])
    ts_wildBS[[hyp]] <- matrix(unlist(erg[hyp,]), ncol = dim, byrow=TRUE)
  }
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out <- matrix(logical(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- reject_multiple(teststats[[ind]], crit_values(t(ts_wildBS[[ind]]), alpha = alpha))
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- reject_multiple(teststats[[ind]], quantile(apply(ts_wildBS[[ind]], 1, max), probs = 1 - alpha))
    
    out
  })
  return(decision)
}

### Functions for groupwise bootstrap
# wrapper for the groupwise bootstrap # the test statistics are calculated in the wrapper
wrap_groupwiseBS <- function(values, tau, c_mat_list, Nres, test_stat_list, n_vec, crit.value.method, M, kM, alpha = 0.05){
  
  # calculate the the test statistic
  group <- values$group
  n_group <- max(group)
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)

  teststats <- test_stat_list$teststats
  rmtl <- test_stat_list$rmtl
  
  erg <- replicate(Nres , expr = {

    # calculate the rmst and variance estimator for the bootstrap sample
    bs_rmtl <- numeric(kM)
    bs_var <- matrix(0,ncol=kM,nrow=kM)
    
    for(i in 1:n_group) {
      # generate a bootstrap sample in group i
      values2 <- values[group == i,][sample(n_vec[i],n_vec[i],replace = TRUE),]
      #calculate for each group
      temp <- RMTL(values2$X, values2$D, tau, M = M)
      bs_rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
      bs_var[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
    }
    
    # calculate the groupwise BS test statistics for all matrices in c_mat_list
    erg_int <- list()
    for( hyp in 1:C){
      erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
        vec <- c_mat %*% (bs_rmtl- rmtl)
        t(vec) %*% MASS::ginv(c_mat %*% (bs_var) %*% t(c_mat)) %*% vec  
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
    out <- matrix(logical(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- reject_multiple(teststats[[ind]], crit_values(t(ts_groupwiseBS[[ind]]), alpha = alpha))
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- reject_multiple(teststats[[ind]], quantile(apply(ts_groupwiseBS[[ind]], 1, max), probs = 1 - alpha))
    
    out
  })
  
  return(decision)
}

### Functions for studentized permutation
# function for the parts of the studentized permutation test statistic
test_stat_perm_bonf <- function(perm_values, tau, mu_hat, n_group, M, kM){
  group <- perm_values$group
  rmtl <- numeric(kM)
  var <- matrix(0, ncol=kM,nrow=kM)
  
  for(i in 1:n_group) {
    values2 <- perm_values[group == i,]
    #calculate for each group
    temp <- RMTL(values2$X, values2$D, tau, M = M)
    rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
    var[((i-1)*M +1):(i*M), ((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
  }

  return(list(vec1 = (rmtl- mu_hat), var = var))
  #calculate the pooled BS test statistic
  #vec <- c_mat %*% (rmtl- mu_hat)
  #return( t(vec) %*% MASS::ginv(c_mat %*% Sigma_hat %*% t(c_mat)) %*% vec )  
  # the factor n is eliminated by the other factors
}

# wrapper for the studentized permutation
wrap_perm_bonf <- function(values, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, alpha = 0.05){
  
  # calculate Sigma_hat and the test statistic
  group <- values$group
  n_group <- max(group)
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)
  
  teststats <- test_stat_list$teststats
  
  #calculate mu_hat for the pooled sample
  temp <- RMTL(values$X, values$D, tau, M = M, var = FALSE)
  mu_hat <- temp[["rmtl"]]
  
  # generate permuted samples
  perm_samples <- replicate(Nres,{
    perm_values <- values[,c("X", "D")]
    perm_values$group <- group[sample(n,n,replace = FALSE)]
    perm_values
  }, simplify = FALSE)
  
  erg <- sapply(perm_samples, function(perm_values){
    perms <- test_stat_perm_bonf(perm_values=perm_values, tau=tau, # wie beim pooled BS
                                    mu_hat=mu_hat, n_group = n_group, M = M, kM = kM)
    # calculate the two pooled BS test statistics for all matrices in c_mat_list
    erg_int <- list()
    for(hyp in 1:C){
      erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
        vec <- c_mat %*% perms$vec1
        return( t(vec) %*% MASS::ginv(c_mat %*% perms$var %*% t(c_mat)) %*% vec )  
      })
    }
    names(erg_int) <- namen
    erg_int
  })
  dim(erg) <- c(C,Nres)
  
  # save the test statistics in arrays in a list
  ts_pooledBS <- list()
  for( hyp in 1:C){
    dim <- length(erg[hyp,][[1]])
    ts_pooledBS[[hyp]] <- matrix(unlist(erg[hyp,]), ncol = dim, byrow=TRUE)
  }
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out <- matrix(logical(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    out_temp <- reject_multiple(teststats[[ind]], apply(ts_pooledBS[[ind]], 2, function(x) quantile(x, prob = 1 - alpha/length(c_mat_list[[ind]]), type = 1)))
    if("inequi" %in% crit.value.method) out[,"inequi"] <- out_temp
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- out_temp
    
    out
  })
  
  return(decision)
}

### Functions for randomization
# wrapper for randomization
wrap_rand <- function(values, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, ginvs, alpha = 0.05){
  
  # calculate Sigma_hat and the test statistic
  group <- values$group
  n_group <- max(group)
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)
  
  teststats <- test_stat_list$teststats
  rmtl <- test_stat_list$rmtl
  
  #calculate the centering mu_hat 
  mu_hat <- rep(colMeans(matrix(rmtl,ncol=n_group)),each=M)
  
  # generate randomized samples
  rand_ergs <- replicate(Nres,{
    rand_values <- values
    my_ind <- (rand_values$D > 0)
    rand_values$D[my_ind] <- sample(M,sum(my_ind),replace=TRUE)
#    test_stat_rand(rand_values=rand_values, tau=tau,
#                       mu_hat=mu_hat, n_group = n_group, M = M, kM = kM)
    test_stat_pooledBS(rand_values, tau=tau, sqrt_Sigma_hat, mu_hat=mu_hat, n_group = n_group, M = M, kM = kM, noinv = TRUE)
  }, simplify = FALSE)
  

  erg <- sapply(rand_ergs, function(temp1){
    # calculate the randomized test statistics for all matrices in c_mat_list
    erg_int <- list()
    for(hyp in 1:C){
      erg_int[[hyp]] <- sapply(1:length(c_mat_list[[hyp]]), function(c_mat_ind){ 
        vec <- c_mat_list[[hyp]][[c_mat_ind]] %*% temp1$vec1
        return( t(vec) %*% MASS::ginv(c_mat_list[[hyp]][[c_mat_ind]] %*% temp1$mat1 %*% t(c_mat_list[[hyp]][[c_mat_ind]])) %*% vec )  
      })
    }
    names(erg_int) <- namen
    erg_int
  })
  dim(erg) <- c(C,Nres)
  
  # save the test statistics in arrays in a list
  ts_rand <- list()
  for( hyp in 1:C){
    dim <- length(erg[hyp,][[1]])
    ts_rand[[hyp]] <- matrix(unlist(erg[hyp,]), ncol = dim, byrow=TRUE)
  }
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out <- matrix(logical(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    out_temp <- reject_multiple(teststats[[ind]], apply(ts_rand[[ind]], 2, function(x) quantile(x, prob = 1 - alpha/length(c_mat_list[[ind]]), type = 1)))
    if("inequi" %in% crit.value.method) out[,"inequi"] <- out_temp
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- out_temp
    
    out
  })
  return(decision)
}

### Functions for p-value permutation
# calculate the p value
p.value <- function(data_mat, teststat){
  data_order <- t( apply(data_mat, 1, sort) )
  pvalues <- numeric(nrow(data_mat))
  B <- ncol(data_mat)
  
  for(l in 1:length(teststat)){
    if(teststat[l] < data_order[l,1]){
      pvalues[l] <- 1
    }else{
      beta <- mean(data_mat[l,] >= teststat[l])
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

# function for p-values of a sample
pvals_gwBS <- function(values, bs_indices, Nres, kM, M, tau, C, c_mat_list, n_group, group){
  
  rmtl <- numeric(kM)
  var <- matrix(0, ncol=kM,nrow=kM)
  
  for(i in 1:n_group) {
    values2 <- values[group == i,]
    #calculate for each group
    temp <- RMTL(values2$X, values2$D, tau, M)
    rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
    var[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
  }
  
  erg_int <- list()
  
  for( hyp in 1:C){
    C2 <- length(c_mat_list[[hyp]])
    erg_int[[hyp]] <- sapply(1:C2, function(c_ind){ 
      c_mat <- c_mat_list[[hyp]][[c_ind]]
      vec <- c_mat %*% rmtl # - c_vec_list[[hyp]][[c_ind]]
      return(t(vec) %*% MASS::ginv(c_mat %*% var %*% t(c_mat)) %*% vec)
    })
    
  }
  names(erg_int) <- names(c_mat_list)
  teststats <- erg_int
  
  # groupwise Bootstrap
  erg <- sapply(1:Nres , function(b){
    
    # calculate the rmst and variance estimator for the bootstrap sample
    bs_rmtl <- numeric(kM)
    bs_var <- matrix(0,ncol=kM,nrow=kM)
    
    for(i in 1:n_group) {
      # generate a bootstrap sample in group i
      values2 <- values[group == i,][bs_indices[[i,b]],]
      #calculate for each group
      temp <- RMTL(values2$X, values2$D, tau, M = M)
      bs_rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
      bs_var[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
    }
    
    # calculate the groupwise BS test statistics for all matrices in c_mat_list
    erg_int <- list()
    for( hyp in 1:C){
      erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
        vec <- c_mat %*% (bs_rmtl- rmtl)
        t(vec) %*% MASS::ginv(c_mat %*% (bs_var) %*% t(c_mat)) %*% vec  
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
  # p-values
  pvals <- sapply(1:C, function(ind){
    out <- matrix(logical(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- p.value(t(ts_groupwiseBS[[ind]]), teststats[[ind]])
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- sapply(teststats[[ind]], function(x) mean(x <= apply(ts_groupwiseBS[[ind]], 1, max)))
    
    out
  }, simplify = FALSE)
  pvals
}


# wrapper for the p-value permutation
wrap_perm <- function(values, tau, c_mat_list, Nres, test_stat_list, crit.value.method, M, kM, alpha = 0.05){
  
  # calculate the the test statistic
  group <- values$group
  n_group <- max(group)
  n <- nrow(values)
  C <- length(c_mat_list)
  namen <- names(c_mat_list)
  
  teststats <- test_stat_list$teststats
  rmtl <- test_stat_list$rmtl
  
  # groupwise bootstrap
  bs_indices <- replicate(Nres, expr = {
    index <- vector(mode = "list", length = n_group)
    for(i in 1:n_group) {
      index[[i]] <- sample(n_vec[i],n_vec[i],replace = TRUE)
    }
    index
  })
  
  # groupwise Bootstrap
  erg <- sapply(1:Nres , function(b){
    
    # calculate the rmst and variance estimator for the bootstrap sample
    bs_rmtl <- numeric(kM)
    bs_var <- matrix(0,ncol=kM,nrow=kM)
    
    for(i in 1:n_group) {
      # generate a bootstrap sample in group i
      values2 <- values[group == i,][bs_indices[[i,b]],]
      #calculate for each group
      temp <- RMTL(values2$X, values2$D, tau, M = M)
      bs_rmtl[((i-1)*M +1):(i*M)] <- temp[["rmtl"]]
      bs_var[((i-1)*M +1):(i*M),((i-1)*M +1):(i*M)] <- temp[["var_rmtl"]]
    }
    
    # calculate the groupwise BS test statistics for all matrices in c_mat_list
    erg_int <- list()
    for( hyp in 1:C){
      erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
        vec <- c_mat %*% (bs_rmtl- rmtl)
        t(vec) %*% MASS::ginv(c_mat %*% (bs_var) %*% t(c_mat)) %*% vec  
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
  # p-values
  pvals <- sapply(1:C, function(ind){
    out <- matrix(logical(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- p.value(t(ts_groupwiseBS[[ind]]), teststats[[ind]])
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- sapply(teststats[[ind]], function(x) mean(x <= apply(ts_groupwiseBS[[ind]], 1, max)))
    
    out
  }, simplify = FALSE)
  
  # generate permuted samples
  perm_samples <- replicate(Nperm,{
    perm_values <- values[,c("X", "D")]
    perm_values$group <- group[sample(n,n,replace = FALSE)]
    perm_values
  }, simplify = FALSE)
  
  pvals_perm <- sapply(perm_samples, function(perm_values){
    lapply(pvals_gwBS(perm_values, bs_indices, Nres, kM, M, tau, 
                      C, c_mat_list, n_group, perm_values$group), function(x) apply(x,2,min))
  })
  dim(pvals_perm) <- c(C,Nperm)
  
  # save the test statistics in arrays in a list
  ts_perm <- list()
  for( hyp in 1:C){
    dim <- length(pvals_perm[hyp,][[1]])
    ts_perm[[hyp]] <- matrix(unlist(pvals_perm[hyp,]), ncol = dim, byrow=TRUE)
    colnames(ts_perm[[hyp]]) <- names(pvals_perm[hyp,][[1]])
  }
  
  # test decisions
  decision <- sapply(1:C, function(ind){
    out <- matrix(logical(length(crit.value.method)*length(teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- (pvals[[ind]][,"inequi"] < quantile(ts_perm[[ind]][,"inequi"], probs = alpha,type=1))
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- (pvals[[ind]][,"equi"] < quantile(ts_perm[[ind]][,"equi"], probs = alpha,type=1))
    
    out
  })
  
  return(decision)
}



