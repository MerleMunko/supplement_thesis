
### useful functions for the Simulation


#Sort object of data_gen
#Input: 
# values:   matrix - created by data_gen

sort_data <- function(values) values[order(values[, 1]), ]

# ## a faster function than table(), which is not really faster...
# # x: numeric vector of survival times
# # y: numeric vector of 0 (= censored) and 1 (= uncensored)
# # output
# # table of x, y 
# fast_table <- function (x,y){
#   x_uniq <- unique(x)
#   out <- matrix(NA, ncol=2, nrow=length(x_uniq))
#   colnames(out) <- c(0,1)
#   rownames(out) <- x_uniq
# 
#   censored <- (y==0)
#   x_cens <- x[censored]
#   x_uncens <- x[!censored]
#   out[,1] <- sapply(x_uniq, function(z) sum(x_cens==z))
#   out[,2] <- sapply(x_uniq, function(z) sum(x_uncens==z))
#   out <- out[order(x_uniq),]
# 
#   out
# }

# function: simple_surv() determines basic survival quantities in a faster 
#           manner than survfit()
# input: time - survival time
#        cens - indicator for censoring (1=not censored, 0=censored)
# output: matrix of dimensions (n, 4) to speed up later calculations

simple_surv <- function(time, cens) {
  n.all <- length(time)
  tab <- table(time, factor(cens, levels = c(0,1)))
  d <- tab[,1] # number of withdrawals in time points
  w <- tab[,2] # number of events in time points
  n <- c(n.all, n.all - cumsum(d + w)) # calculate risk set
  n <- n[-length(n)]
  s <- cumprod(1 - (w /  n)) # calculate Kaplan-Meier-Estimator
  cbind(as.numeric(row.names(tab)), s, w, n)
}

# function: RMST() estimates Restricted Mean Survival Time and its variance 
#           for a sample of tupel (time, cens), given a time point tau
# input: time - survival time
#        cens - indicator for censoring (1=censored, 0=not censored)
#        tau - restriction time point
# output: rmst - estimated restricted mean survival time for tau
#         var_rmst - estimated variance of the rmst

RMST <- function(time, cens, tau){
  #n <- length(time) # number of observation in the beginning of study
  survtab <- simple_surv(time, cens) # fit convenient model for quantities
  
  t_i <- survtab[,1] <= tau # identify time points t_i <= tau
  t_sel <- survtab[,1][t_i] # select relevent time points
  S_km <- survtab[,2][t_i] # calculate Kaplan-Meier estimator
  
  w.factor <- diff(c(0, t_sel,tau)) # width of the area under S(t) from t_0=0
  rmst <- sum(w.factor * c(1, S_km)) # calculate area under S(t) up to t_i=tau
  
  area <- rev(cumsum(rev(diff(c(t_sel,tau)) * S_km)))^2 # calculate areas under S(t) from t_i to tau for
  d <- survtab[,3][t_i] # determine number of events
  Y <- survtab[,4][t_i] # determine number of individuals under risk
  var_rmst <- sum(( d/Y/(Y-d)) * area, na.rm = TRUE) # calculate final variance without 
                                                      # the factor n_i, because this factor 
                                                      # is eliminated by the factor in front 
                                                      # of the variances (i.e. n/n_i) and in 
                                                      # front of the test statistic (i.e. n^{1/2})
                                                      # We add na.rm = TRUE to the sum to easily handle
                                                      # the case Y = d (in the case the last summand is 0/0)
  c(rmst = rmst, var_rmst = var_rmst)
}

#wrapper for RMST several groups in data
#Input:
# values    - matrix created by dat_gen, sorted by sort_data (?)
#Output:    Standard deviation estimated with interval formula for each group

test_stat <- function(values, tau, c_mat) {
  group <- values$group
  n_group <- max(group)
  rmst <- numeric(n_group)
  var <- numeric(n_group)
  
  for(i in 1:n_group) {
    values2 <- values[group == i,]
    #calculate for each group
    temp <- RMST(values2$time, values2$status, tau)
    rmst[i] <- temp["rmst"]
    var[i] <- temp["var_rmst"]
  }
  return( t(c_mat %*% rmst) %*% 
            MASS::ginv(c_mat %*% diag(var) %*% t(c_mat)) %*% 
            (c_mat %*% rmst))     # the factor n is eliminated by the other factors
  # appearing in the formula for the test statistic and
  # the factor n_i in the variance formula
}





#wrapper to get the multiple test statistics directly from the data 
#for potentially more than one partitionized hypothesis matrix
#Input
# values:   data with columns time, status and group
#Output
# list of (multiple) test statistics and of the variance estimators

wrap_multiple_test_stat <- function(values, c_mat_list, tau){
  
  group <- values$group
  n_group <- max(group)
  rmst <- numeric(n_group)
  var <- numeric(n_group)
  
  for(i in 1:n_group) {
    values2 <- values[group == i,]
    #calculate for each group
    temp <- RMST(values2$time, values2$status, tau)
    rmst[i] <- temp["rmst"]
    var[i] <- temp["var_rmst"]
  }
  
  erg_int <- list()
  
  for( hyp in 1:length(c_mat_list)){
    erg_int[[hyp]] <- sapply(c_mat_list[[hyp]], function(c_mat){ 
      vec <- c_mat %*% rmst
      return(t(vec) %*% MASS::ginv(c_mat %*% diag(var) %*% t(c_mat)) %*% vec)
      })
    
  }
  names(erg_int) <- names(c_mat_list)
  return(list(teststats = erg_int, var = var, rmst=rmst))
}


# function for getting the rows of a matrix as list
rows2list <- function(mat) apply(mat, 1, function(x)t(x), simplify = FALSE)


# function for calculating the multiple quantiles
calc_quantiles <- function(points, alpha = 0.05){
  
  # find \tilde{beta} firstly
  func <- function(beta){
    q <- apply(points, 2, function(x) quantile(x, probs = 1-beta,type=1))
    mean(apply(t(points) > q, 2, any))
  }
  
  # vector of possible \tilde{beta}s
  B <- nrow(points)
  beta <- (0:(B-1))/B
  beta <- beta[(beta >= floor(B*alpha/ncol(points))/B)  & (beta <= alpha)]
  
  laenge <- length(beta)
  i <- 0
  
  while(func(beta[i+1]) <= alpha){
    i <- i + 1
    if(i+1 > laenge) break
  }
  
  tilde_beta <- beta[i]
  
  return(apply(points, 2, function(x) quantile(x, probs = 1-tilde_beta,type=1)))
  
}


# half partition search, quicker for higher degree of dependence
crit_values2 <- function(data_mat, alpha = 0.05){
  # the input is a matrix data_mat of dimenson n_obs x dim_obs containing the B (resampling) observations X_b of dimension dim_obs 
  # for our specific purpose X_i is the vector of the test statistics based on the different contrast vectors, i.e. dim_obs = r in the pdf
  # each column represent one observation (i.e. one resampling iteration)
  # the output are the critical values for each dimension (large X_b lead to a rejection)
  # such that the overall level is alpha
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
  data_order[,n - j_low ] # critical values
}

# function checks whether at least one hypothesis is rejected
reject <- function(teststat, calculated_quantiles) any(teststat > calculated_quantiles)

# function checks which of the hypotheses is rejected
reject_multiple <- function(teststat, calculated_quantiles) (teststat > calculated_quantiles)


# function for getting the global hypothesis matrix from a list of single hypothesis matrices
global_mat <- function(c_mat, k) matrix(unlist(c_mat), ncol = k, byrow=TRUE)

# wrapper for the asymptotic approach for contrast vectors
wrap_asymptotics <- function(test_stat_list, c_mat_list, k, crit.value.method, alpha = 0.05){
  # calculate Sigma_hat
  Sigma_hat <- diag(test_stat_list$var)
  # if all matrices are partitionized in vectors
  if(all(sapply(c_mat_list, function(hyp_mat) all(sapply(hyp_mat, function(x) nrow(x)) == 1)))){
  # calculate the scaling matrices for all contrast matrices
  D_list <- lapply(c_mat_list, function(c_mat){
    diag(lapply(c_mat, function(vec){
      MASS::ginv(sqrt(vec %*% Sigma_hat %*% t(vec))) #using the g-inverse to avoid problems with zero variances
    }))})
  # calculate the quantiles
  quant <- sapply(1:C, function(D_ind){
    D <- D_list[[D_ind]]
    H <- global_mat(c_mat_list[[D_ind]], k)
    # equicoordinate normal quantiles
    sigma <- (D%*%H%*%Sigma_hat%*%t(H)%*%D)
    # since errors occur when any of diag(sigma) = 0, we do not consider these components
    my_index <- (diag(sigma) != 0)
    (qmvnorm(1 - alpha , tail="both.tails", sigma=sigma[my_index, my_index])$quantile)^2
  })
  # save the test decision
  l <- length(crit.value.method)
  out <- lapply(1:C, function(c_ind){
    out_temp <- t(matrix(rep(test_stat_list$teststats[[c_ind]] > quant[c_ind],l), nrow=l, byrow=TRUE))
    colnames(out_temp) <- crit.value.method
    rownames(out_temp) <- names(test_stat_list$teststats[[c_ind]])
    out_temp
  })
  names(out) <- names(c_mat_list)
  return(out)
}else{
  # generate Nres multivariate random vectors
  random_numbers <- mvtnorm::rmvnorm(Nres , sigma=Sigma_hat)

  decision <- sapply(1:C, function(ind){
    hyp_mat <- c_mat_list[[ind]]
    # determine values for approximating the limiting distribution
    random_values <- t(sapply(hyp_mat, function(mat) apply(random_numbers, 1, function(z) t(mat%*%z)%*%MASS::ginv(mat%*%Sigma_hat%*%t(mat))%*%mat%*%z)))
    out <- matrix(numeric(length(crit.value.method)*length(test_stat_list$teststats[[ind]])),ncol=length(crit.value.method))
    colnames(out) <- crit.value.method
    
    if("inequi" %in% crit.value.method) out[,"inequi"] <- reject_multiple(test_stat_list$teststats[[ind]], crit_values2(random_values))
    #reject_multiple(teststats[[ind]], calc_quantiles(ts_pooledBS[[ind]]))
    if("equi" %in% crit.value.method) out[,"equi"] <- reject_multiple(test_stat_list$teststats[[ind]], quantile(apply(random_values, 2, max), probs = 1 - alpha))
    
    out
  })
  names(decision) <- names(c_mat_list)
  
  return(decision)
  
}
}


# wrapper for the asymptotic approach for the global hypothesis
# Input
# global_teststats:   global Wald-type test statistics
# global_c_mat_list:  list of the global hypothesis matrices
# crit.value.method:  method for calculating the critical value
# Output
# test decisions of the global asymptotic tests
wrap_asymptotics_global <- function(global_teststats, global_c_mat_list, crit.value.method, alpha = 0.05){
  
  C <- length(global_c_mat_list)
  # calculate the quantiles
  quant <- sapply(global_c_mat_list, function(H) qchisq(1 - alpha,df=qr(H)$rank))
  # test decision
  l <- length(crit.value.method)
  out <- lapply(1:C, function(c_ind) t(matrix(rep(any(global_teststats[c_ind] > quant[c_ind]), l),nrow=l,byrow=TRUE)))
  names(out) <- names(global_c_mat_list)
  out
}

# wrapper for the asymptotic approach for multiple hypotheses with bonferroni correction
# Input
# teststats:   local Wald-type test statistics
# c_mat_list:  list of the hypothesis matrices
# crit.value.method:  method for calculating the critical value
# Output
# test decisions of the multiple asymptotic tests
wrap_asymptotics_bonf <- function(teststats, c_mat_list, crit.value.method, alpha = 0.05){
  
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




