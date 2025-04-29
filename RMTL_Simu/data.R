
### Generate the data
### The following quantities need to be specified in advance
# alt_ind: The subgroup(s) following the alternative distribution
# cens_dis: The censoring distribution, see below in the code for all options
# n_vec: the vector of sample sizes per subgroup
# Distr: The survival distribution setting, see below for all options
# M: Number of risks
# D_probs: A matrix of k columns with the probabilities for the different risks


# The output of this file is the function data_gen(n_vec, n = sum(n_vec)) generating the survival data

n <- sum(n_vec)

# Function to define the indices for the alternative
if( alt_ind[1] == 0 | all(alt_ind == -(1:k))){
  alt_vec <- function(x, ind = alt_ind){ -(1:sum(x))}
  alt_ind <- -(1:k)
}else{
  alt_vec <- function(x, ind = alt_ind){
    x0 <- c(0,cumsum(x))
    out <- NULL
    for(i in sort(ind)){
      out <- c(out,(x0[i]+1):x0[i+1])
    }
    out
  }
}


#### Censoring
r_list <- list()
if( cens_dis == "weib_eq"){
  for(i in 1:k){
    r_list[[i]] <- function(x) rweibull(x,3,10)
  }
}

if( cens_dis == "weib_eq+high"){
  for(i in 1:k){
    r_list[[i]] <- function(x) rweibull(x,3,5)
  }
}

if( cens_dis == "weib_uneq+high"){
  for(i in 1:ceiling(k/4)){
    r_list[[4*i-3]] <- function(x) rweibull(x,0.5,15)
    r_list[[4*i-2]] <- function(x) rweibull(x,0.5,10)
    r_list[[4*i-1]] <- function(x) rweibull(x,1,8)
    r_list[[4*i]]   <- function(x) rweibull(x,1,10)
  }
}

if( cens_dis == "weib_uneq+low"){
  for(i in 1:ceiling(k/4)){
    r_list[[4*i-3]]   <- function(x) rweibull(x,1,20)
    r_list[[4*i-2]] <- function(x) rweibull(x,3,10)
    r_list[[4*i-1]] <- function(x) rweibull(x,1,15)
    r_list[[4*i]] <- function(x) rweibull(x,3,20)
  }
}

if( cens_dis == "unif_eq"){
  for(i in 1:k){
    r_list[[i]] <- function(x) runif(x,0,25)
  }
}

if( cens_dis == "unif_uneq"){
  for(i in 1:ceiling(k/4)){
    r_list[[4*i-3]] <- function(x) runif(x,0,15)
    r_list[[4*i-2]] <- function(x) runif(x,0,20)
    r_list[[4*i-1]] <- function(x) runif(x,0,25)
    r_list[[4*i]]   <- function(x) runif(x,0,30)
  }
}

if( cens_dis == "exp_eq"){
  for(i in 1:k){
    r_list[[i]] <- function(x) rexp(x,0.05)
  }
}

if( cens_dis == "exp_uneq"){
  for(i in 1:ceiling(k/4)){
    r_list[[4*i-3]] <- function(x) rexp(x,0.03)
    r_list[[4*i-2]] <- function(x) rexp(x,0.04)
    r_list[[4*i-1]] <- function(x) rexp(x,0.05)
    r_list[[4*i]]   <- function(x) rexp(x,0.06)
  }
}


#### Survival times (continuous case)
if(discreet == FALSE){
if( distr == "pwExp"){
  tau <- 10
  lambda0 <- 0.2
  rate1 <- 0.5
  rate2 <- 0.05
  
  # survival function and pseudo random number generation of the piecewise exponential hazard
  pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
  pweR <- function(n, salt, h1, h2) {
    pre_salt <- rexp(n, rate = h1-h2)
    pre_salt <- ifelse(pre_salt <= salt, pre_salt, Inf)
    post_salt <- rexp(n, rate = h2)
    apply(cbind(pre_salt, post_salt), 1, min)
  }
  # calculation of the RMST for a piecewise constant hazard set up with variable
  # discontinuity point and fixed tau
  pwRMST <- function(x, end_tau) {
    integrate(function (t) pweS(t, salt  = x, h1 = rate1, h2 = rate2), 
              lower = 0, upper = end_tau)$value
  }
  
  # calculate RMST for Exp(lambda0)
  expRMST <- integrate(function(x) 1 - pexp(x, lambda0), lower = 0, upper = tau)$value
  
  delta0 <- delta 
  # Determine discontinuity point for a RMST difference of delta
  salt <- optimize(function(s) (delta - expRMST + pwRMST(s, tau))^30, 
                   lower = 0.001, upper = tau, tol = .Machine$double.eps)$minimum
  #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
  #lower = 0.01, upper = tau)$minimum
  
  delta_diff <- delta0
  delta_ratio <- expRMST / (expRMST - delta)
  
  
  data_gen_function <- function(x){ 
    out <- numeric(sum(x))
    out[alt_vec(x)] <- pweR(sum(x[alt_ind]), salt, rate1, rate2)
    out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
    out
  }
}


if( distr == "exp_early"){
  tau <- 10
  lambda0 <- 0.2
  dep_point <- 2 # splitting point of the curves
  
  if( delta == 0){
    delta_diff <- 0
    delta_ratio <- 1
    
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- rexp(sum(x[alt_ind]), lambda0)
      out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
      out
    }
  }else{
    # survival function of the piecewise exponential hazard
    pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
    
    # calculation of the RMST for a piecewise constant hazard set up with variable
    # discontinuity point and fixed tau
    pwRMST <- function(x, end_tau) {
      integrate(function (t) pweS(t, salt  = dep_point, h1 = x, h2 = lambda0), 
                lower = 0, upper = tau)$value
    }
    
    # calculate RMST for Exp(lambda0)
    expRMST <- integrate(function(x) 1 - pexp(x, lambda0), lower = 0, upper = tau)$value
    
    delta0 <- delta 
    # Determine h1 for a RMST difference of delta
    rate2 <- optimize(function(s) (delta - expRMST + pwRMST(s, tau))^30, 
                      lower = 0.001, upper = 200, tol = .Machine$double.eps)$minimum
    #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
    #lower = 0.01, upper = tau)$minimum
    
    delta_diff <- delta0
    delta_ratio <- expRMST / (expRMST - delta)
    
    #Pseudo random generator 
    if(lambda0 < rate2){
      pweR <- function(n){
        pre_salt <- rexp(n, rate = rate2 - lambda0)
        pre_salt <- ifelse(pre_salt <= dep_point, pre_salt, Inf)
        post_salt <- rexp(n, rate = lambda0)
        apply(cbind(pre_salt, post_salt), 1, min)
      }
    }
    if(lambda0 > rate2){
      pweR <- function(n){
        pre_salt <- rexp(n, rate = rate2)
        post_salt <- dep_point + rexp(n, rate = lambda0 - rate2)
        apply(cbind(pre_salt, post_salt), 1, min)
      }
    }
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- pweR(sum(x[alt_ind]))
      out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
      out
    }
  }
  
}

if( distr == "exp_late"){
  tau <- 10
  lambda0 <- 0.2
  dep_point <- 2 # splitting point of the curves
  
  if( delta == 0){
    delta_diff <- 0
    delta_ratio <- 1
    
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- rexp(sum(x[alt_ind]), lambda0)
      out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
      out
    }
  }else{
    # survival function of the piecewise exponential hazard
    pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
    
    # calculation of the RMST for a piecewise constant hazard set up with variable
    # discontinuity point and fixed tau
    pwRMST <- function(x, end_tau) {
      integrate(function (t) pweS(t, salt  = dep_point, h1 = lambda0, h2 = x), 
                lower = 0, upper = tau)$value
    }
    
    # calculate RMST for Exp(1)
    expRMST <- integrate(function(x) 1 - pexp(x, lambda0), lower = 0, upper = tau)$value
    
    delta0 <- delta 
    # Determine h2 for a RMST difference of delta
    rate2 <- optimize(function(s) (delta - expRMST + pwRMST(s, tau))^30, 
                      lower = 0.001, upper = 200, tol = .Machine$double.eps)$minimum
    #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
    #lower = 0.01, upper = tau)$minimum
    
    delta_diff <- delta0
    delta_ratio <- expRMST / (expRMST - delta)
    
    #Pseudo random generator 
    if(lambda0 > rate2){
      pweR <- function(n){
        pre_salt <- rexp(n, rate = lambda0 - rate2)
        pre_salt <- ifelse(pre_salt <= dep_point, pre_salt, Inf)
        post_salt <- rexp(n, rate = rate2)
        apply(cbind(pre_salt, post_salt), 1, min)
      }
    }
    if(lambda0 < rate2){
      pweR <- function(n){
        pre_salt <- rexp(n, rate = lambda0)
        post_salt <- dep_point + rexp(n, rate = rate2 - lambda0)
        apply(cbind(pre_salt, post_salt), 1, min)
      }
    }
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- pweR(sum(x[alt_ind]))
      out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
      out
    }
  }
}


if( distr == "exp_prop"){
  # Exp(1) vs Exp(lambda), where lambda is chosen such that the RMST difference equals delta
  tau <- 10
  lambda0 <- 0.2
  #  expRMST <- function(lambda) { (1 - exp(-tau)) - (1 - exp(-lambda*tau))/lambda  }
  # 
  #  if( delta == 0){ 
  #    lambda <- 1
  #  }else{
  #  lambda <- optimize(function(s) (expRMST(s) - delta)^30, 
  #                   lower = 0.001, upper = 10)$minimum
  #  }
  
  if( delta == 0){ 
    lambda <- lambda0
    delta_diff <- delta
    delta_ratio <- 1
  }else{
    lambda <- optimize(function(lambda)( delta - (1 - exp(-tau*lambda0))/lambda0 + (1 - exp(-lambda*tau))/lambda    )^30, 
                       lower = 0.001, upper = 20, tol = .Machine$double.eps)$minimum
    
    delta_diff <- delta
    delta_ratio <- (1 - exp(-tau*lambda0))/lambda0 / ( (1 - exp(-lambda*tau))/lambda )
  }  
  
  data_gen_function <- function(x){ 
    out <- numeric(sum(x))
    out[alt_vec(x)] <- rexp(sum(x[alt_ind]), rate = lambda)
    out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), rate = lambda0)
    out
  }
}



#if( distr == "Weib_pwExp"){
#  shape <- 2
#  scale <- 7
#  tau <- 10
#  rate1 <- 0.15
#  rate2 <- 0.02
  
#  pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
#  pweR <- function(n, salt, h1, h2) {
#    pre_salt <- rexp(n, rate = h1-h2)
#    pre_salt <- ifelse(pre_salt <= salt, pre_salt, Inf)
#    post_salt <- rexp(n, rate = h2)
#    apply(cbind(pre_salt, post_salt), 1, min)
#  }
#  # calculation of the RMST for a piecewise constant hazard set up with variable
#  # discontinuity point and fixed tau
#  pwRMST <- function(x) {
#    integrate(function (t) pweS(t, salt  = x, h1 = rate1, h2 = rate2), 
#              lower = 0, upper = tau)$value
#  }
#  #tau <- p_tau <- optimize(function(x) (integrate(function(x)  pweS(x, const, rate1, rate2), lower = 0, upper = x)$value -
#  #                                       integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = x)$value)^30,
#  #                         lower = 0, upper = 5)$minimum
#  
#  # calculate RMST for Weibull
#  weibRMST <- integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = tau)$value
#  delta0 <- delta 
#  
#  delta_diff <- delta0
#  delta_ratio <- weibRMST / (weibRMST - delta)
#  
#  # Determine discontinuity point for a RMST difference of delta
#  salt <- optimise(function(s) (delta - weibRMST + pwRMST(s))^30, lower = 0.001, upper = tau)$minimum
#  #salt <- optimise(function(s) ( (1 - delta) * weibRMST - pwRMST(s))^30, lower = 0.1, upper = tau)$minimum
#  
#  
#  data_gen_function <- function(x){ 
#    out <- numeric(sum(x))
#    out[alt_vec(x)] <- pweR(sum(x[alt_ind]), salt, rate1, rate2)
#    out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
#    out
#  }
#}

if( distr == "Weib_scale"){
  shape <- 3
  shape2 <- 1.5
  scale <- 8
  tau <- 10
  
  # weibRMST is equal to 6.95
  
  weibRMST <- integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = tau)$value
  weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(x, shape2, s), lower = 0, upper = tau)$value
  
  
  delta0 <- delta 
  # Determine discontinuity point for a RMST difference of delta
  scale2 <- optimise(function(s) (delta - weibRMST + weibRMST2(s))^30, lower = 0.001, upper = 100, tol = .Machine$double.eps)$minimum
  #scale2 <- optimise(function(s) ( (1 - delta) * weibRMST - weibRMST2(s))^30, lower = 0.1, upper = tau)$minimum
  
  delta_diff <- delta0
  delta_ratio <- weibRMST / (weibRMST - delta)
  
  data_gen_function <- function(x){ 
    out <- numeric(sum(x))
    out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
    out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
    out
  }
}

if( distr == "Weib_prop"){
  shape <- 3
  shape2 <- 3
  scale <- 8
  tau <- 10
  
  # weibRMST is equal to 6.95
  
  weibRMST <- integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = tau)$value
  weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(x, shape2, s), lower = 0, upper = tau)$value
  
  
  delta0 <- delta 
  # Determine discontinuity point for a RMST difference of delta
  if(delta != 0) scale2 <- optimise(function(s) (delta - weibRMST + weibRMST2(s))^30, lower = 0.001, upper = 100, tol = .Machine$double.eps)$minimum
  if(delta == 0) scale2 <- scale
  #scale2 <- optimise(function(s) ( (1 - delta) * weibRMST - weibRMST2(s))^30, lower = 0.1, upper = tau)$minimum
  
  delta_diff <- delta0
  delta_ratio <- weibRMST / (weibRMST - delta)
  
  data_gen_function <- function(x){ 
    out <- numeric(sum(x))
    out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
    out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
    out
  }
}


if( distr == "Weib_shape"){
  shape <- 3
  scale <- 8
  scale2 <- 14
  tau <- 10
  
  # weibRMST is equal to 6.95
  
  weibRMST <- integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = tau)$value
  weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(x, s, scale2), lower = 0, upper = tau)$value
  
  
  delta0 <- delta 
  
  delta_diff <- delta0
  delta_ratio <- weibRMST / (weibRMST - delta)
  # Determine discontinuity point for a RMST difference of delta
  #salt <- optimise(function(s) (delta - weibRMST + pwRMST(s))^30, lower = 0.001, upper = tau)$minimum
  shape2 <- optimise(function(s) ( delta -  weibRMST + weibRMST2(s))^30, lower = 0.1, upper = 100, tol = .Machine$double.eps)$minimum
  
  
  data_gen_function <- function(x){ 
    out <- numeric(sum(x))
    out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
    out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
    out
  }
}

if( distr == "Weib_late"){
  shape <- 3
  scale <- 8
  tau <- 10
  
  # weibRMST is equal to 6.95
  
  weibRMST <- integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = tau)$value
  weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(x, shape* s, scale/s), lower = 0, upper = tau)$value
  
  
  delta0 <- delta 
  
  delta_diff <- delta0
  delta_ratio <- weibRMST / (weibRMST - delta)
  # Determine discontinuity point for a RMST difference of delta
  #salt <- optimise(function(s) (delta - weibRMST + pwRMST(s))^30, lower = 0.001, upper = tau)$minimum
  if(delta != 0) s <- optimise(function(s) ( delta -  weibRMST + weibRMST2(s))^30, lower = 0.1, upper = 100, tol = .Machine$double.eps)$minimum
  if(delta == 0) s <- 1
  shape2 <- shape * s
  scale2 <- scale / s
  
  data_gen_function <- function(x){ 
    out <- numeric(sum(x))
    out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
    out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
    out
  }
}

# if( distr == "Weib_early"){
#   shape <- 0.8
#   scale <- 5
#   tau <- 10
#   
#   # weibRMST is equal to 7.17
#   
#   weibRMST <- integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = tau)$value
#   weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(x, shape - s, scale - 5* s), lower = 0, upper = tau)$value
#   
#   
#   delta0 <- delta 
#   
#   delta_diff <- delta0
#   delta_ratio <- weibRMST / (weibRMST - delta)
#   # Determine discontinuity point for a RMST difference of delta
#   #salt <- optimise(function(s) (delta - weibRMST + pwRMST(s))^30, lower = 0.001, upper = tau)$minimum
#   s <- optimise(function(s) ( delta -  weibRMST + weibRMST2(s))^30, lower = -0.2, upper = 0.3)$minimum
#   shape2 <- shape  - s 
#   scale2 <- scale - 5*s
#   
#   data_gen_function <- function(x){ 
#     out <- numeric(sum(x))
#     out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
#     out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
#     out
#   }
# }


# if( distr == "logn+exp"){
#   tau <- 10
#   mean <- 1.5
#   sd <- 1
#   
#   RMSTlogn <- integrate(function(t) 1-pnorm((log(t) - mean)/sd),0,10)$value
#   
#   
#   # calculate RMST for Exp(1)
#   expRMST <- function(lambda){1/lambda - 1/lambda * exp(-lambda*10) }
#   
#   delta0 <- delta 
#   # Determine discontinuity point for a RMST difference of delta
#   lambda <- optimize(function(s) (delta + expRMST(s) - RMSTlogn)^30, 
#                      lower = 0.001, upper = tau)$minimum
#   #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
#   #lower = 0.01, upper = tau)$minimum
#   
#   delta_diff <- delta0
#   delta_ratio <- RMSTlogn/ (RMSTlogn - delta)
#   
#   data_gen_function <- function(x){ 
#     out <- numeric(sum(x))
#     out[alt_vec(x)] <-  rexp(sum(x[alt_ind]), lambda) 
#     out[-alt_vec(x)] <- exp(rnorm(sum(x[-alt_ind]), mean = mean, sd = sd))
#     out
#   }
# }
# 





if( distr == "logn"){
  tau <- 10
  mean <- 2
  sd <- 0.5
  
  RMSTlogn <- integrate(function(t) 1-pnorm( (log(t) - mean)/sd),0,10)$value
  RMSTlogn2 <- function(m) integrate(function(t) 1-pnorm( (log(t) - m)/sd ),0,10)$value
  
  
  delta0 <- delta 
  # Determine discontinuity point for a RMST difference of delta
  if(delta != 0){
    m2 <- optimize(function(s) (delta + RMSTlogn2(s) - RMSTlogn)^30, 
                 lower = 0.001, upper = tau, tol = .Machine$double.eps)$minimum
  }else{
    m2 <- mean
  }
  #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
  #lower = 0.01, upper = tau)$minimum
  
  delta_diff <- delta0
  delta_ratio <- RMSTlogn/ (RMSTlogn - delta)
  
  data_gen_function <- function(x){ 
    out <- numeric(sum(x))
    out[alt_vec(x)] <-  exp(rnorm(sum(x[alt_ind]), mean = m2, sd = sd))
    out[-alt_vec(x)] <- exp(rnorm(sum(x[-alt_ind]), mean = mean, sd = sd))
    out
  }
}


# my special distributions for different survivals with same RMST
if( distr == "Weib_diff"){ # only wors if alt_ind are the last indices
  shape <- 3
  shape2 <- 1.5
  scale <- 8
  scale3 <- 14
  tau <- 10
  
  # weibRMST is equal to 6.95
  
  weibRMST <- integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = tau)$value
  weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(x, shape2, s), lower = 0, upper = tau)$value
  weibRMST3 <- function(s) integrate(function(x) 1 - pweibull(x, s, scale3), lower = 0, upper = tau)$value
  
  
  delta0 <- delta 
  # Determine discontinuity point for a RMST difference of delta
  scale2 <- optimise(function(s) (0 - weibRMST + weibRMST2(s))^30, lower = 0.001, upper = 100, tol = .Machine$double.eps)$minimum
  #scale2 <- optimise(function(s) ( (1 - delta) * weibRMST - weibRMST2(s))^30, lower = 0.1, upper = tau)$minimum
  shape3 <- optimise(function(s) ( 0 -  weibRMST + weibRMST3(s))^30, lower = 0.1, upper = 100, tol = .Machine$double.eps)$minimum
  if(delta != 0){
    shape4 <- optimise(function(s) ( delta -  weibRMST + weibRMST3(s))^30, lower = 0.1, upper = 100, tol = .Machine$double.eps)$minimum
  }else{ shape4 <- shape3 }
  
  delta_diff <- delta0
  delta_ratio <- weibRMST / (weibRMST - delta)
  
  data_gen_function <- function(x){ 
    out <- numeric(sum(x))
    
    for(i in 1:ceiling((k-length(alt_ind))/3)){
      out[alt_vec(x, ind = 3*i-2)] <- rweibull(sum(x[3*i-2]), shape,  scale)
      out[alt_vec(x, ind = 3*i-1)] <- rweibull(sum(x[3*i-1]), shape2, scale2)
      out[alt_vec(x, ind = 3*i)]   <- rweibull(sum(x[3*i])  , shape3, scale3)
    }
    out[alt_vec(x, ind = alt_ind)] <- rweibull(sum(x[alt_ind]), shape4, scale3)
    out
  }
}



if( distr == "pwExp_diff"){
  tau <- 10
  lambda0 <- 0.2
  rate1 <- 0.5
  rate2 <- 0.05
  rate3 <- 0.3
  rate4 <- 0.1
  rate5 <- 1.5
  rate6 <- 0.01
  
  # survival function and pseudo random number generation of the piecewise exponential hazard
  pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
  pweR <- function(n, salt, h1, h2) {
    pre_salt <- rexp(n, rate = h1-h2)
    pre_salt <- ifelse(pre_salt <= salt, pre_salt, Inf)
    post_salt <- rexp(n, rate = h2)
    apply(cbind(pre_salt, post_salt), 1, min)
  }
  # calculation of the RMST for a piecewise constant hazard set up with variable
  # discontinuity point and fixed tau
  pwRMST <- function(x, end_tau, rate1, rate2) {
    integrate(function (t) pweS(t, salt  = x, h1 = rate1, h2 = rate2), 
              lower = 0, upper = end_tau)$value
  }
  
  # calculate RMST for Exp(lambda0)
  expRMST <- integrate(function(x) 1 - pexp(x, lambda0), lower = 0, upper = tau)$value
  
  delta0 <- delta 
  # Determine discontinuity point for a RMST difference of delta
  salt1 <- optimize(function(s) (delta - expRMST + pwRMST(s, tau, rate1, rate2))^30, 
                    lower = 0.001, upper = tau, tol = .Machine$double.eps)$minimum
  salt2 <- optimize(function(s) (0 - expRMST + pwRMST(s, tau, rate3, rate4))^30, 
                    lower = 0.001, upper = tau, tol = .Machine$double.eps)$minimum
  salt3 <- optimize(function(s) (0 - expRMST + pwRMST(s, tau, rate5, rate6))^30, 
                    lower = 0.001, upper = tau, tol = .Machine$double.eps)$minimum
  #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
  #lower = 0.01, upper = tau)$minimum
  
  delta_diff <- delta0
  delta_ratio <- expRMST / (expRMST - delta)
  
  
  data_gen_function <- function(x){ 
    out <- numeric(sum(x))
    for(i in 1:ceiling((k-length(alt_ind))/3)){
      out[alt_vec(x, ind = 3*i-2)] <- rexp(sum(x[3*i-2]), lambda0)
      out[alt_vec(x, ind = 3*i-1)] <- pweR(sum(x[3*i-1]), salt2, rate3, rate4)
      out[alt_vec(x, ind = 3*i)]   <- pweR(sum(x[3*i])  , salt3, rate5, rate6)
    }
    out[alt_vec(x)] <- pweR(sum(x[alt_ind]), salt1, rate1, rate2)
    out
  }
}

################################################################################
#################################### data.gen - functon #######################
###############################################################################

data_gen <- function(n_vec, n){
  status <- numeric(n)
  group <- rep(1:k, n_vec)
  X <- data_gen_function( n_vec )   # realisations for all groups
  C <- D <- numeric(n)
  C[1:n_vec[1]] <-  r_list[[1]]( n_vec[1])
  D[1:n_vec[1]] <-  sample(M,n_vec[1],replace=TRUE,prob = D_probs[1,])
  for( i in 2:k){
    C[sum(n_vec[1:(i-1)]) + (1 : n_vec[i]) ] <- r_list[[i]]( n_vec[i])
    D[sum(n_vec[1:(i-1)]) + (1 : n_vec[i]) ] <- sample(1:M,n_vec[i],replace=TRUE,prob = D_probs[i,])
  }
  cens_ind <- which(X <= C)
  C[cens_ind] <- X[cens_ind]
  status[cens_ind] <- D[cens_ind]
  return( data.frame(X = C, D = status, group = group) )
}

}

#### Survival times (discreet case)
if(discreet == TRUE){
  
  #### Survival times
  
  if( distr == "pwExp"){
    tau <- 10
    lambda0 <- 0.2
    rate1 <- 0.5
    rate2 <- 0.05
    
    # survival function and pseudo random number generation of the piecewise exponential hazard
    pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
    pweR <- function(n, salt, h1, h2) {
      pre_salt <- rexp(n, rate = h1-h2)
      pre_salt <- ifelse(pre_salt <= salt, pre_salt, Inf)
      post_salt <- rexp(n, rate = h2)
      (apply(cbind(pre_salt, post_salt), 1, min))
    }
    # calculation of the RMST for a piecewise constant hazard set up with variable
    # discontinuity point and fixed tau
    pwRMST <- function(x, end_tau) {
      sum(sapply(seq(0,ceiling(tau-1),by=1), function(t) pweS(t, salt  = x, h1 = rate1, h2 = rate2))) +
        (tau-floor(tau))*pweS(floor(tau), salt  = x, h1 = rate1, h2 = rate2)
    }
    
    # calculate RMST for Exp(lambda0)
    expRMST <- integrate(function(x) 1 - pexp(floor(x), lambda0), lower = 0, upper = tau)$value
    
    delta0 <- delta 
    # Determine discontinuity point for a RMST difference of delta
    salt <- optimize(function(s) (delta - expRMST + pwRMST(s, tau))^30, 
                     lower = 0.001, upper = tau, tol = .Machine$double.eps)$minimum
    #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
    #lower = 0.01, upper = tau)$minimum
    
    delta_diff <- delta0
    delta_ratio <- expRMST / (expRMST - delta)
    
    
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- pweR(sum(x[alt_ind]), salt, rate1, rate2)
      out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
      out
    }
  }
  
  
  if( distr == "exp_early"){
    tau <- 10
    lambda0 <- 0.2
    dep_point <- 2 # splitting point of the curves
    
    if( delta == 0){
      delta_diff <- 0
      delta_ratio <- 1
      
      data_gen_function <- function(x){ 
        out <- numeric(sum(x))
        out[alt_vec(x)] <- rexp(sum(x[alt_ind]), lambda0)
        out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
        out
      }
    }else{
      # survival function of the piecewise exponential hazard
      pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
      
      # calculation of the RMST for a piecewise constant hazard set up with variable
      # discontinuity point and fixed tau
      pwRMST <- function(x, end_tau) {
        sum(sapply(seq(0,ceiling(tau-1),by=1), function(t) pweS(t, salt  = dep_point, h1 = x, h2 = lambda0))) +
          (tau-floor(tau))*pweS(floor(tau), salt  = dep_point, h1 = x, h2 = lambda0)
      }
      
      # calculate RMST for Exp(lambda0)
      expRMST <- integrate(function(x) 1 - pexp(floor(x), lambda0), lower = 0, upper = tau)$value
      
      delta0 <- delta 
      # Determine h1 for a RMST difference of delta
      rate2 <- optimize(function(s) (delta - expRMST + pwRMST(s, tau))^30, 
                        lower = 0.001, upper = 1, tol = .Machine$double.eps)$minimum
      #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
      #lower = 0.01, upper = tau)$minimum
      
      delta_diff <- delta0
      delta_ratio <- expRMST / (expRMST - delta)
      
      #Pseudo random generator 
      if(lambda0 < rate2){
        pweR <- function(n){
          pre_salt <- rexp(n, rate = rate2 - lambda0)
          pre_salt <- ifelse(pre_salt <= dep_point, pre_salt, Inf)
          post_salt <- rexp(n, rate = lambda0)
          apply(cbind(pre_salt, post_salt), 1, min)
        }
      }
      if(lambda0 > rate2){
        pweR <- function(n){
          pre_salt <- rexp(n, rate = rate2)
          post_salt <- dep_point + rexp(n, rate = lambda0 - rate2)
          apply(cbind(pre_salt, post_salt), 1, min)
        }
      }
      data_gen_function <- function(x){ 
        out <- numeric(sum(x))
        out[alt_vec(x)] <- pweR(sum(x[alt_ind]))
        out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
        out
      }
    }
    
  }
  
  if( distr == "exp_late"){
    tau <- 10
    lambda0 <- 0.2
    dep_point <- 2 # splitting point of the curves
    
    if( delta == 0){
      delta_diff <- 0
      delta_ratio <- 1
      
      data_gen_function <- function(x){ 
        out <- numeric(sum(x))
        out[alt_vec(x)] <- rexp(sum(x[alt_ind]), lambda0)
        out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
        out
      }
    }else{
      # survival function of the piecewise exponential hazard
      pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
      
      # calculation of the RMST for a piecewise constant hazard set up with variable
      # discontinuity point and fixed tau
      pwRMST <- function(x, end_tau) {
        sum(sapply(seq(0,ceiling(tau-1),by=1), function(t) pweS(t, salt  = dep_point, h1 = lambda0, h2 = x))) +
          (tau-floor(tau))*pweS(floor(tau), salt  = dep_point, h1 = lambda0, h2 = x)
      }
      
      # calculate RMST for Exp(1)
      expRMST <- integrate(function(x) 1 - pexp(floor(x), lambda0), lower = 0, upper = tau)$value
      
      delta0 <- delta 
      # Determine h2 for a RMST difference of delta
      rate2 <- optimize(function(s) (delta - expRMST + pwRMST(s, tau))^30, 
                        lower = 0.001, upper = 1, tol = .Machine$double.eps)$minimum
      #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
      #lower = 0.01, upper = tau)$minimum
      
      delta_diff <- delta0
      delta_ratio <- expRMST / (expRMST - delta)
      
      #Pseudo random generator 
      if(lambda0 > rate2){
        pweR <- function(n){
          pre_salt <- rexp(n, rate = lambda0 - rate2)
          pre_salt <- ifelse(pre_salt <= dep_point, pre_salt, Inf)
          post_salt <- rexp(n, rate = rate2)
          apply(cbind(pre_salt, post_salt), 1, min)
        }
      }
      if(lambda0 < rate2){
        pweR <- function(n){
          pre_salt <- rexp(n, rate = lambda0)
          post_salt <- dep_point + rexp(n, rate = rate2 - lambda0)
          apply(cbind(pre_salt, post_salt), 1, min)
        }
      }
      data_gen_function <- function(x){ 
        out <- numeric(sum(x))
        out[alt_vec(x)] <- pweR(sum(x[alt_ind]))
        out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), lambda0)
        out
      }
    }
  }
  
  
  if( distr == "exp_prop"){
    # Exp(1) vs Exp(lambda), where lambda is chosen such that the RMST difference equals delta
    tau <- 10
    lambda0 <- 0.2
    #  expRMST <- function(lambda) { (1 - exp(-tau)) - (1 - exp(-lambda*tau))/lambda  }
    # 
    #  if( delta == 0){ 
    #    lambda <- 1
    #  }else{
    #  lambda <- optimize(function(s) (expRMST(s) - delta)^30, 
    #                   lower = 0.001, upper = 10)$minimum
    #  }
    
    if( delta == 0){ 
      lambda <- lambda0
      delta_diff <- delta
      delta_ratio <- 1
    }else{
      # survival function of the piecewise exponential hazard
      pweS <- function(t, h1) exp(-t*h1)
      
      # calculation of the RMST for a piecewise constant hazard set up with variable
      # lambda and fixed tau
      pwRMST <- function(x, end_tau) {
        sum(sapply(seq(0,ceiling(tau-1),by=1), function(t) pweS(t, h1 = x))) +
          (tau-floor(tau))*pweS(floor(tau), h1 = x)
      }
      
      # calculate RMST for Exp(1)
      expRMST <- integrate(function(x) 1 - pexp(floor(x), lambda0), lower = 0, upper = tau)$value
      
      lambda <- optimize(function(s) (delta - expRMST + pwRMST(s, tau))^30, 
                         lower = 0.001, upper = 1, tol = .Machine$double.eps)$minimum
      
      delta_diff <- delta
      delta_ratio <- (1 - exp(-tau*lambda0))/lambda0 / ( (1 - exp(-lambda*tau))/lambda )
    }  
    
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- rexp(sum(x[alt_ind]), rate = lambda)
      out[-alt_vec(x)] <- rexp(sum(x[-alt_ind]), rate = lambda0)
      out
    }
  }
  
  
  # 
  # if( distr == "Weib_pwExp"){
  #   shape <- 2
  #   scale <- 7
  #   tau <- 10
  #   rate1 <- 0.15
  #   rate2 <- 0.02
  #   
  #   pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
  #   pweR <- function(n, salt, h1, h2) {
  #     pre_salt <- rexp(n, rate = h1-h2)
  #     pre_salt <- ifelse(pre_salt <= salt, pre_salt, Inf)
  #     post_salt <- rexp(n, rate = h2)
  #     apply(cbind(pre_salt, post_salt), 1, min)
  #   }
  #   # calculation of the RMST for a piecewise constant hazard set up with variable
  #   # discontinuity point and fixed tau
  #   pwRMST <- function(x) {
  #     integrate(function (t) pweS(t, salt  = x, h1 = rate1, h2 = rate2), 
  #               lower = 0, upper = tau)$value
  #   }
  #   #tau <- p_tau <- optimize(function(x) (integrate(function(x)  pweS(x, const, rate1, rate2), lower = 0, upper = x)$value -
  #   #                                       integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = x)$value)^30,
  #   #                         lower = 0, upper = 5)$minimum
  #   
  #   # calculate RMST for Weibull
  #   weibRMST <- integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = tau)$value
  #   delta0 <- delta 
  #   
  #   delta_diff <- delta0
  #   delta_ratio <- weibRMST / (weibRMST - delta)
  #   
  #   # Determine discontinuity point for a RMST difference of delta
  #   salt <- optimise(function(s) (delta - weibRMST + pwRMST(s))^30, lower = 0.001, upper = tau)$minimum
  #   #salt <- optimise(function(s) ( (1 - delta) * weibRMST - pwRMST(s))^30, lower = 0.1, upper = tau)$minimum
  #   
  #   
  #   data_gen_function <- function(x){ 
  #     out <- numeric(sum(x))
  #     out[alt_vec(x)] <- pweR(sum(x[alt_ind]), salt, rate1, rate2)
  #     out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
  #     out
  #   }
  # }
  
  if( distr == "Weib_scale"){
    shape <- 3
    shape2 <- 1.5
    scale <- 8
    tau <- 10
    
    # weibRMST is equal to 6.95
    
    weibRMST <- integrate(function(x) 1 - pweibull(floor(x), shape, scale), lower = 0, upper = tau)$value
    weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(floor(x), shape2, s), lower = 0, upper = tau)$value
    
    
    delta0 <- delta 
    # Determine discontinuity point for a RMST difference of delta
    scale2 <- optimise(function(s) (delta - weibRMST + weibRMST2(s))^30, lower = 0.001, upper = 100, tol = .Machine$double.eps)$minimum
    #scale2 <- optimise(function(s) ( (1 - delta) * weibRMST - weibRMST2(s))^30, lower = 0.1, upper = tau)$minimum
    
    delta_diff <- delta0
    delta_ratio <- weibRMST / (weibRMST - delta)
    
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
      out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
      out
    }
  }
  
  if( distr == "Weib_prop"){
    shape <- 3
    shape2 <- 3
    scale <- 8
    tau <- 10
    
    # weibRMST is equal to 6.95
    
    weibRMST <- integrate(function(x) 1 - pweibull(floor(x), shape, scale), lower = 0, upper = tau)$value
    weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(floor(x), shape2, s), lower = 0, upper = tau)$value
    
    
    delta0 <- delta 
    # Determine discontinuity point for a RMST difference of delta
    if(delta != 0) scale2 <- optimise(function(s) (delta - weibRMST + weibRMST2(s))^30, lower = 0.001, upper = 100, tol = .Machine$double.eps)$minimum
    if(delta == 0) scale2 <- scale
    #scale2 <- optimise(function(s) ( (1 - delta) * weibRMST - weibRMST2(s))^30, lower = 0.1, upper = tau)$minimum
    
    delta_diff <- delta0
    delta_ratio <- weibRMST / (weibRMST - delta)
    
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
      out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
      out
    }
  }
  
  
  if( distr == "Weib_shape"){
    shape <- 3
    scale <- 8
    scale2 <- 14
    tau <- 10
    
    # weibRMST is equal to 6.95
    
    weibRMST <- integrate(function(x) 1 - pweibull(floor(x), shape, scale), lower = 0, upper = tau)$value
    weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(floor(x), s, scale2), lower = 0, upper = tau)$value
    
    
    delta0 <- delta 
    
    delta_diff <- delta0
    delta_ratio <- weibRMST / (weibRMST - delta)
    # Determine discontinuity point for a RMST difference of delta
    #salt <- optimise(function(s) (delta - weibRMST + pwRMST(s))^30, lower = 0.001, upper = tau)$minimum
    shape2 <- optimise(function(s) ( delta -  weibRMST + weibRMST2(s))^30, lower = 0.1, upper = 100, tol = .Machine$double.eps)$minimum
    
    
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
      out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
      out
    }
  }
  
  if( distr == "Weib_late"){
    shape <- 3
    scale <- 8
    tau <- 10
    
    # weibRMST is equal to 6.95
    
    weibRMST <- integrate(function(x) 1 - pweibull(floor(x), shape, scale), lower = 0, upper = tau)$value
    weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(floor(x), shape* s, scale/s), lower = 0, upper = tau)$value
    
    
    delta0 <- delta 
    
    delta_diff <- delta0
    delta_ratio <- weibRMST / (weibRMST - delta)
    # Determine discontinuity point for a RMST difference of delta
    #salt <- optimise(function(s) (delta - weibRMST + pwRMST(s))^30, lower = 0.001, upper = tau)$minimum
    if(delta != 0) s <- optimise(function(s) ( delta -  weibRMST + weibRMST2(s))^30, lower = 0, upper = 2, tol = .Machine$double.eps)$minimum
    if(delta == 0) s <- 1
    
    shape2 <- shape * s
    scale2 <- scale / s
    
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
      out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
      out
    }
  }
  
  # if( distr == "Weib_early"){
  #   shape <- 0.8
  #   scale <- 5
  #   tau <- 10
  #   
  #   # weibRMST is equal to 7.17
  #   
  #   weibRMST <- integrate(function(x) 1 - pweibull(floor(x), shape, scale), lower = 0, upper = tau)$value
  #   weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(floor(x), shape - s, scale - 5* s), lower = 0, upper = tau)$value
  #   
  #   
  #   delta0 <- delta 
  #   
  #   delta_diff <- delta0
  #   delta_ratio <- weibRMST / (weibRMST - delta)
  #   # Determine discontinuity point for a RMST difference of delta
  #   #salt <- optimise(function(s) (delta - weibRMST + pwRMST(s))^30, lower = 0.001, upper = tau)$minimum
  #   s <- optimise(function(s) ( delta -  weibRMST + weibRMST2(s))^30, lower = -0.2, upper = 0.3)$minimum
  #   shape2 <- shape  - s 
  #   scale2 <- scale - 5*s
  #   
  #   data_gen_function <- function(x){ 
  #     out <- numeric(sum(x))
  #     out[alt_vec(x)] <- rweibull(sum(x[alt_ind]), shape2, scale2)
  #     out[-alt_vec(x)] <- rweibull(sum(x[-alt_ind]), shape, scale)
  #     out
  #   }
  # }
  
  
  # if( distr == "logn+exp"){
  #   tau <- 10
  #   mean <- 1.5
  #   sd <- 1
  #   
  #   RMSTlogn <- integrate(function(t) 1-pnorm((log(t) - mean)/sd),0,10)$value
  #   
  #   
  #   # calculate RMST for Exp(1)
  #   expRMST <- function(lambda){1/lambda - 1/lambda * exp(-lambda*10) }
  #   
  #   delta0 <- delta 
  #   # Determine discontinuity point for a RMST difference of delta
  #   lambda <- optimize(function(s) (delta + expRMST(s) - RMSTlogn)^30, 
  #                      lower = 0.001, upper = tau)$minimum
  #   #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
  #   #lower = 0.01, upper = tau)$minimum
  #   
  #   delta_diff <- delta0
  #   delta_ratio <- RMSTlogn/ (RMSTlogn - delta)
  #   
  #   data_gen_function <- function(x){ 
  #     out <- numeric(sum(x))
  #     out[alt_vec(x)] <-  rexp(sum(x[alt_ind]), lambda) 
  #     out[-alt_vec(x)] <- exp(rnorm(sum(x[-alt_ind]), mean = mean, sd = sd))
  #     out
  #   }
  # }
  # 
  
  
  
  
  
  if( distr == "logn"){
    tau <- 10
    mean <- 2
    sd <- 0.5
    
    RMSTlogn <- integrate(function(t) 1-pnorm( (log(floor(t)) - mean)/sd),0,10)$value
    RMSTlogn2 <- function(m) integrate(function(t) 1-pnorm( (log(floor(t)) - m)/sd ),0,10)$value
    
    
    delta0 <- delta 
    # Determine discontinuity point for a RMST difference of delta
    if(delta != 0) {
      m2 <- optimize(function(s) (delta + RMSTlogn2(s) - RMSTlogn)^30, 
                   lower = 0.001, upper = tau, tol = .Machine$double.eps)$minimum
    }else{ m2 <- mean }
    #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
    #lower = 0.01, upper = tau)$minimum
    
    delta_diff <- delta0
    delta_ratio <- RMSTlogn/ (RMSTlogn - delta)
    
    data_gen_function <- function(x){ 
      out <- numeric(sum(x))
      out[alt_vec(x)] <-  exp(rnorm(sum(x[alt_ind]), mean = m2, sd = sd))
      out[-alt_vec(x)] <- exp(rnorm(sum(x[-alt_ind]), mean = mean, sd = sd))
      out
    }
  }
  
  
  # # my special distributions for different survivals with same RMST
  # if( distr == "Weib_diff"){ # only works if alt_ind are the last indices
  #   shape <- 3
  #   shape2 <- 1.5
  #   scale <- 8
  #   scale3 <- 14
  #   tau <- 10
  #   
  #   # weibRMST is equal to 6.95
  #   
  #   weibRMST <- integrate(function(x) 1 - pweibull(x, shape, scale), lower = 0, upper = tau)$value
  #   weibRMST2 <- function(s) integrate(function(x) 1 - pweibull(x, shape2, s), lower = 0, upper = tau)$value
  #   weibRMST3 <- function(s) integrate(function(x) 1 - pweibull(x, s, scale3), lower = 0, upper = tau)$value
  #   
  #   
  #   delta0 <- delta 
  #   # Determine discontinuity point for a RMST difference of delta
  #   scale2 <- optimise(function(s) (0 - weibRMST + weibRMST2(s))^30, lower = 0.001, upper = 100)$minimum
  #   #scale2 <- optimise(function(s) ( (1 - delta) * weibRMST - weibRMST2(s))^30, lower = 0.1, upper = tau)$minimum
  #   shape3 <- optimise(function(s) ( 0 -  weibRMST + weibRMST3(s))^30, lower = 0.1, upper = 100)$minimum
  #   if(delta != 0){
  #     shape4 <- optimise(function(s) ( delta -  weibRMST + weibRMST3(s))^30, lower = 0.1, upper = 100)$minimum
  #   }else{ shape4 <- shape3 }
  #   
  #   delta_diff <- delta0
  #   delta_ratio <- weibRMST / (weibRMST - delta)
  #   
  #   data_gen_function <- function(x){ 
  #     out <- numeric(sum(x))
  #     
  #     for(i in 1:ceiling((k-length(alt_ind))/3)){
  #       out[alt_vec(x, ind = 3*i-2)] <- rweibull(sum(x[3*i-2]), shape,  scale)
  #       out[alt_vec(x, ind = 3*i-1)] <- rweibull(sum(x[3*i-1]), shape2, scale2)
  #       out[alt_vec(x, ind = 3*i)]   <- rweibull(sum(x[3*i])  , shape3, scale3)
  #     }
  #     out[alt_vec(x, ind = alt_ind)] <- rweibull(sum(x[alt_ind]), shape4, scale3)
  #     out
  #   }
  # }
  # 
  # 
  # 
  # if( distr == "pwExp_diff"){
  #   tau <- 10
  #   lambda0 <- 0.2
  #   rate1 <- 0.5
  #   rate2 <- 0.05
  #   rate3 <- 0.3
  #   rate4 <- 0.1
  #   rate5 <- 1.5
  #   rate6 <- 0.01
  #   
  #   # survival function and pseudo random number generation of the piecewise exponential hazard
  #   pweS <- function(t, salt, h1, h2) exp(-ifelse(t <= salt, h1*t, h1 * salt + h2 * (t - salt)))
  #   pweR <- function(n, salt, h1, h2) {
  #     pre_salt <- rexp(n, rate = h1-h2)
  #     pre_salt <- ifelse(pre_salt <= salt, pre_salt, Inf)
  #     post_salt <- rexp(n, rate = h2)
  #     apply(cbind(pre_salt, post_salt), 1, min)
  #   }
  #   # calculation of the RMST for a piecewise constant hazard set up with variable
  #   # discontinuity point and fixed tau
  #   pwRMST <- function(x, end_tau, rate1, rate2) {
  #     integrate(function (t) pweS(t, salt  = x, h1 = rate1, h2 = rate2), 
  #               lower = 0, upper = end_tau)$value
  #   }
  #   
  #   # calculate RMST for Exp(lambda0)
  #   expRMST <- integrate(function(x) 1 - pexp(x, lambda0), lower = 0, upper = tau)$value
  #   
  #   delta0 <- delta 
  #   # Determine discontinuity point for a RMST difference of delta
  #   salt1 <- optimize(function(s) (delta - expRMST + pwRMST(s, tau, rate1, rate2))^30, 
  #                     lower = 0.001, upper = tau)$minimum
  #   salt2 <- optimize(function(s) (0 - expRMST + pwRMST(s, tau, rate3, rate4))^30, 
  #                     lower = 0.001, upper = tau)$minimum
  #   salt3 <- optimize(function(s) (0 - expRMST + pwRMST(s, tau, rate5, rate6))^30, 
  #                     lower = 0.001, upper = tau)$minimum
  #   #salt <- optimize(function(s) ( (1 - delta) * expRMST - pwRMST(s, tau))^30, 
  #   #lower = 0.01, upper = tau)$minimum
  #   
  #   delta_diff <- delta0
  #   delta_ratio <- expRMST / (expRMST - delta)
  #   
  #   
  #   data_gen_function <- function(x){ 
  #     out <- numeric(sum(x))
  #     for(i in 1:ceiling((k-length(alt_ind))/3)){
  #       out[alt_vec(x, ind = 3*i-2)] <- rexp(sum(x[3*i-2]), lambda0)
  #       out[alt_vec(x, ind = 3*i-1)] <- pweR(sum(x[3*i-1]), salt2, rate3, rate4)
  #       out[alt_vec(x, ind = 3*i)]   <- pweR(sum(x[3*i])  , salt3, rate5, rate6)
  #     }
  #     out[alt_vec(x)] <- pweR(sum(x[alt_ind]), salt1, rate1, rate2)
  #     out
  #   }
  # }
  
  ################################################################################
  #################################### data.gen - functon #######################
  ###############################################################################
  
  data_gen <- function(n_vec, n){
    status <- numeric(n)
    group <- rep(1:k, n_vec)
    # round the data
    X <- ceiling(data_gen_function( n_vec ))   # realisations for all groups
    C <- D <- numeric(n)
    C[1:n_vec[1]] <-  r_list[[1]]( n_vec[1])
    D[1:n_vec[1]] <-  sample(M,n_vec[1],replace=TRUE,prob = D_probs[1,])
    for( i in 2:k){
      C[sum(n_vec[1:(i-1)]) + (1 : n_vec[i]) ] <- r_list[[i]]( n_vec[i])
      D[sum(n_vec[1:(i-1)]) + (1 : n_vec[i]) ] <- sample(1:M,n_vec[i],replace=TRUE,prob = D_probs[i,])
    }
    cens_ind <- which(X <= C)
    C[cens_ind] <- X[cens_ind]
    status[cens_ind] <- D[cens_ind]
    return( data.frame(X = C, D = status, group = group) )
  }
  
  
}