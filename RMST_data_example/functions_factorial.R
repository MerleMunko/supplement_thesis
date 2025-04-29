#Null hypothesis matrix

#classical null hypothesis
#Input:
# n_group:   - integer, number of groups
#Output:
# Returns the matrix testing the Nullhypothesis for the respective number of 
# groups
null_mat_x1 <- function(n_group_a, n_group_b){
  n_group  <- n_group_a * n_group_b
  return(diag(n_group) - matrix(1, n_group, n_group)/n_group)
}
#null_mat_clas(5)



#No main effect of factor a
#Input:
# n_group_a:  - integer, number of levels in factor a
# n_group_b:  - integer, number of levels in factor b
#Output:
# Returns the matrix testing the Nullhypothesis 
null_mat_A <- function(n_group_a, n_group_b) {
  pa <- diag(n_group_a) - matrix(1, n_group_a, n_group_a)/n_group_a
  jbb <- matrix(1,1, n_group_b)/n_group_b
  ha <- pa %x% jbb
  
  return(ha)
}

#No main effect of factor B
#Input:
# n_group_a:  - integer, number of levels in factor a
# n_group_b:  - integer, number of levels in factor b
#Output:
# Returns the matrix testing the Nullhypothesis 
null_mat_B <- function(n_group_a, n_group_b) {
  pb <- diag(n_group_b) - matrix(1, n_group_b, n_group_b)/n_group_b
  jaa <- matrix(1, 1, n_group_a)/n_group_a
  ha <- jaa %x% pb
  
  return(ha)
}


#null_mat_fac(2, 2)

#No interaction effect
#Input:
# n_group_a:  - integer, number of levels in factor a
# n_group_b:  - integer, number of levels in factor b
#Output:
# Returns the matrix testing the Nullhypothesis 
null_mat_AB <- function(n_group_a, n_group_b){
  pa <- diag(n_group_a) - matrix(1, n_group_a, n_group_a)/n_group_a
  pb <- diag(n_group_b) - matrix(1, n_group_b, n_group_b)/n_group_b
  hab <- pa %x% pb

  return(hab)
}

#projector
#Input:
# h       - Matrix
#Output: projected matrix

project <- function(h) {
  return(t(h) %*% MASS:::ginv(h %*% t(h)) %*% h)
}



#wrapper for RMST several groups in data
#Input:
# values    - matrix created by dat_gen, sorted by sort_data and with added 
#             columns by KME. Should contain several groups
# group:    - integer vector containing the group of the observations. Default
#             is the third column of the values, the groups drawn by data_gen
#Output:    Standard deviation estimated with interval formula for each group

test_stat <- function(values, group, tau, c_mat) {
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





#wrapper to get the teststatistics directly from the sorted data
#Input
# values:   matrix. The data to start with
# group:    integer vector containing the group of the observations. Default is
#           the third column of the values, the groups drawn by data_gen
wrap_sim <- function(values, group, c_mat_list, hyp_list, tau){
  
  erg_int <- list()
  
  for( hyp in hyp_list){
  c_mat <- c_mat_list[[hyp]]
    erg_int[[hyp]] <- test_stat(values = values, group = group, c_mat = c_mat_list[[hyp]], tau = tau)
  }
  out <- c( unlist(erg_int))
  names(out) <- c( paste0(hyp_list)) 
  return(out)
}


#wrap_sim(values = test_dat,
#         c_mat = null_mat_clas(3))



#Permutations
#Input:
# values    - Matrix, Data to be entered in Simulation
# n_perm    - Integer, Number of permutations
# c_mat     - matrix, test matrix
#Output
# Vector with the quantiles of the test-statistics and number of permutations
# 

perm_fun <- function(values, n_perm, c_mat_list = c_mat_list, hyp_list = hyp_list, tau) {
  values2 <- sort_data(values)
  group_org <- values2[, 3]
  group_new <- replicate(n_perm, sample(group_org))
  test_stat_erg <- matrix(apply(group_new, 2, 
                         function(x) wrap_sim(values = values2, group = x, 
                                           c_mat_list = c_mat_list, hyp_list = hyp_list, tau = tau) ), nrow = length(hyp_list))
  row.names(test_stat_erg) <- hyp_list
  q <- list()
  for( hyp in hyp_list){
    q1 <- unname(quantile(test_stat_erg[hyp, ], 0.95, na.rm = TRUE)) 
    q[[hyp]] <- c(q1)
    names(q[[hyp]]) <- paste0("q")
}
  return(list(q = q, test_stat_erg = test_stat_erg ) )   
}


#wrapper fuer testoutput
sim_test <- function(values, n_perm, c_mat_list, hyp_list, chi_quant, group_s, a, b, tau) {
  #values = sort_data(do.call(fun, args = list(n_vec = group_s,  censp = censp_s)))
  erg_stat <- wrap_sim(values, group = values[,3], c_mat_list = c_mat_list, hyp_list = hyp_list, tau = tau)
  out <- list() 
    erg_perm <- perm_fun(values, n_perm, c_mat_list = c_mat_list, hyp_list = hyp_list, tau = tau)
    for(hyp in hyp_list){ 
      q_perm <- erg_perm$test_stat_erg
      t_perm <- mean(erg_stat[hyp] <= q_perm[hyp, ], na.rm = TRUE)
      t_chi <- 1-pchisq(erg_stat[hyp], df = qr(c_mat_list[[hyp]])$rank )
      
      out1 <- c("perm" = t_perm, "chi" = unname(t_chi))
      out[[hyp]] <- out1
    }
   return(out)
}
