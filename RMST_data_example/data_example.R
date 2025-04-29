
library(survival)
library(multcomp)
library(muhaz)

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

source("functions_factorial.R")
source("../RMST_Simu/functions_multiple.R")
source("../RMST_Simu/resampling_multiple.R")
source("p.value_functions.R")

RMSTANOVA <- function(my_data, c_mat_list, tau, crit.value.method = "inequi", Nres = 1999){
  
  # determine the number of groups and number of (global) hypotheses
  k <- length(unique(my_data$group))
  C <- length(c_mat_list)
  
  # print the censoring rates
  print(paste("Censoring rate in group", 1:k, ":", sapply(1:k, function(i) mean(1 - my_data$status[my_data$group == i]))))
  
  # calculate the test statistics, global hypothesis matrices and global test statistics
  test_stat_list <- wrap_multiple_test_stat(my_data, c_mat_list, tau)
  global_c_mat_list <- lapply(c_mat_list, function(c_mat) global_mat(c_mat, k))
  global_teststats  <- lapply(global_c_mat_list, function(mat) test_stat(my_data, tau, mat))
  
  out <- list()
  
  # Global Hypothesis: Asymptotics
  out[["asymptotic_global"]] <- wrap_asymptotics_global(global_teststats, global_c_mat_list, crit.value.method)
  
  # Global Hypothesis: Permutation
  out[["permutation"]] <- as.list(wrap_perm(my_data, tau, global_c_mat_list, Nres, global_teststats, crit.value.method))
  
  ## Multiple Testing
  # Asymptotics wrap_asymptotics(test_stat_list, c_mat_list, k, crit.value.method)
  Sigma_hat <- diag(test_stat_list$var)
  random_numbers <- (rmvnorm(Nres , sigma=Sigma_hat))
  random_values <- lapply(c_mat_list, function(c_mat) (t(sapply(c_mat, function(mat) apply(random_numbers, 1, function(z) t(mat%*%z)%*%MASS::ginv(mat%*%Sigma_hat%*%t(mat))%*%mat%*%z)))))
  out[["asymptotic"]] <- lapply(1:C, function(c_ind) p.value(data_mat =  random_values[[c_ind]], teststat=test_stat_list$teststats[[c_ind]]))
  
  # wild Bootstrap test, Rademacher multipliers
  out[["wild, Rademacher"]] <- wrap_wildBS(my_data, tau, c_mat_list, Nres, test_stat_list, n_vec, multiplier = "Rademacher", crit.value.method)

  # wild Bootstrap test, Gaussian multipliers
  out[["wild, Gaussian"]] <- wrap_wildBS(my_data, tau, c_mat_list, Nres, test_stat_list, n_vec, multiplier = "Gaussian", crit.value.method)

  # groupwise Bootstrap test
  out[["groupwise"]] <- wrap_groupwiseBS(my_data, tau, c_mat_list, Nres, test_stat_list, n_vec, crit.value.method)

  ## Bonferroni approaches
  # Global Hypothesis: Asymptotics
  out[["asymptotic_bonf"]] <- as.list(wrap_asymptotics_bonf(test_stat_list$teststats, c_mat_list, crit.value.method))
  
  # Global Hypothesis: Permutation
  out[["permutation_bonf"]] <- wrap_perm_bonf(my_data, tau, c_mat_list, Nres, test_stat_list$teststats, crit.value.method)
  
  
  for(i in 1:length(out)){
    names(out[[i]]) <- paste("Matrix", 1:length(out[[i]]))
    for(j in 1:length(out[[i]])) if(length(out[[i]][[j]]) == length(c_mat_list[[j]])) names(out[[i]][[j]]) <- names(c_mat_list[[j]])
  }
  
  out
}

#######################################################
##########      Hay fever             #################
#######################################################

# load the data
dat <- read.csv("gabriel.csv", header = TRUE)
## IMPORTANT NOTE: The data set in this repository is NOT the original data set 
## due to ethical restrictions and informed consent limitations.
## Here, we just consider a similar example data set.
## Data access for the original data set can be requested from Jon Genuneit.

# remove NAs
dat <- dat[!(is.na(dat$q11_4_2x)),]
dat <- dat[(dat$g1r07 != "M"), ]


### 2 levels for faktor A (farm/non-farm)
my_data <- data.frame(time = as.numeric(dat$q11_4_2x), 
                      status = as.numeric(dat$q11_4_1),  
                      trt = as.numeric(dat$g1r07), 
                      sex = as.numeric(dat$sex)  )


group <- numeric(length(my_data$status))
gr_count <- 1
group_s <- numeric()
for(i in unique(my_data$trt)){
  for(j in unique(my_data$sex)){
    ind <- (my_data$trt == i) & (my_data$sex == j)
    group_s[gr_count] <- sum(ind)
    group[ind] <- gr_count
    gr_count <- gr_count + 1
  }
}

my_data$group <- as.factor(group)

#pdf("Surv+haz.pdf", width=12.5,height=6.25)
#tiff("Surv+haz.tiff", width=12.5,height=6.25, units = "in", res = 600)
par(mfrow=c(1,2))
surv_obj <- Surv(my_data$time, event = my_data$status)
plot(survfit(surv_obj ~ my_data$group), 
     ylim = c(0.7,1), col = c("grey1", "grey", "grey1", "grey"), 
     lty = c("solid", "solid", "dashed", "dashed"),
     xlab = "Time in years", ylab = "Survival probability")
legend("bottomleft", col = c("grey1", "grey", "grey1", "grey"), 
       lty = c("solid", "solid", "dashed", "dashed"),
       c("boys growing up on a farm", 
         "girls growing up on a farm", 
         "boys not growing up on a farm", 
         "girls not growing up on a farm"), bty = "n")

plot(survfit(surv_obj ~ my_data$group), 
     ylim = c(0,0.35), col = c("grey1", "grey", "grey1", "grey"), 
     lty = c("solid", "solid", "dashed", "dashed"),
     xlab = "Time in years", ylab = "Cumulative hazard",
     cumhaz = TRUE)
legend("topleft", col = c("grey1", "grey", "grey1", "grey"), 
       lty = c("solid", "solid", "dashed", "dashed"),
       c("boys growing up on a farm", 
         "girls growing up on a farm", 
         "boys not growing up on a farm", 
         "girls not growing up on a farm"), bty = "n")
#dev.off()

tau <- 15
for(i in levels(my_data$group)){
  haz <- muhaz(my_data$time, my_data$status, subset = (my_data$group==i), max.time = tau)
  if(i == levels(my_data$group)[1]) plot(haz, ylim = c(0,0.04), col="grey1")
  else lines(haz, col = c("grey1", "grey", "grey1", "grey")[which(i==levels(my_data$group))], 
             lty = c("solid", "solid", "dashed", "dashed")[which(i==levels(my_data$group))])
}
legend("topleft", col = c("grey1", "grey", "grey1", "grey"), 
       lty = c("solid", "solid", "dashed", "dashed"), 
       c("boys growing up on a farm", 
         "girls growing up on a farm", 
         "boys not growing up on a farm", 
         "girls not growing up on a farm"), bty = "n")

# Cox proportional hazard model
my_data$trt <- as.factor(my_data$trt)
my_data$sex <- as.factor(my_data$sex)
(cox_model <- coxph(Surv(time, status) ~ trt * sex, data=my_data))
cox.zph(cox_model)
# prop. hazard assumption is rejected

n_vec <- c(table(group))
k <- length(n_vec)
a <- 2; b <- 2
HA <- null_mat_A(a, b)
HB <- null_mat_B(a, b)
HAB <- null_mat_AB(a, b)

# Define the partitionized hypothesis matrix
# test simultaneously whether there is an effect of factor A, B or an interaction effect
c_mat1 <- list(FaktorA = HA, FaktorB = HB, Interaktion = HAB)
#c_mat2 <- list(t(c(-1,0,1,0)), t(c(0,-1,0,1)))
c_mat_list <- list(c_mat1)
C <- length(c_mat_list)

my_data$group <- as.numeric(my_data$group)
RNGkind("L'Ecuyer-CMRG")
set.seed(1)
(RMST_hayfever <- RMSTANOVA(my_data, c_mat_list, tau = 15, crit.value.method = "inequi",Nres = 19999))
# in contrast to the cox model, we now see a significance for the sex
# and this even simultaneously with the other tests!

# show the table as in the thesis
sapply(RMST_hayfever[-c(1,2)], function(x) if(is.list(x[[1]])){
  return(round(x[[1]][[1]],3))
}else{ 
  return(round(x[[1]],3))
})

# print the RMSTs and the variance estimators
for(i in 1:4){
  daten <- my_data[my_data$group == i,]
  print(RMST(daten$time, daten$status, tau))
}
