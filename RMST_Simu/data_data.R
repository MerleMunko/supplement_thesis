# generate the data based on the data example


setwd(path)

# load the data
dat <- read.csv("gabriel.csv", header = TRUE)

# remove NAs
dat <- dat[!(is.na(dat$q11_4_2x)),]
dat <- dat[(dat$g1r07 != "M"), ]

### mit 2 Faktorstufen bei Faktor A (farm/non-farm)
my_data <- data.frame(time = as.numeric(dat$q11_4_2x), status = as.numeric(dat$q11_4_1),  trt = as.numeric(dat$g1r07), sex = as.numeric(dat$sex)  )

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

tau <- 15
k <- 4

rsurv <- function(m, fit_obj){
  
  x <- runif(m)
  y <- numeric(m)

  ind <- (x >= min(fit_obj$surv))
  y[ind] <- sapply(x[ind], function(z) min(fit_obj$time[fit_obj$surv < z]))
  y[!ind] <- max(fit_obj$time)
  
  return(y)
}

# create the needed fit_obj for the censoring
C_obj <- list()
for(i in 1:k) C_obj[[i]] <- survfit(Surv(my_data$time[my_data$group==i], 
                                         event = 1-my_data$status[my_data$group==i]) ~ 1)

if( alt_ind[1] == 0 | all(alt_ind == -(1:k))){
  
  # create the needed fit_obj
  fit_obj_all <- survfit(Surv(my_data$time, event = my_data$status) ~ 1)
  
  data_gen <- function(n_vec, n){
    
    group <- rep(1:k, n_vec)
    
    T <- rsurv(m = n, fit_obj_all)
    C <- c()
    for(i in 1:k){
      C <- c(C,rsurv(m = n_vec[i], C_obj[[i]]))
    }
    X <- status <- numeric(n)
    for(i in 1:n){ 
       X[i] <- min(T[i],C[i]) 
       status[i] <- as.numeric(X[i] == T[i])
    }
    
    return( data.frame(time = X, status = status, group = group) )
  }
}else{
  # create the needed fit_obj
  fit_obj <- list()
  for(i in 1:k) fit_obj[[i]] <- survfit(Surv(my_data$time[my_data$group==i], 
                                        event = my_data$status[my_data$group==i]) ~ 1)
  
  data_gen <- function(n_vec, n){
    
    X <- status <- numeric(n)
    group <- rep(1:k, n_vec)
    n_cumsum <- c(0,cumsum(n_vec))
    
    for(i in 1:k){
      T <- rsurv(n_vec[i], fit_obj[[i]])
      C <- rsurv(n_vec[i], C_obj[[i]])
      TC <- cbind(T,C)
      X[(1+n_cumsum[i]):n_cumsum[i+1]] <- apply(TC,1,min)
      status[(1+n_cumsum[i]):n_cumsum[i+1]] <- as.numeric(X[(1+n_cumsum[i]):n_cumsum[i+1]] == T)
      }
    
    return( data.frame(time = X, status = status, group = group) )
  }
  
}





