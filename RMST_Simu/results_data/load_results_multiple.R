
### View the results for the local hypotheses

# clear the workspace
#rm(list=ls())

# Set the working directory
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# get all result file names
names_error <- list.files(paste0(path))
names2 <- gsub('.{6}$', "", names_error)

# load the data of Setting 1
load(paste0(path,names_error[1]))

# create a result array for the errorrates
error <- vector(mode = "list", length = length(results_multiple))
for(j in 1:length(results_multiple)){
  error[[j]] <- array(NA, dim = c(length(names2),dim(results_multiple[[j]])),
                      dimnames = c(list(names2),dimnames(results_multiple[[j]])))
}

for(i in 1:length(names_error)){for(j in 1:length(results_multiple)){
  # load the results
  load(paste0(path,names_error[i]))
  # save the multiple results
  error[[j]][i,,,] <- results_multiple[[j]]
}}
names(error) <- names(results_multiple)
dimnames(error[[2]])[[3]] <- 1:2

# only consider the methods of interest 
# do not consider global approaches!
#methods <- c("asymptotic","wild, Rademacher","wild, Gaussian","groupwise")
methods <- c("asymptotic","groupwise")
coordinates <- c("inequi")
#merken <- dim(error[[1]])
error_new <- lapply(error, function(x) x[,methods,,coordinates])
#dim(error_new) <- c(merken[1],length(methods),merken[3],length(coordinates))
names(error_new) <- names(error)
error <- error_new


#######################
###      Tables     ###
#######################

# view the error rates
#for(j in 1:length(results_multiple)){
#  View(as.data.frame(error[[j]]))
#}
gesplittet <- (strsplit(names2, split = "_"))
delta <- sapply(gesplittet, function(x) x[2])
n <- sapply(gesplittet, function(x) paste(x[3], collapse=" "))
distribution <- sapply(gesplittet, function(x) paste(x[5:(length(x)-2)], collapse=" "))
censoring_distribution <- sapply(gesplittet, function(x) paste(x[(length(x)-1):length(x)], collapse = " "))
for(d in unique(delta)[-1]){for(j in 1:length(error)){
  #identify the alternative hypotheses
  if(j == 1) hypothesis <- "FaktorB"
  if(j == 2) hypothesis <- c(1,2)
  for(n_temp in unique(n)){
      #for(hyp in hypothesis){
      error_table <- error[[j]][n== n_temp & delta==d,,hypothesis]
      error_table[!is.na(as.numeric(error_table))] <- sprintf("%.3f",as.numeric(error_table[!is.na(as.numeric(error_table))]))
      error_table <- data.frame( distribution[n== n_temp & delta==d], censoring_distribution[n== n_temp & delta==d], error_table)
      colnames(error_table) <- c( "distribution", "censoring distribution", rep(methods, length(hypothesis)))
      print(xtable(error_table,
                   digits = 3,
                   caption = paste("Rejection rates for the",
                                   names(error)[j],
                                   "contrast matrix with $\\delta =$",
                                   d,
                                   "and with sample size",
                                   n_temp)),
            include.rownames = FALSE, sanitize.text.function = function(x) {x})
    }}
  }
  #}



#######################
###     H1 plots    ###
#######################

err_alt <- lapply(error, function(x) x[as.numeric(sapply(strsplit(names2, split ="_"), function(x) x[2])) != 0,,,drop=FALSE])
names(err_alt) <- names(error)
par(mfrow=c(1,1),mar = c(8.1,4.1,2.1,2.1))
for(j in 1:length(err_alt)){
  my_dataframe <- as.data.frame(err_alt[[j]][,,])
  #ord_vec <- c(matrix(c(1:(ncol(my_dataframe)/2),(ncol(my_dataframe)/2+1):(ncol(my_dataframe))),ncol=ncol(my_dataframe)/2,byrow=T))
  boxplot(my_dataframe,
          las = 2,
          col = "palegreen1",
          ylim = c(0,max(err_alt[[j]])+0.03), ylab ="rejection rate")
  #abline(h = 0.05, col="red")
  #legend("topleft", c("equicoordinate", "inequicoordinate"), col = c("lightskyblue1", "palegreen1"), pch = 16, bty = "n")
}


#######################
###   Power plots   ###
#######################


names4 <- substr(names2, 14, nchar(names2))
x <- unique(as.numeric(substr(names2, 10,12)))

palette <- c("red", "orange", "black", "green", "blue", "purple")
for(nam in unique(names4)){
  ind <- which(names4 == nam)
  
  
  for(k in 1:length(error)){
    f_set <- which(as.character(dimnames(error[[k]])[[3]]) %in% c("FaktorA", 1, 2))
    par(mfrow = c(1,length(f_set)))
    for(f in f_set){
    m <- 1
    if(length(error[[k]][ind,1,f])==length(x)){
      plot(x,error[[k]][ind,1,f],type="l",col=palette[1], ylim=c(0,1),
           ylab = "rejection rate", xlab = expression(delta), main = dimnames(error[[k]])[[3]][f])
      for(l in 2:dim(error[[k]])[2]){
        m <- m+1
        lines(x,error[[k]][ind,l,f],col=palette[m])
      }
      legend("topleft",col=palette, dimnames(error[[k]])[[2]], pch=16, bty = "n")}
  }}}
# achtung: nur fuer inequi coordinaten
# auch lieber eine Tabelle hier...

##########################
###   Power barplots   ###
##########################

#rm(list=ls())
# only consider the methods of interest 
# do not consider global approaches!
methods <- c("asymptotic","wild, Rademacher","wild, Gaussian","groupwise")
coordinates <- c("inequi")
# treffe eine auswahl an plots
#auswahl <- 1:90
auswahl <- 1:5


# load all decision_arr objects
# Set the working directory
path <- "C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results_data/decision_list/"
setwd(path)
# get all result file names
names_dec <- list.files(path)
x <- unique(as.numeric(substr(names_dec, 10,12)))
x_char <- sprintf("%.1f", x)
# save the objects in a list
load(names_dec[1])
# save the data in an array
dec <- vector(mode = "list", length = nrow(decision_arr))
names(dec) <- dimnames(decision_arr)[[1]]
for(j in 1:length(dec)){
  v1 <- (dimnames(decision_arr[j,1][[1]])[[2]])
  v_list <- sapply(do.call("c", lapply(seq_along(v1), function(i) combn(v1, i, FUN = list))), paste, collapse = " and ")
  dec[[j]] <- array(NA, dim = c(length(names_dec),length(methods),2^(dim(decision_arr[j,1][[1]])[2])-1,length(coordinates)),
                      dimnames = list(names_dec,methods,v_list,coordinates))
}
for(j in 1:length(dec)){
  # define all possible combinations of rejected hypotheses
  v1 <- 1:(dim(decision_arr[j,1][[1]])[2])
  v_list <- do.call("c", lapply(seq_along(v1), function(i) combn(v1, i, FUN = list)))
  
  for(i in 1:(length(auswahl))){ for(xi in x_char){
  # load the data
  name_auswahl <- gsub("0.0", xi, names_dec[auswahl[i]])
  load(name_auswahl)
    for(meth in methods){ for(coor in coordinates){
  decisions <- simplify2array(decision_arr[j,])[meth,,coor,]
  results_dec <- colMeans(sapply(1:length(v_list), function(v){
    apply(decisions, 2, function(x) all.equal(which(as.logical(x)), v_list[[v]])==TRUE)
    }))
  dec[[j]][which(names_dec == name_auswahl),meth,,coor] <- results_dec
}}}}}

names4 <- substr(names_dec, 14, nchar(names_dec)-6)


for(nam in auswahl){
  ind <- which(names4 == substr(names_dec[nam], 14, nchar(names_dec[nam])-6))
  for(k in 1:length(dec)){
    
    par(mfrow = c(2,ceiling((length(methods)+1)/2)))
    palette <- rainbow(dim(dec[[k]])[3])
    for(meth in methods){
    
      box_heights <- t(apply(dec[[k]][ind[order(x)],meth,,1], 1, cumsum))
      barplot(box_heights[,ncol(box_heights)], col = palette[ncol(box_heights)], 
              border = FALSE, names.arg = sort(x), xlab = expression(delta), ylim = c(0,1), main = meth)
      for(column in (ncol(box_heights)-1):1){
      barplot(box_heights[,column], add=TRUE, col = palette[column], border = FALSE, names.arg = sort(x))
      }
    }
    
#    if(dim(dec[[k]])[3]<= 15){
      plot.new()
      text(0.5,1, labels = "Which hypotheses\n are rejected? ")
      if(is.null(dimnames(dec[[k]])[[3]])){
        text <- c(1:2, "1 and 2")
      }else{ text <- dimnames(dec[[k]])[[3]]}
      legend(0,0.9,col=palette, text, pch=16, bty = "n")
#    }
    # if(dim(dec[[k]])[3]> 15){
    #   par(mfrow = c(1,1))
    #   plot.new()
    #   text(0.5,1, labels = "Which hypotheses\n are rejected? ")
    #   legend(0,0.9,col=palette, paste(dimnames(dec[[k]])[[3]]), pch=16, bty = "n")
    #   par(mfrow = c(1,length(methods)+1))
    #   }
  }}

# achtung: nur fuer inequi coordinaten 
# lieber tabelle... und ein Beispiel.