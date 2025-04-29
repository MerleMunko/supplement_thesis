
### View the results for the local hypotheses
library(xtable)
library(abind)

# clear the workspace
#rm(list=ls())

# Set the working directory
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# should the plots being saved (or just displayed?)
save_plots <- FALSE
# as pdf?
save_plots_pdf <- TRUE

# get all result file names
names_error <- list.files(paste0(path))
names2 <- gsub('.{6}$', "", names_error)

# load the data of Setting 1
load(paste0(path,names_error[1]))
anz_nonbonfs <- dim(results_multiple[[1]])[1]

# create a result array for the errorrates
error <- vector(mode = "list", length = length(results_multiple))
for(j in 1:length(results_multiple)){
  error[[j]] <- array(NA, dim = c(length(names2),dim(results_multiple[[j]])),
                      dimnames = list(names2,c(dimnames(results_multiple[[j]])[[1]]),
                                      dimnames(results_multiple[[j]])[[2]], 
                                      dimnames(results_multiple[[j]])[[3]]))
}

for(i in 1:length(names_error)){for(j in 1:length(results_multiple)){
  # load the results
  load(paste0(path,names_error[i]))
  # save the multiple results
  error[[j]][i,1:anz_nonbonfs,,] <- results_multiple[[j]]
  
  
}}
names(error) <- names(results_multiple)

# only consider the methods of interest 
# do not consider global approaches!
#methods <- c("asymptotic","wild, Rademacher","wild, Gaussian","groupwise")
#methods <- c("asymptotic","groupwise")
methods <- c("asymptotic","groupwise", "asymptotic_bonf","permutation_bonf")
coordinates <- c("inequi")
color <- c("darkgrey")
#merken <- dim(error[[1]])
error_new <- lapply(error, function(x) x[,methods,,coordinates,drop = FALSE])
#dim(error_new) <- c(merken[1],length(methods),merken[3],length(coordinates))
names(error_new) <- names(error)
error <- error_new

# combine the equal settings, i.e. exp early, late & prop under the null
# and Weib late & prop under the null
gesplittet <- (strsplit(names2, split = "_"))
delta <- sapply(gesplittet, function(x) x[2])
n <- sapply(gesplittet, function(x) paste(x[3:4], sep=" ", collapse = "_"))
distribution <- sapply(gesplittet, function(x) paste(x[5:(length(x)-2)], sep=" ", collapse = "_"))
censoring_distribution <- sapply(gesplittet, function(x) paste(x[(length(x)-1):length(x)], sep=" ", collapse = "_"))

names2 <- dimnames(error[[1]])[[1]]

#######################
###      Tables     ###
#######################

# view the error rates
#for(j in 1:length(results_multiple)){
#  View(as.data.frame(error[[j]]))
#}
gesplittet <- (strsplit(names2, split = "_"))
delta <- sapply(gesplittet, function(x) x[2])
n <- sapply(gesplittet, function(x) paste(x[3:4], collapse=" "))
distribution <- sapply(gesplittet, function(x) paste(x[5:(length(x)-2)], collapse=" "))
censoring_distribution <- as.factor(sapply(gesplittet, function(x) paste(x[(length(x)-1):length(x)], collapse = " ")))
# rename the censoring distribution
levels(censoring_distribution) <- c("equal", "unequal, high", "unequal, low")
capture.output({
  for(d in unique(delta)[-1]){for(j in 1:length(error)){
    #identify the alternative hypotheses
    if(j == 1) hypothesis <- "4 - 1"
    if(j == 2) hypothesis <- c("4 - 1", "4 - 2", "4 - 3")
    if(j == 3) hypothesis <- dimnames(error[[j]])[[3]] # all hypotheses are wrong
    for(coor in 1:dim(error[[j]])[4]){
      for(n_temp in unique(n)){
        for(hyp in hypothesis){
          error_table <- error[[j]][n== n_temp & delta==d,,hyp,coor]
          error_table[!is.na(as.numeric(error_table))] <- sprintf("%.3f",as.numeric(error_table[!is.na(as.numeric(error_table))]))
          error_table <- data.frame( distribution[n== n_temp & delta==d], censoring_distribution[n== n_temp & delta==d], error_table)
          colnames(error_table) <- c( "distribution", "censoring distribution", gsub(",","", gsub("_"," ",methods)))
          tab <- xtable(error_table,
                        digits = 3,
                        caption = paste("Rejection rates for" ,
                                        "hypothesis",
                                        "$\\mathcal H_{0,",
                                        which(hyp == dimnames(error[[j]])[[3]]),
                                        "}$ of the",
                                        paste0(gsub("GrandMean","Grand-mean",names(error)[j]),"-type"),
                                        "contrast matrix with",
                                        #dimnames(error[[j]])[[4]][coor],
                                        #"coordinates, 
                                        "$\\delta =$",
                                        d,
                                        "and with sample size",
                                        n_temp))
          align(tab) <- c("r", "l", "l", rep("c",length(align(tab))-3) )
          print(tab, include.rownames = FALSE, sanitize.text.function = function(x) {x},
                sanitize.colnames.function = function(x) {x})
        }}
    }
  }
  }}, file = "C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results_new/Multiple_Tables.txt")


#######################
###     H1 plots    ###
#######################

err_alt <- lapply(error, function(x) x[as.numeric(sapply(strsplit(names2, split ="_"), function(x) x[2])) != 0,,,,drop=FALSE])
names(err_alt) <- names(error)
for(j in 1:length(err_alt)){
  if(j <= 2) f_set <- which(dimnames(err_alt[[j]])[[3]] %in% c("4 - 1", "4 - 2", "4 - 3", "4"))
  if(j == 3) f_set <- 1:4
  if(save_plots) tiff(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results_new/",names(err_alt)[j],"_H1_rejection_rate.tiff"), 
                     width = ifelse(j < 3, 7,10), height = 7, units = "in", res = 600)
  if(save_plots_pdf) pdf(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results_new/",names(err_alt)[j],"_H1_rejection_rate.pdf"), 
                      width = ifelse(j < 3,ifelse(j<2,3.5,4.5),5),height = ifelse(j < 3,ifelse(j<2,3.5,3),2.5))
  par(mfrow=c(1,length(f_set)),mar = c(8.1,4.1,2.1,2.1))
  for(f in f_set){
    my_dataframe <- as.data.frame(err_alt[[j]][,,f,])
    if(length(coordinates) == 2){
      ord_vec <- c(matrix(c(1:(ncol(my_dataframe)/2),(ncol(my_dataframe)/2+1):(ncol(my_dataframe))),ncol=ncol(my_dataframe)/2,byrow=T))
    } else{
      ord_vec <- 1:(ncol(my_dataframe))
    }
    boxplot(my_dataframe[,ord_vec],
            main = bquote(paste(H[paste(0, ",", .(f))])),
            las = 2, names = rep(dimnames(err_alt[[j]][,,,])[[2]],each=length(coordinates)),
            col = ifelse(ord_vec <= ncol(my_dataframe)/2, color[1], color[length(color)]),
            ylim = c(0,max(err_alt[[j]], na.rm = TRUE)+0.03), ylab ="rejection rate")
    #abline(h = 0.05, col="red")
    if(identical(coordinates, c("equi", "inequi"))) legend("bottomleft", c("equicoordinate", "inequicoordinate"), col = color, pch = 16, bty = "n")
  }
  if(save_plots | save_plots_pdf) dev.off()
}

#######################
###   H1 Asymptotic ###
#######################

my_ind <- 3
split_liste <- strsplit(dimnames(err_alt[[1]])[[1]], split = "_")
split_mat <- sapply(split_liste, function(x) if(length(x) > 7){c(x[2:4],paste(x[5],x[6]),x[7:8])}else{x[-1]})



for(j in 1:length(error)){
  
  if(save_plots) tiff(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results_new/",names(err_alt)[j],"_H1_asymptotic.tiff"), 
                     width = ifelse(j < 3, 7,10), height = ifelse(j < 3, 7,10), units = "in", res = 600)
  if(save_plots_pdf) pdf(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results_new/",names(err_alt)[j],"_H1_asymptotic.pdf"), 
                      width = ifelse(j < 3, 7,10), height = ifelse(j < 3, 7,10))
  
  f_set <- which(dimnames(err_alt[[j]])[[3]] %in% c("4 - 1", "4 - 2", "4 - 3", "4"))
  if(j <= 2) f_set <- which(dimnames(err_alt[[j]])[[3]] %in% c("4 - 1", "4 - 2", "4 - 3", "4"))
  if(j == 3) f_set <- 1:4
  laenge <- length(unique(split_mat[my_ind,]))
  par(mfrow=c(length(f_set),laenge),mar = c(8.1,4.1,2.1,2.1))
  
  for(f in f_set){
    for(setup in unique(split_mat[my_ind,])){
      my_err <-  err_alt[[j]][split_mat[my_ind,]==setup,,,,drop = FALSE]
      my_dataframe <- as.data.frame(my_err[,,f,])
      if(length(coordinates) == 2){
        ord_vec <- c(matrix(c(1:(ncol(my_dataframe)/2),(ncol(my_dataframe)/2+1):(ncol(my_dataframe))),ncol=ncol(my_dataframe)/2,byrow=T))
      } else{
        ord_vec <- 1:(ncol(my_dataframe))
      }
      boxplot(my_dataframe[,ord_vec],
              main = bquote(paste(H[paste(0, ",", .(f))], " with ", .(setup), " sample sizes")),
              las = 2, names = rep(dimnames(my_err[,,,])[[2]],each=length(coordinates)),
              col = ifelse(ord_vec <= ncol(my_dataframe)/2, color[1], color[length(color)]),
              ylim = c(0,max(err_alt[[j]], na.rm = TRUE)+0.1), ylab ="rejection rate")
      #abline(h = 0.05, col="red")
      if(identical(coordinates, c("equi", "inequi"))) if(setup == unique(split_mat[my_ind,])[laenge]) legend("topright", c("equicoordinate", "inequicoordinate"), col = color, pch = 16, bty = "n")
    } 
  }
  if(save_plots | save_plots_pdf) dev.off()
}



