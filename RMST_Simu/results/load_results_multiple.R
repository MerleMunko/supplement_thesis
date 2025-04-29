
### View the results for the local hypotheses
library(xtable)
library(abind)

# clear the workspace
#rm(list=ls())

# Set the working directory
path <- paste0(dirname(rstudioapi::getSourceEditorContext()$path), 
               "/multiple_errorrate/")
setwd(path)

# should the plots being saved (or just displayed?)
save_plots <- FALSE
# in pdf?
save_plots_pdf <- TRUE

# get all result file names
names_error <- list.files(paste0(path))
names2 <- gsub('.{6}$', "", names_error)

# load the data of Setting 1
load(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results/bonf/multiple_errorrate/",names_error[1]))
anz_bonfs <- dim(results_multiple[[1]])[1]
names_bonfs <- dimnames(results_multiple[[1]])[[1]]
load(paste0(path,names_error[1]))
anz_nonbonfs <- dim(results_multiple[[1]])[1]

# create a result array for the errorrates
error <- vector(mode = "list", length = length(results_multiple))
for(j in 1:length(results_multiple)){
  error[[j]] <- array(NA, dim = c(length(names2),dim(results_multiple[[j]])+c(anz_bonfs,0,0)),
                      dimnames = list(names2,c(dimnames(results_multiple[[j]])[[1]], names_bonfs),
                                   dimnames(results_multiple[[j]])[[2]], 
                                   dimnames(results_multiple[[j]])[[3]]))
}

for(i in 1:length(names_error)){for(j in 1:length(results_multiple)){
  # load the results
  load(paste0(path,names_error[i]))
  # save the multiple results
  error[[j]][i,1:anz_nonbonfs,,] <- results_multiple[[j]]
  
  # load the bonferroni results
  try({
  load(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results/bonf/multiple_errorrate/",names_error[i]))
  # save the multiple results
  error[[j]][i,(anz_nonbonfs+1):(anz_nonbonfs+anz_bonfs),,] <- results_multiple[[j]]
  })
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
for(j in 1:length(error)){
  for(n_i in rev(unique(n))){
    for(cens_i in rev(unique(censoring_distribution))){
      # weib
      namesj <- dimnames(error[[j]])[[1]]
      merge_ind <- which(namesj %in% paste0("multiple_0.0_", n_i, c("_Weib_late_", "_Weib_prop_"), cens_i   ))
      merge <- colMeans(error[[j]][merge_ind,,,])
      dim(merge) <- dim(error[[j]][merge_ind[1],,,])
      error[[j]][merge_ind[1],,,] <- merge
      error[[j]] <- error[[j]][-merge_ind[-1],,,,drop=FALSE]
      #error[[j]] <- abind(merge , error[[j]], along = 1)
      dimnames(error[[j]])[[1]][merge_ind[1]] <- paste0("multiple_0.0_", n_i, "_Weib_late,prop_", cens_i   )
    }
  } 
for(n_i in rev(unique(n))){
  for(cens_i in rev(unique(censoring_distribution))){
    # exp
    namesj <- dimnames(error[[j]])[[1]]
    merge_ind <- which(namesj %in% paste0("multiple_0.0_", n_i, c("_exp_early_", "_exp_late_", "_exp_prop_"), cens_i   ))
    merge <- colMeans(error[[j]][merge_ind,,,])
    dim(merge) <- dim(error[[j]][merge_ind[1],,,])
    error[[j]][merge_ind[1],,,] <- merge
    error[[j]] <- error[[j]][-merge_ind[-1],,,,drop=FALSE]
    #error[[j]] <- abind(merge , error[[j]], along = 1)
    dimnames(error[[j]])[[1]][merge_ind[1]] <- paste0("multiple_0.0_", n_i, "_exp_early,late,prop_", cens_i   )
    }
}
  
  }
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
}}, file = "C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results/Multiple_Tables.txt")


#######################
###     H1 plots    ###
#######################

err_alt <- lapply(error, function(x) x[as.numeric(sapply(strsplit(names2, split ="_"), function(x) x[2])) != 0,,,,drop=FALSE])
names(err_alt) <- names(error)
for(j in 1:length(err_alt)){
  if(j <= 2) f_set <- which(dimnames(err_alt[[j]])[[3]] %in% c("4 - 1", "4 - 2", "4 - 3", "4"))
  if(j == 3) f_set <- 1:4
  if(save_plots) tiff(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results/",names(err_alt)[j],"_H1_rejection_rate.tiff"), 
                      width = ifelse(j < 3, 7,10), height = 7, units = "in", res = 600)
  if(save_plots_pdf) pdf(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results/",names(err_alt)[j],"_H1_rejection_rate.pdf"), 
                      width = ifelse(j < 3, 7,10))
  
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
  
  if(save_plots) tiff(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results/",names(err_alt)[j],"_H1_asymptotic.tiff"), 
                      width = 4.5, height = ifelse(j < 3, ifelse(j<2,2,6),8), units = "in", res = 600)
  if(save_plots_pdf) pdf(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results/",names(err_alt)[j],"_H1_asymptotic.pdf"), 
                      width = 4.5, height = ifelse(j < 3, ifelse(j<2,2,6),8))
  
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
              main = bquote(paste(H[paste(0, ",", .(f))], ", ", .(setup))),
              las = 2, names = rep(dimnames(my_err[,,,])[[2]],each=length(coordinates)),
              col = ifelse(ord_vec <= ncol(my_dataframe)/2, color[1], color[length(color)]),
              ylim = c(0,max(err_alt[[j]], na.rm = TRUE)+0.1), ylab ="rejection rate")
      #abline(h = 0.05, col="red")
      if(identical(coordinates, c("equi", "inequi"))) if(setup == unique(split_mat[my_ind,])[laenge]) legend("topright", c("equicoordinate", "inequicoordinate"), col = color, pch = 16, bty = "n")
    } 
  }
  if(save_plots | save_plots_pdf) dev.off()
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
    f_set <- which(dimnames(error[[k]])[[3]] %in% c("4 - 1", "4 - 2", "4 - 3", "4"))
    par(mfrow = c(1,length(f_set)))
    for(f in f_set){
    m <- 1
    if(length(error[[k]][ind,1,f,2])==length(x)){
      plot(x,error[[k]][ind,1,f,2],type="l",col=palette[1], ylim=c(0,1),
           ylab = "rejection rate", xlab = expression(delta), main = dimnames(error[[k]])[[3]][f])
      for(l in 2:dim(error[[k]])[2]){
        m <- m+1
        lines(x,error[[k]][ind,l,f,2],col=palette[m])
      }
      legend("topleft",col=palette, paste("\n", dimnames(error[[k]])[[2]],"\n"), pch=16, bty = "n")}
  }}}
# achtung: nur fuer inequi coordinaten
# auch lieber eine Tabelle hier...

##########################
###   Power barplots   ###
##########################

#rm(list=ls())
# only consider the methods of interest 
# do not consider global approaches!
methods <- c("asymptotic","wild, Rademacher","wild, Gaussian","groupwise","asymptotic_bonf","permutation_bonf")
coordinates <- c("equi", "inequi")
# treffe eine auswahl an plots
#auswahl <- 1:90
auswahl <- 1


# load all decision_arr objects
# Set the working directory
path <- "C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results/decision_list/"
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
  # load the bonferroni data
  load(paste0("C:/Users/munko/Documents/Survival/RMST/RMST_RCode/results/bonf/decision_list/",name_auswahl))
  bonf_decision_arr <- decision_arr
  # load the non-bonferroni data
  load(name_auswahl)
  nonbonf_decision_arr <- decision_arr
  nonbonf_methods <- dimnames(nonbonf_decision_arr[[1]])[[1]]
  
  for(meth in methods){ 
    if(meth %in% nonbonf_methods){
      decision_arr <- nonbonf_decision_arr
    }else{
      decision_arr <- bonf_decision_arr
      }
    for(coor in coordinates){
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
    
      box_heights <- t(apply(dec[[k]][ind[order(x)],meth,,2], 1, cumsum))
      barplot(box_heights[,ncol(box_heights)], col = palette[ncol(box_heights)], 
              border = FALSE, names.arg = sort(x), xlab = expression(delta), ylim = c(0,1), main = meth)
      for(column in (ncol(box_heights)-1):1){
      barplot(box_heights[,column], add=TRUE, col = palette[column], border = FALSE, names.arg = sort(x))
      }
    }
    
#    if(dim(dec[[k]])[3]<= 15){
      plot.new()
      text(0.5,1, labels = "Which hypotheses are rejected? ")
      legend(0,0.9,col=palette, paste(dimnames(dec[[k]])[[3]]), pch=16, bty = "n")
#    }
    # if(dim(dec[[k]])[3]> 15){
    #   par(mfrow = c(1,1))
    #   plot.new()
    #   text(0.5,1, labels = "Which hypotheses are rejected? ")
    #   legend(0,0.9,col=palette, paste(dimnames(dec[[k]])[[3]]), pch=16, bty = "n")
    #   par(mfrow = c(1,length(methods)+1))
    #   }
  }}

# achtung: nur fuer inequi coordinaten 
# lieber tabelle... und ein Beispiel.