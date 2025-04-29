
### View the results

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
save_plots_pdf <- FALSE

# get all result file names
names_error <- list.files(paste0(path,"/errorrate"))
names2 <- gsub('.{6}$', "", names_error)

# load the data of Setting 1
load(paste0(path,"/errorrate/",names_error[1]))


# only consider the methods of interest
methods <- c("asymptotic", "asymptotic_bonf", "perm_bonf")
coordinates <- c("inequi")
color <- c("lightgrey")
color2 <- c(rgb(0,0,0,0), rgb(0,0,0,0))
pch2 <- c(18,16,15)
length_coor <- length(coordinates)

# create a result array for the error- and cens.-rates
error <- array(NA, dim = c(length(names2),length(methods)
                           ,length_coor,dim(results)[3]))
dimnames(error) <- list(names2, methods, coordinates,c("2x2","Dunnett","Tukey"))


for(i in 1:length(names_error)){
  # load the results
  load(paste0(path,"/errorrate/",names_error[i]))
  error[i,1:length(methods),,] <- results[methods,coordinates,]
}
merken <- dim(error)

dimnames(error)[[2]] <- c("asymptotic", "asymptotic_bonf", "permutation_bonf")

grenzen <- qbinom(c(0.025,0.975),2000,0.05)/2000
if(length_coor == 2){
  ord_vec <- c(matrix(c(1:dim(error)[2],(dim(error)[2]+1):(2*dim(error)[2])),ncol=dim(error)[2],byrow=T))
}else{
  ord_vec <- c(1:dim(error)[2])
  
}
err_null <- error[as.numeric(sapply(strsplit(names2, split ="_"), function(x) x[1])) == 0,,,,drop=FALSE]

#######################
### Asymptotik plot ###
#######################

split_liste <- strsplit(dimnames(err_null)[[1]], split = "_")
split_mat <- sapply(split_liste, function(x) if(length(x) > 7){c(x[1:3],paste(x[4],x[5]),x[6:7],x[8])}else{x})
split_mat2 <- apply(split_mat, 2, function(x) paste(x[2],x[3]))

my_ind <- 1 #3
for(j in 1:dim(error)[4]){
  if(save_plots) tiff(paste0(dimnames(error)[[4]][j],"_asymptotic_plots.tiff"), 
                      width = 6, height = 4, units = "in", res = 600)
  if(save_plots_pdf) pdf(paste0(dimnames(error)[[4]][j],"_asymptotic_plots.pdf"), 
                         width = 8, height = 6)
  par(mfrow=c(2,3#1,length(unique(split_mat2)) #split_mat[my_ind,]
    ),mar = c(8.1,4.1,2.1,2.1), xpd = FALSE)
  for(setup in unique(split_mat2 #split_mat[my_ind,]
    )){
    my_err <-  err_null[split_mat2 #split_mat[my_ind,]
                        ==setup,,,,drop=FALSE]
    
    boxplot(as.data.frame(my_err[,,,j])[,ord_vec],main=paste0(setup, " samples"),
            las = 2, names = rep(dimnames(error)[[2]],each=length_coor),
            col = ifelse(ord_vec <= dim(error)[2], color[1], color[length(color)]),
            ylim = c(0,max(err_null, na.rm = TRUE)), ylab ="rejection rate")
    abline(h = 0.05, lty = "dotted")
    abline(h = grenzen[1], lty="dashed")
    abline(h = grenzen[2], lty="dashed")
    if(length_coor > 1) legend("topleft", c("equicoordinate", "inequicoordinate"), col = color, pch = 16, bty = "n")
    
  }
  if(save_plots | save_plots_pdf) dev.off()
}



# the same under the alternative
err_null <- error[as.numeric(sapply(strsplit(names2, split ="_"), function(x) x[1])) == 1.5,,,,drop=FALSE]

#######################
### Asymptotik plot ###
#######################

split_liste <- strsplit(dimnames(err_null)[[1]], split = "_")
split_mat <- sapply(split_liste, function(x) if(length(x) > 7){c(x[1:3],paste(x[4],x[5]),x[6:7],x[8])}else{x})
split_mat2 <- apply(split_mat, 2, function(x) paste(x[2],x[3]))

my_ind <- 1 #3
for(j in 1:dim(error)[4]){
  if(save_plots) tiff(paste0(dimnames(error)[[4]][j],"_power_asymptotic_plots.tiff"), 
                      width = 8, height = 6, units = "in", res = 600)
  if(save_plots_pdf) pdf(paste0(dimnames(error)[[4]][j],"_power_asymptotic_plots.pdf"), 
                         width = 8, height = 6)
  par(mfrow=c(2,3#1,length(unique(split_mat2)) #split_mat[my_ind,]
  ),mar = c(8.1,4.1,2.1,2.1), xpd = FALSE)
  for(setup in unique(split_mat2 #split_mat[my_ind,]
  )){
    my_err <-  err_null[split_mat2 #split_mat[my_ind,]
                        ==setup,,,,drop=FALSE]
    
    boxplot(as.data.frame(my_err[,,,j])[,ord_vec],main=paste0(setup, " samples"),
            las = 2, names = rep(dimnames(error)[[2]],each=length_coor),
            col = ifelse(ord_vec <= dim(error)[2], color[1], color[length(color)]),
            ylim = c(0,max(err_null, na.rm = TRUE)), ylab ="rejection rate")
    #abline(h = 0.05, lty = "dotted")
    #abline(h = grenzen[1], lty="dashed")
    #abline(h = grenzen[2], lty="dashed")
    if(length_coor > 1) legend("topleft", c("equicoordinate", "inequicoordinate"), col = color, pch = 16, bty = "n")
    
  }
  if(save_plots | save_plots_pdf) dev.off()
}

