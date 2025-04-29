
### View the results

library(xtable)

# clear the workspace
#rm(list=ls())

# Set the working directory
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# get all result file names
names_cens <- list.files(paste0(path,"/censrate"))
names1 <- gsub('.{6}$', "", names_cens)
names3 <- sapply(strsplit(names1, split ="_"), 
                 function(x)gsub('.{1}$', "",paste0(x[-c(1,3)],collapse = "",sep="_")))
names_error <- list.files(paste0(path,"/errorrate"))
names2 <- gsub('.{6}$', "", names_error)

# load the data of Setting 1
load(paste0(path,"/censrate/",names_cens[1]))
load(paste0(path,"/errorrate/",names_error[1]))

# create a result array for the error- and cens.-rates
cens <- vector(mode = "list", length = length(names2))
names(cens) <- names3
error <- array(NA, dim = c(length(names2),dim(results)))
dimnames(error) <- list(names2, dimnames(results)[[1]], dimnames(results)[[2]],dimnames(results)[[3]])

for(i in 1:length(names_error)){
  # load the results
  load(paste0(path,"/censrate/",names_cens[i]))
  load(paste0(path,"/errorrate/",names_error[i]))
  
  
  error[i,,,] <- results
  cens[[i]] <- censp_s
}

# only consider the methods of interest
methods <- c("asymptotic_global","permutation","asymptotic","wild, Rademacher","wild, Gaussian","groupwise")
coordinates <- c("inequi")
merken <- dim(error)
error_new <- error[,methods,,coordinates]
dim(error_new) <- c(merken[1],length(methods),merken[3],length(coordinates))
dimnames(error_new) <- list(dimnames(error)[[1]], methods, dimnames(error)[[3]], coordinates)
error <- error_new

#######################
###      Tables     ###
#######################

# view the censoring rates
cens_table <- t(as.data.frame(cens[unique(names3)]))
#rownames(cens_table) <- NULL
colnames(cens_table) <- paste("group", 1:ncol(cens_table))
gesplittet <- (strsplit(unique(names3), split = "_"))
delta <- sapply(gesplittet, function(x) x[1])
distribution <- sapply(gesplittet, function(x) paste(x[2:(length(x)-2)], collapse=" "))
censoring_distribution <- sapply(gesplittet, function(x) paste(x[(length(x)-1):length(x)], collapse = " "))
cens_table <- data.frame(delta, distribution, censoring_distribution, round(cens_table,digits=2))
print(xtable(cens_table,
             digits = 2,
             caption = "Censoring rates for the different settings"),
      include.rownames = FALSE)
# view the error rates
#View(as.data.frame(error))
grenzen <- qbinom(c(0.025,0.975),5000,0.05)/5000
gesplittet <- (strsplit(names2, split = "_"))
delta <- sapply(gesplittet, function(x) x[1])
n <- sapply(gesplittet, function(x) paste(x[2], collapse=" "))
distribution <- sapply(gesplittet, function(x) paste(x[3:(length(x)-2)], collapse=" "))
censoring_distribution <- sapply(gesplittet, function(x) paste(x[(length(x)-1):length(x)], collapse = " "))
for(d in unique(delta)){for(j in 1:dim(error)[3]){
  for(coor in 1:dim(error)[4]){
    for(n_temp in unique(n)){
    error_table <- error[n== n_temp & delta==d,,j,coor]
    error_table[error_table >= grenzen[1] & error_table <= grenzen[2]] <- paste("\\textcolor{red}{",sprintf("%.3f",error_table[error_table >= grenzen[1] & error_table <= grenzen[2]]),"}")
    error_table[!is.na(as.numeric(error_table))] <- sprintf("%.3f",as.numeric(error_table[!is.na(as.numeric(error_table))]))
    error_table <- data.frame( distribution[n== n_temp & delta==d], censoring_distribution[n== n_temp & delta==d], error_table)
    colnames(error_table) <- c("distribution", "censoring distribution", gsub("_"," ",methods))
    print(xtable(error_table,
                 digits = 3,
                 caption = paste("Rejection rates for the",
                                 dimnames(error)[[3]][j],
                                 "contrast matrix with",
                                 dimnames(error)[[4]][coor],
                                 "coordinates, $\\delta =$",
                                 d,
                                 "and with sample size",
                                 n_temp)),
          include.rownames = FALSE, sanitize.text.function = function(x) {x})
    }}
  }
}


#######################
###     H0 plots    ###
#######################

#ord_vec <- c(matrix(c(1:dim(error)[2],(dim(error)[2]+1):(2*dim(error)[2])),ncol=dim(error)[2],byrow=T))
err_null <- error[as.numeric(sapply(strsplit(names2, split ="_"), function(x) x[1])) == 0,,,]
dim(err_null)[4] <- 1
for(j in 1:dim(error)[3]){
  tiff(paste0(dimnames(error)[[3]][j],"_rejection_rate.tiff"), res = 800)
  par(mfrow=c(1,1),mar = c(8.1,4.1,2.1,2.1), xpd = FALSE)
boxplot(as.data.frame(err_null[,,j,]), #[,ord_vec],
        las = 2, names = dimnames(error)[[2]],
        col = "palegreen1",
        ylim = c(0,max(err_null)+0.03), ylab ="rejection rate")
abline(h = 0.05, col="red")
abline(h = grenzen[1],col="red", lty="dashed")
abline(h = grenzen[2],col="red", lty="dashed")
#legend("topleft", c("equicoordinate", "inequicoordinate"), col = c("lightskyblue1", "palegreen1"), pch = 16, bty = "n")
dev.off()
}


#######################
###   Power plots   ###
#######################

par(mfrow = c(1,2))
names4 <- substr(names2, 5, nchar(names2))
x <- unique(as.numeric(substr(names2, 1,3)))

palette <- c("red", "orange", "black", "green", "blue", "purple")
for(nam in unique(names4)){
  ind <- which(names4 == nam)
  for(k in 1:dim(error)[3]){
  m <- 1
  if(length(error[ind,1,k,1])==length(x)){
  plot(x,error[ind,1,k,1],type="l",col=palette[1], ylim=c(0,1),
       ylab = "rejection rate", xlab = expression(delta), main = dimnames(error)[[3]][k])
  for(l in 2:dim(error)[2]){
    m <- m+1
    lines(x,error[ind,l,k,1],col=palette[m])
  }
  legend("topleft",col=palette, dimnames(error)[[2]], pch=16, bty = "n")}
  }}

