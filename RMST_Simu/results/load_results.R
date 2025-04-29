
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
save_plots_pdf <- TRUE

# get all result file names
names_cens <- list.files(paste0(path,"/censrate"))
names1 <- gsub('.{6}$', "", names_cens)
names3 <- sapply(strsplit(names1, split ="_"), 
                 function(x)gsub('.{1}$', "",paste0(x[-c(1,3,4)],collapse = "",sep="_")))
names_error <- list.files(paste0(path,"/errorrate"))
names2 <- gsub('.{6}$', "", names_error)

# load the data of Setting 1
load(paste0(path,"/bonf/errorrate/",names_error[1]))
anz_bonfs <- dim(results)[1]
names_bonfs <- dimnames(results)[[1]]
load(paste0(path,"/censrate/",names_cens[1]))
load(paste0(path,"/errorrate/",names_error[1]))
anz_nonbonfs <- dim(results)[1]

# create a result array for the error- and cens.-rates
cens <- vector(mode = "list", length = length(names2))
names(cens) <- names3
error <- array(NA, dim = c(length(names2),dim(results)+c(anz_bonfs,0,0))) # add the number of bonferroni methods
dimnames(error) <- list(names2, c(dimnames(results)[[1]], names_bonfs), 
                        dimnames(results)[[2]],dimnames(results)[[3]])

for(i in 1:length(names_error)){
  # load the results
  load(paste0(path,"/censrate/",names_cens[i]))
  load(paste0(path,"/errorrate/",names_error[i]))
  
  error[i,1:anz_nonbonfs,,] <- results
  cens[[i]] <- censp_s
  
  # add bonferroni methods
  try({load(paste0(path,"/bonf/errorrate/",names_error[i]))
  error[i,(anz_nonbonfs+1):(anz_nonbonfs+anz_bonfs),,] <- results})

}

# only consider the methods of interest
# methods <- c("asymptotic_global","permutation","asymptotic","wild, Rademacher","wild, Gaussian","groupwise")
methods <- c("asymptotic_global","permutation","asymptotic","wild, Rademacher","wild, Gaussian","groupwise",
             "asymptotic_bonf","permutation_bonf")
coordinates <- c("inequi")
color <- c("darkgrey")
length_coor <- length(coordinates)
merken <- dim(error)
error_new <- error[,methods,,coordinates]
dim(error_new) <- c(merken[1],length(methods),merken[3],length_coor)
dimnames(error_new) <- list(dimnames(error)[[1]], methods, dimnames(error)[[3]], coordinates)
error <- error_new

# combine the equal settings, i.e. exp early, late & prop under the null
# and Weib late & prop under the null
gesplittet <- (strsplit(names2, split = "_"))
delta <- sapply(gesplittet, function(x) x[1])
n <- sapply(gesplittet, function(x) paste(x[2:3], sep=" ", collapse = "_"))
distribution <- sapply(gesplittet, function(x) paste(x[4:(length(x)-2)], sep=" ", collapse = "_"))
censoring_distribution <- sapply(gesplittet, function(x) paste(x[(length(x)-1):length(x)], sep=" ", collapse = "_"))

for(n_i in rev(unique(n))){
  for(cens_i in rev(unique(censoring_distribution))){
    # weib
    merge_ind <- which(names2 %in% paste0("0.0_", n_i, c("_Weib_late_", "_Weib_prop_"), cens_i   ))
    merge <- colMeans(error[merge_ind,,,])
    dim(merge) <- dim(error[merge_ind[1],,,])
    error[merge_ind[1],,,] <- merge
    error <- error[-merge_ind[-1],,,,drop=FALSE]
    #error <- abind(merge , error, along = 1)
    dimnames(error)[[1]][merge_ind[1]] <- paste0("0.0_", n_i, "_Weib_late,prop_", cens_i   )
    names2 <- dimnames(error)[[1]]
  }
}

for(n_i in rev(unique(n))){
  for(cens_i in rev(unique(censoring_distribution))){
    # exp
    merge_ind <- which(names2 %in% paste0("0.0_", n_i, c("_exp_early_", "_exp_late_", "_exp_prop_"), cens_i   ))
    merge <- colMeans(error[merge_ind,,,,drop=FALSE])
    error <- error[-merge_ind,,,,drop=FALSE]
    error <- abind(merge , error, along = 1)
    dimnames(error)[[1]][dimnames(error)[[1]] == ""] <- paste0("0.0_", n_i, "_exp_early,late,prop_", cens_i   )
    names2 <- dimnames(error)[[1]]
  }}

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
n <- sapply(gesplittet, function(x) paste(x[2:3], collapse=" "))
distribution <- sapply(gesplittet, function(x) paste(x[4:(length(x)-2)], collapse=" "))
censoring_distribution <- as.factor(sapply(gesplittet, function(x) paste(x[(length(x)-1):length(x)], collapse = " ")))
# rename the censoring distribution
levels(censoring_distribution) <- c("equal", "unequal, high", "unequal, low")
for(d in unique(delta)){capture.output({
for(j in 1:dim(error)[3]){
  for(coor in 1:dim(error)[4]){
    for(n_temp in unique(n)){
    error_table <- error[n== n_temp & delta==d,,j,coor]
    error_table[error_table >= grenzen[1] & error_table <= grenzen[2]] <- paste0("\\textbf{",sprintf("%.3f",error_table[error_table >= grenzen[1] & error_table <= grenzen[2]]),"}")
    error_table[!is.na(as.numeric(error_table))] <- sprintf("%.3f",as.numeric(error_table[!is.na(as.numeric(error_table))]))
    error_table <- data.frame( distribution[n== n_temp & delta==d], censoring_distribution[n== n_temp & delta==d], error_table)
    colnames(error_table) <- c("distribution", "censoring \\newline distribution", gsub(",","",gsub("_"," ",methods)))
    tab <- xtable(error_table,
                  digits = 3,
                  caption = paste("Rejection rates for the",
                                  paste0(gsub("GrandMean","Grand-mean",dimnames(error)[[3]][j]),"-type"),
                                  "contrast matrix with",
                                  #dimnames(error)[[4]][coor],
                                  #"coordinates, 
                                  "$\\delta =$",
                                  d,
                                  "and sample size",
                                  n_temp))
    align(tab) <- c("r", "l",  "p{0.6in}", rep(">{\\centering\\arraybackslash}p{0.45in}",length(align(tab))-3) )
    print(tab, include.rownames = FALSE, sanitize.text.function = function(x) {x},
          sanitize.colnames.function = function(x) {x})
    }}
  }
}, file = paste0("Tables",d,".txt"))}

# barplot highest rejection rates
#barplot(table(factor(methods[apply(matrix(aperm(error[delta==d,,,], c(2,1,3,4)), ncol = 8, byrow = TRUE),1, which.max )], levels = methods)))

######################
###     H0 plots    ###
#######################

if(length_coor == 2){
  ord_vec <- c(matrix(c(1:dim(error)[2],(dim(error)[2]+1):(2*dim(error)[2])),ncol=dim(error)[2],byrow=T))
}else{
  ord_vec <- c(1:dim(error)[2])
  
}
err_null <- error[as.numeric(sapply(strsplit(names2, split ="_"), function(x) x[1])) == 0,,,,drop=FALSE]
for(j in 1:dim(error)[3]){
  if(save_plots) tiff(paste0(dimnames(error)[[3]][j],"_rejection_rate.tiff"), 
                      units = "in", width =7, height = 7, res = 600)
  if(save_plots_pdf) pdf(paste0(dimnames(error)[[3]][j],"_rejection_rate.pdf"))
  par(mfrow=c(1,1),mar = c(8.1,4.1,2.1,2.1), xpd = FALSE)
  boxplot(as.data.frame(err_null[,,j,])[,ord_vec],
        las = 2, names = rep(dimnames(error)[[2]],each=length_coor),
        col = ifelse(ord_vec <= dim(error)[2], color[1], color[length(color)]),
        ylim = c(0,max(err_null, na.rm = TRUE)), ylab ="rejection rate")
  abline(h = 0.05, lty = "dotted")
  abline(h = grenzen[1], lty="dashed")
  abline(h = grenzen[2], lty="dashed")
  if(length_coor > 1) legend("topleft", c("equicoordinate", "inequicoordinate"), col = color, pch = 16, bty = "n")
  if(save_plots | save_plots_pdf) dev.off()
}



#######################
### Asymptotik plot ###
#######################

split_liste <- strsplit(dimnames(err_null)[[1]], split = "_")
split_mat <- sapply(split_liste, function(x) if(length(x) > 6){c(x[1:3],paste(x[4],x[5]),x[6:7])}else{x})

my_ind <- 3
for(j in 1:dim(error)[3]){
  if(save_plots) tiff(paste0(dimnames(error)[[3]][j],"_asymptotic_plots.tiff"), 
                      width = 6, height = 4, units = "in", res = 600)
  if(save_plots_pdf) pdf(paste0(dimnames(error)[[3]][j],"_asymptotic_plots.pdf"), 
                      width = 6, height = 4)
  par(mfrow=c(1,length(unique(split_mat[my_ind,]))),mar = c(8.1,4.1,2.1,2.1), xpd = FALSE)
  for(setup in unique(split_mat[my_ind,])){
    my_err <-  err_null[split_mat[my_ind,]==setup,,,,drop=FALSE]
    
    boxplot(as.data.frame(my_err[,,j,])[,ord_vec],main=setup,
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



## alter Asymptotik Plot
#ohne_groesse <- sapply(strsplit(names2[as.numeric(sapply(strsplit(names2, split ="_"), function(x) x[1])) == 0], split ="_"), function(x)x[-3])
groesse <- sapply(strsplit(names2[as.numeric(sapply(strsplit(names2, split ="_"), function(x) x[1])) == 0], split ="_"), function(x)x[3])
#par(mfrow=c(1,1),mar = c(8.1,4.1,2.1,2.1))
for(k in 1:dim(error)[3]){
  #plot(c(60,120,240), numeric(3), ylim = c(0,max(err_null)+0.1), type="n", ylab = "Ablehnrate",
  #     xlab = "n")
  
#for(j in 1:length(unique(ohne_groesse))){
  #namen <- paste("0.0" , unique(ohne_groesse)[[j]][2], c("small", "medium", "large"), paste(  unique(ohne_groesse)[[j]][-c(1,2)], collapse = " "))
  #namen <- gsub( " ", "_", namen)
  #welche <- namen %in% names2
  #X <- matrix(NA, ncol = 2, nrow = sum(welche))
  #X[,1] <- c(60,120,240)[welche]
  #error[namen]
  #farben <- rainbow(dim(error)[2])
  #pdf(paste0(dimnames(error)[[3]][k],"_asymptotic_plots.pdf"), width = 10)
  par(mfrow=c(3,2))
  for(g in 1:dim(error)[2]){
    boxplot(list(small = err_null[groesse == "small",g,k,1], 
                 medium = err_null[groesse == "medium",g,k,1], 
                 large = err_null[groesse == "large",g,k,1],
            small = err_null[groesse == "small",g,k,2], 
            medium = err_null[groesse == "medium",g,k,2], 
            large = err_null[groesse == "large",g,k,2]),
            col = c(rep(color[1],3),rep(color[2],3)),main = dimnames(error)[[2]][g],
            ylim = c(min(err_null),max(err_null)), ylab ="rejection rate", cex.axis = 1.5, cex.lab = 1.5)
    abline(h = 0.05, lty=2)
    text(x = 2, y = max(err_null)-0.02, "equicoordinate", col = color[1], cex = 1.5)
    text(x = 5, y = max(err_null)-0.02, "inequicoordinate", col = color[2], cex = 1.5)
    # boxplot(list(small = err_null[groesse == "small",g,k,2], 
    #              medium = err_null[groesse == "medium",g,k,2], 
    #              large = err_null[groesse == "large",g,k,2]),
    #         col = "palegreen1", main = dimnames(error)[[2]][g],
    #         ylim = c(min(err_null),max(err_null)), ylab ="rejection rate")
    # abline(h = 0.05, lty=2)
    #X[,2] <- error[namen[welche],g,k,1]
    #lines(X, col = farben[g])
    #X[,2] <- error[namen[welche],g,k,2]
    #lines(X, col = farben[g])
  }
#dev.off()
}
  #abline(h = 0.05, lty=2)
  #legend("topright", pch=16, col = farben, dimnames(error)[[2]], cex=0.7)
  #}


#######################
###   Power plots   ###
#######################

par(mfrow = c(1,3))
names4 <- substr(names2, 5, nchar(names2))
x <- unique(as.numeric(substr(names2, 1,3)))

palette <- c("red", "orange", "black", "green", "blue", "purple")
for(nam in unique(names4)){
  ind <- which(names4 == nam)
  for(k in 1:dim(error)[3]){
  m <- 1
  if(length(error[ind,1,k,2])==length(x)){
  plot(x,error[ind,1,k,2],type="l",col=palette[1], ylim=c(0,1),
       ylab = "rejection rate", xlab = expression(delta), main = dimnames(error)[[3]][k])
  for(l in 2:dim(error)[2]){
    m <- m+1
    lines(x,error[ind,l,k,2],col=palette[m])
  }
  legend("topleft",col=palette, paste("\n", dimnames(error)[[2]],"\n"), pch=16, bty = "n")}
  }}

# achtung: nur fuer coordinates ,2 aktuell
# hier vielleicht lieber eine Tabelle.