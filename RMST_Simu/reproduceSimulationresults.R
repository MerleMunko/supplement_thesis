# File to reproduce all results

# set working directory
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# set the number of cores for parallel computing
cores <- 16
# cores <- 1 # on Windows
## WARNING: The results are only the same as in our simulations for 16 cores

# source simulation parameters and start to simulate
for(set1 in c(1,2,3,4,6,7,8,9,11)){
  for(set2 in 1:3){
    setting <- c(set1, set2)
    # the main simulation
    source("source.R")
    # the bonferroni-adjusted methods of the main simulation
    source("source_bonf.R")
    # the additional simulation for analyzing asymptotic behavior
    source("source_large.R")
  }
}

# the simulation under non-exchangeability
source("new_Simu.R")
# the simulation based on the data example
source("data_Simulation.R")
