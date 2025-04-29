
### Example file how to start a simulation ###

# Which setting? # vector of length 2, e.g.,
setting <- c(1, 1)

# set working directory
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# set the number of cores for parallel computing
cores <- 16
# cores <- 1 # on Windows

# source simulation parameters and start to simulate
# the main simulation
source("source.R")
# the bonferroni-adjusted methods of the main simulation
source("source_bonf.R")
# the additional simulation for analyzing asymptotic behavior
source("source_large.R")

