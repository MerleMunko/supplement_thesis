
### Example file how to start a simulation ###

# Which setting? # value in {1,...,9}, e.g.,
setting <- 1

# set the number of cores for parallel computing
cores <- 16
# cores <- 1 # on Windows
## WARNING: The results are only the same as in our simulations for 16 cores

# set working directory
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# source simulation parameters and start to simulate
source("source.R")