# File to reproduce all results

# set working directory
path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

# set the number of cores for parallel computing
cores <- 16
# cores <- 1 # on Windows
## WARNING: The results are only the same as in our simulations for 16 cores

# source simulation parameters and start to simulate
for(setting in 1:9){
  # reset random number generation
  RNGversion("4.4.1")
  RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  
  source("source.R")
}
