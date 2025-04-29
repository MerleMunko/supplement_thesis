library(GFDrmtl)
library(mstate)
data("ebmt2")

# data preprocessing
ebmt2 <- ebmt2[ebmt2$tcd != "Unknown",]
ebmt2$tcd <- droplevels(ebmt2$tcd)
my_data <- ebmt2[,c("time","status")]
my_data$status[my_data$status > 2] <- 3
my_data$group[ebmt2$tcd == "No TCD" & ebmt2$match == "No gender mismatch"] <- 1
my_data$group[ebmt2$tcd == "TCD" & ebmt2$match == "No gender mismatch"] <- 2
my_data$group[ebmt2$tcd == "No TCD" & ebmt2$match == "Gender mismatch"] <- 3
my_data$group[ebmt2$tcd == "TCD" & ebmt2$match == "Gender mismatch"] <- 4
colnames(my_data) <- c("X","D","group")

## asymptotic
asymptotic <- RMTL.test(
  time = my_data$X,
  status = my_data$D,
  group = my_data$group,
  hyp_mat = "2by2 cause-wisely",
  M = 3,
  tau = 120,
  method = "asymptotic",
  stepwise = FALSE,
  alpha = 0.05,
  Nres = 19999,
  seed = 1
)

## asymptotic_bonf
# define the hypothesis matrices
M <- 3
HA <- cbind(diag(M), diag(M), -diag(M), -diag(M))
HB <- cbind(diag(M), -diag(M), diag(M), -diag(M))
HAB <- cbind(diag(M), -diag(M), -diag(M), diag(M))
Y <- rbind(HA, HB, HAB)
X <- data.frame(t(Y))
hyp_mat <- as.list(X)
hyp_mat <- lapply(hyp_mat, t)
# go through all hypothesis matrices and calculate p-values
asymptotic_bonf_unadjusted <- sapply(hyp_mat, function(mat){
  RMTL.test(
    time = my_data$X,
    status = my_data$D,
    group = my_data$group,
    hyp_mat = list(mat),
    hyp_vec = list(0),
    M = 3,
    tau = 120,
    method = "asymptotic",
    stepwise = FALSE,
    alpha = 0.05,
    Nres = 19999,
    seed = 1
  )$p.value
})
# adjust with Bonferroni
asymptotic_bonf <- p.adjust(asymptotic_bonf_unadjusted, method = "bonferroni")


## permutation_bonf
permutation <- RMTL.test(
  time = my_data$X,
  status = my_data$D,
  group = my_data$group,
  hyp_mat = "2by2 cause-wisely",
  M = 3,
  tau = 120,
  method = "permutation",
  stepwise = FALSE,
  alpha = 0.05,
  Nres = 19999,
  seed = 1
)

## show estimated RMTLs as in the table
matrix(round(RMTL.test(
  time = my_data$X,
  status = my_data$D,
  group = my_data$group,
  hyp_mat = list(diag(12)),
  hyp_vec = list(numeric(12)),
  M = 3,
  method = "asymptotic",
  tau = 120
)$res[,"estimator"][[1]],3), ncol = 3, byrow = TRUE)

## show resulting p-values as in the table
p_tab <- round(rbind(asymptotic$p.value, asymptotic_bonf, permutation$p.value), 
               3)
colnames(p_tab) <- paste("H", rep(c("A", "B", "AB"), each = 3), rep(1:3, 3))
rownames(p_tab) <- c("asymptotic", "asymptotic_bonf", "permutation_bonf")
p_tab
