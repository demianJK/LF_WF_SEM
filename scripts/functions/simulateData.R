### this function simulates data ...
# - in LF 
# - separately for the between and within level 
# - from population covariance matrices 

simulateData <- function(N, n, g, p, var_B, var_W) { 
  
  ## (a) variances
  
  var_w_All <- c()
  var_b_All <- c()
  for (i in 1:p) {
    var_w_All[i] <- paste0("x", i, "~~", var_W, "*", "x", i)
    var_b_All[i] <- paste0("x", i, "~~", var_B, "*", "x", i)
  }
  var_w_All <- paste(var_w_All, collapse = "; ")
  var_b_All <- paste(var_b_All, collapse = "; ")
  
  
  ## (b) covariances 
  
  # The correlation of the two variables is the same at both levels (.3).
  # We get the covariance at each level by transforming the formula for the correlation.
  # corr_x1x2 = cov_x1x2 / (sd_x1 * var_x2)    | * (sd_x1 * sd_x2) and sd_x1 = sd_x2 thus | * var_x1
  # corr_x1x2 * var_x1 = cov_x1x2
  
  cov_w <- 0.3 * var_W 
  cov_b <- 0.3 * var_B
  covs_w <- c()
  covs_b <- c()
  count = 0
  for (i in 1:p) {
    for (j in 1:p) {
      if (i != j & j > i) {
        count <- count + 1
        covs_w[count] <- paste("x", i, "~~", cov_w, "*", "x", j, sep = "")
        covs_b[count] <- paste("x", i, "~~", cov_b, "*", "x", j, sep = "")
      }
    }
  }
  covs_w <- paste(covs_w, collapse = "; ")
  covs_b <- paste(covs_b, collapse = "; ")
  
  
  ## (c) means
  
  # means are 0 per default
  
  
  # put all parts (a), (b), (c) together
  popModel_W <- paste(var_w_All, covs_w, sep = ";")
  popModel_B <- paste(var_b_All, covs_b, sep = ";")
  
  # sample data
  sample_B <- simulateData(popModel_B, sample.nobs = g, model.type = "lavaan")
  sample_W <- simulateData(popModel_W, sample.nobs = N, model.type = "lavaan")
  groups <- rep(1:g, each = n) # group numbers ("j" in Fig. 1)
  data <- sample_W # create data frame with the same dimensions
  data[,] <- 0 # .. and clear all entries
  for (j in unique(groups)) {
    for (i in min(which(groups == j)):max(which(groups == j)))
      data[i,] <- sample_W[i,] + sample_B[j,]
  }
  data$groups <- as.factor(groups)
  data$persons <- rep(1:n, g) # unit numbers ("i" in Fig. 1)
  data <- cbind(data[,(p+1):(p+2)], data[,1:p]) # rearrange columns

  return(data)
}