# function to estimate S_pw (=Sigma_W), using covW(), S_B, using cov(), and Sigma_B using S_pw and S_B

covMulti <- function(data, groups){
  # data = (total) raw data (ie between and within) WITHOUT grouping variable
  # nb of rows ought to respond to total nb of L1 units (i.e. nb of rows = N), cols to METRIC observed variables
  # groups = vector coding groups
  
  library(Morpho)
  
  ### sanity checks:
  
  if (!is.vector(data) & !is.matrix(data) & !is.data.frame(data)){
    stop("data must be a vector, matrix or data frame")
  }
  
  if (is.vector(data)){ # univariate
    if (length(groups) != length(data)) {
      stop("length of groups must correspond to the number of rows in data")
    }
  } else { # multivariate
    if (length(groups) != nrow(data)) {
      stop("length of groups must correspond to the number of rows in data")
    }
  }
  
  if (!is.factor(groups)) { # so no error thrown when using covW()
    groups <- as.factor(groups)
  }
  
  tmp <- table(groups)
  if (min(tmp) != max(tmp)){
    stop("Function works only with balanced group sizes")
  }
  
  ### data preparation
  data <- as.data.frame(data) # necessary for covW()
  n <- min(tmp) # nb of individuals within each group (N_c)
  g <- length(unique(groups)) # nb of groups (c)
  N <- n * g
  p <- ncol(data) # nb of observed variables
  
  
  ### estimation
  
  #  S_pw = Sigma_W 
  
  if (p == 1){ # univariate (duplicate variable otherwise covW() wont work)
    data2 <- cbind(data, data)
    out_covW <- Morpho::covW(data2, groups)
    S_pw <- out_covW[1,1]
    group_means <- attr(out_covW, "means")[,1]
    attr(S_pw, "means") <- NULL
  } else { # multivariate
    S_pw <- Morpho::covW(data, groups)
    group_means <- attr(S_pw, "means") 
    attr(S_pw, "means") <- NULL
  }
  
  
  # S_B and Sigma_B

  grand_means <- colMeans(data, na.rm = TRUE)
  datac <- t(t(group_means) - grand_means) # group means - grand means
  S_B <- lavaan:::lav_matrix_crossprod(datac * n, datac)/(g - 1)
  Sigma_B <- (S_B - S_pw) / n
  
  out <- list(S_pw = S_pw, 
              S_B = S_B, 
              Sigma_B = Sigma_B, 
              group_means = group_means
  )
  return(out)
}
