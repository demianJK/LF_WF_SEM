### this function analyses data in the LF and WF approaches with an intercept-only model ###

analyseData <- function(data_LF, data_WF, p, n, g, center){
  
  ### LF approach
  
  ## 1) set model specification (same for within and between)
  
  # variances
  Var <- c()
  for (i in 1:p) {
    Var[i] <- paste0("x", i, "~~", "x", i)
  }
  Var <- paste(Var, collapse="; ")
  
  # covariances
  Cov <- c()
  count=0
  for(i in 1:p){   
    for(j in 1:p){
      if(i != j & j > i){
        count = count+1
        Cov[count] <- paste("x", i, "~~", "x", j, sep="") 
      }
    }
  }
  Cov <- paste(Cov, collapse="; ")
  
  model_LF_W <- paste0("Level: 1\n", paste(Var, Cov, sep=";"))
  model_LF_B <- paste0("Level: 2\n", paste(Var, Cov, sep=";"))
  model_LF <- paste(model_LF_W, model_LF_B, sep="\n")
  
  # 2) fit model
  
  fit_LF <- try(sem(model = model_LF, # for multilevel SEM, use sem(), multilevel meanstructure: within ints=0, between ints estimated, see https://www.lavaan.ugent.be/tutorial/multilevel.html
                    data = data_LF,
                    cluster = "groups", # based on user-specified model
                    auto.fix.first=FALSE), # only relevant for factor models
                silent = TRUE)
  
  if(inherits(fit_LF,"try-error")){
    fit_LF <- NA 
  }
  
  
  ### WF approach
  
  ## 1) set model specification (different for within and between)
  
  ## within (p * n)
  
  # variances (p * n)
  tmp2 <- c()
  resid_w <- c() # residual variances (p * n - equality among n blocks)
  tmp3 <- c()
  for (j in 1:p){
    for (i in 1:n){
      tmp2[i] <- paste0("x", j, "_", i)
      tmp3[i] <- paste0(tmp2[i], "~~Vx", j, "_w*", tmp2[i])
    }
    resid_w[j] <- paste(tmp3, collapse="; ")
  }
  resid_w <- paste(resid_w, collapse="; ")
  
  # covariances (p-pc=c*n with c-wise equality constraints)
  resid_cov <- c() # manifest correlations (p * n - equality among n blocks)
  count <- 0
  for (i in 1:n){ # n-unit-wise ordering
    for(j in 1:p){   
      for(m in 1:p){
        if(j != m & m > j){
          count <- count + 1
          resid_cov[count] <- paste0("x", j, "_", i, "~~", "Cx", j, m, "_w*", "x", m, "_", i)
        }
      } 
    }
  }
  resid_cov <- paste(resid_cov, collapse="; ")
  
  # "If variables are both at the within- and the between-level, the intercepts at the within-level should be fixed to zero." (p.701f, Barendse & Rosseel, 2020)
  fac_int_w <- c()
  tmp <- c()
  count <- 0
  for (j in 1:p){
    for (i in 1:n){
      count <- count + 1
      tmp[ind] <- paste0("x", j, "_", i, "~0*1")
    }
  }
  fac_int_w <- paste(tmp, collapse = "; ")
  
  model_WF_W <- paste(resid_w, resid_cov, fac_int_w, sep = "; ")
  
  ## between (p)
  
  fac_b <- c() # latent factor = random intercepts (p - loadings fixed to 1)
  tmp <- c()
  for (j in 1:p){
    for (i in 1:n){
      tmp[i] <- paste0("1*x", j, "_", i)
    }
    fac_b[j] <- paste0("fx", j, "=~", paste(tmp, collapse="+"))
  }
  fac_b <- paste(fac_b, collapse="; ")
  
  # variances
  fac_var_b <- c() # factor variance (p)
  fac_int_b <- c() # to estimate means
  for (j in 1:p){
    fac_var_b[j] <- paste0("fx", j, "~~fx", j)
    fac_int_b[j] <- paste0("fx", j, "~1")
  }
  fac_var_b <- paste(fac_var_b, collapse="; ")
  fac_int_b <- paste(fac_int_b, collapse="; ")
  
  # covariances
  fac_cov_b <- c() # correlations between factors
  count <- 0
  for(j in 1:p){   
    for(m in 1:p){
      if(j != m & m > j){
        count <- count + 1
        fac_cov_b[count] <- paste0("fx", j, "~~", "fx", m)
      }
    } 
  }
  fac_cov_b <- paste(fac_cov_b, collapse = "; ")
  
  model_WF_B <- paste(fac_b, fac_var_b, fac_cov_b, fac_int_b, sep="; ")
  
  model_WF <- paste(model_WF_W, model_WF_B, sep="; ")
  
  ## 2) fit model
  
  fit_WF <- try(sem(model = model_WF, # Barendse & Rossel (2020) use sem() as well (see Appendix A, p.718)
                    data = data_WF),
                silent = TRUE)
  
  if(inherits(fit_WF,"try-error")){
    fit_WF <- NA 
  }
  
  return(list(fit_LF=fit_LF, 
              fit_WF=fit_WF
  ) )
}

# Barendse, M. T., & Rosseel, Y. (2020). Multilevel Modeling in the ‘Wide Format’ Approach with Discrete Data: A Solution for Small Cluster Sizes. Structural Equation Modeling: A Multidisciplinary Journal, 27(5), 696–721. https://doi.org/10.1080/10705511.2019.1689366
