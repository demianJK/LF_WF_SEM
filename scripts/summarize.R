### summarize results of Monte Carlo study

# (install and) load packages
packages <- c("lavaan", 
              "tidyr", 
              "dplyr",
              "xlsx",
              "DescTools"
)
newPackages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(newPackages)){install.packages(newPackages)}
lapply(packages, require, character.only = TRUE)

# set paths and load data
# setwd("C:/Users/Julia/Documents/projects/MSEMcovreg_deprecated")
pathO <- "objects/"
pathF <- "data/fit/"
pathD <- "data/raw/"
source("scripts/functions/projectInfo.R") # infos about MC study
pI <- projectInfo()
source("scripts/functions/covMulti.R") # to compute within- and between-cluster covariance matrix in LF

TableConds <- readRDS(paste0(pathO, "SimConds.rds"))
ncond <- nrow(TableConds)
nrep <- pI$nrep
approaches <- pI$approaches
napproach <- length(approaches)

# if correlations are fixed, they are set here
cB <- TRUE
cW <- TRUE
if (!"cor_B" %in% colnames(TableConds)){
  cor_B <- 0.3
  cB <- FALSE
}
if (!"cor_W" %in% colnames(TableConds)){
  cor_W <- 0.3
  cW <- FALSE
}

# Do you want to center the data in LF and WF?
center <- FALSE

## Data -> ListStats ######################################################################################################################################################
# put all relevant info from single rds files lavaan objects into one list

critsList <- c("Sigma_W_ns", # within-cluster sample covariance matrix
               "Sigma_W_allZero",
               "Sigma_W_def", # 
               "Sigma_W_kappa",
               "Sigma_B_ns", # between-cluster sample covariance matrix
               "Sigma_B_allZero",
               "Sigma_B_def",
               "Sigma_B_kappa",
               # note that for WF, there is only one covariance matrix S_T, whose stats are saved in the Sigma_B slots
               "conv", # Did the model converge?
               "time", # How long did it take till the model converged?
               "W", # within-cluster model parameters
               "B", # between-cluster model parameters
               "coverage_W", # TRUE/FALSE whether confidence interval encompasses population parameter
               "coverage_B",
               "ICC_S", # only for LF (based on lavaans (!) S_pw and S_b)
               "ICC_Sigma" # only for LF (based on Sigma_W and Sigma_B)
               ) 

# create list to save results of each replication
ListStats <- setNames(as.list(approaches), approaches) 
ListStats <- lapply(ListStats, function(x){
  setNames(vector("list", length(critsList)), critsList)
})
for (i in 1:napproach){ # initializing necessary
  for (j in 1:length(critsList)){
    ListStats[[i]][[j]] <- setNames(vector("list", ncond), 1:ncond)
    for (k in 1:ncond){
      ListStats[[i]][[j]][[k]] <- setNames(vector("list", nrep), 1:nrep)
    }
  }
}

errorMessages <- setNames(as.list(approaches), approaches)
for (i in 1:napproach){ # initializing necessary
  errorMessages[[i]]<- setNames(vector("list", ncond), 1:ncond)
  for (k in 1:ncond){
    errorMessages[[i]][[k]] <- setNames(vector("list", nrep), 1:nrep)
  }
}

for (rep in 1:nrep){    
  for (cond in 1:ncond){
    p <- TableConds$p[cond]
    n <- TableConds$n[cond]
    g <- TableConds$g[cond]
    popVar_B <- TableConds$ICC[cond]
    popVar_W <- (1 - TableConds$ICC[cond])
    if (cB){
      popCov_B <- TableConds$cor_B[cond] * popVar_B
    } else{
      popCov_B <- cor_B * popVar_B
    }
    if (cW){
      popCov_W <- TableConds$cor_W[cond] * popVar_W
    } else{
      popCov_W <- cor_W * popVar_W
    }
    # relevant for indexing parameters
    pc <- p * (p + 1) / 2 # number of unique covariance matrix elements (in LF)
    c <- pc-p # numb of unique covariances
    popParams_W <- c(rep(popVar_W, p), rep(popCov_W, c))
    popParams_B <- c(rep(popVar_B, p), rep(popCov_B, c))
    fitList <- readRDS(paste0(pathF, "fit_C", cond, "_R", rep, ".rds"))
    data <- readRDS(paste0(pathD, "data_C", cond, "_R", rep, ".rds"))
    
    # prepare data (reformating and centering)
    data_WF <- pivot_wider(data, names_from = "persons", values_from = 3:(p+2))
    if (center == TRUE){  # grand-mean centering
      data_LF_centered <- as.data.frame(scale(data[,1:p], scale = FALSE))
      data_LF <- cbind(data[,(p+1):(p+2)], data_LF_centered)
      data_WF_centered <- as.data.frame(scale(data_WF[,2:(p*n+1)], scale = FALSE))
      data_WF <- cbind(data_WF[,1], data_WF_centered)
    } else {
      data_LF <- data
    }
      
    
    for (i in 1:napproach){ 
      fit <- fitList[[(paste0("fit_", approaches[i]))]]
      
      ### matrix properties
      # tolerance for zero (1e-08) and definiteness checks taken from matrixcalc functions
      
      if ("LF" %in% approaches[i]){
        covs <- covMulti(data_LF[,3:(3+p-1)], data_LF[,1])
        #lavInspect(fit, "h1")

        # covs$S_pw == fit@SampleStats@YLp[[1]][[2]][["Sigma.W"]]
        eigen_W <- eigen(covs$S_pw)$values
        for (m in 1:p){
          if (abs(eigen_W[m]) < 1e-08){
            eigen_W[m] <- 0
          }
        }
        ListStats[[approaches[i]]][["Sigma_W_ns"]][[cond]][[rep]] <- !any(eigen_W == 0)
        ListStats[[approaches[i]]][["Sigma_W_allZero"]][[cond]][[rep]] <- all(eigen_W == 0)
        if ( !any(eigen_W <= 0) ) { def <- "pd"
        } else if ( !any(eigen_W <  0) ) { def <- "psd" # matrix with all zero eigenvalues is psd
        } else if ( !any(eigen_W >= 0) ) { def <- "nd"
        } else if ( !any(eigen_W >  0) ) { def <- "nsd"
        } else if ( (any(eigen_W >  0) && any(eigen_W < 0)) ) { def <- "id" }
        ListStats[[approaches[i]]][["Sigma_W_def"]][[cond]][[rep]] <- def
        ListStats[[approaches[i]]][["Sigma_W_kappa"]][[cond]][[rep]] <- max(abs(eigen_W))/min(abs(eigen_W))


        Sigma_B <- covs$Sigma_B
        # in lavaan, negative variances, and related covariances are set to zero; we emulate that
        for (k in 1:p){
          if (diag(Sigma_B)[k] < 0){
            Sigma_B[,k] <- 0
            Sigma_B[k,] <- 0
          }
        }
        # Sigma_B == fit@SampleStats@YLp[[1]][[2]][["Sigma.B"]]
        eigen_B <- eigen(Sigma_B)$values
        for (m in 1:p){
          if (abs(eigen_B[m]) < 1e-08){
            eigen_B[m] <- 0
          }
        }
        ListStats[[approaches[i]]][["Sigma_B_ns"]][[cond]][[rep]] <- !any(eigen_B == 0)
        ListStats[[approaches[i]]][["Sigma_B_allZero"]][[cond]][[rep]] <- all(eigen_B == 0)
        if ( !any(eigen_B <= 0) ) { def <- "pd"
        } else if ( !any(eigen_B <  0) ) { def <- "psd"
        } else if ( !any(eigen_B >= 0) ) { def <- "nd"
        } else if ( !any(eigen_B >  0) ) { def <- "nsd"
        } else if ( (any(eigen_B >  0) && any(eigen_B < 0)) ) { def <- "id" }
        ListStats[[approaches[i]]][["Sigma_B_def"]][[cond]][[rep]] <- def
        ListStats[[approaches[i]]][["Sigma_B_kappa"]][[cond]][[rep]] <- max(abs(eigen_B))/min(abs(eigen_B))

      } else { # WF biased

        Sigma_T <- round(cov(data_WF[,2:(p*n+1)]) * (g-1)/g, 3) # here the covariance matrix is rounded in lavaan
        # Sigma_T == lavInspect(fit, "sampstat")$cov
        eigen_T <- eigen(Sigma_T)$values
        for (m in 1:p){
          if (abs(eigen_T[m]) < 1e-08){
            eigen_T[m] <- 0
          }
        }
        ListStats[[approaches[i]]][["Sigma_B_ns"]][[cond]][[rep]] <- !any(eigen_T == 0)
        ListStats[[approaches[i]]][["Sigma_B_allZero"]][[cond]][[rep]] <- all(eigen_T == 0)
        if ( !any(eigen_T <= 0) ) { def <- "pd"
        } else if ( !any(eigen_T <  0) ) { def <- "psd"
        } else if ( !any(eigen_T >= 0) ) { def <- "nd"
        } else if ( !any(eigen_T >  0) ) { def <- "nsd"
        } else if ( (any(eigen_T >  0) && any(eigen_T < 0)) ) { def <- "id" }
        ListStats[[approaches[i]]][["Sigma_B_def"]][[cond]][[rep]] <- def
        ListStats[[approaches[i]]][["Sigma_B_kappa"]][[cond]][[rep]] <- max(abs(eigen_T))/min(abs(eigen_T))

      }

      ### model properties
      
      ListStats[[approaches[i]]][["conv"]][[cond]][[rep]] <- FALSE
      ListStats[[approaches[i]]][["time"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["W"]][[cond]][[rep]] <- NA
      ListStats[[approaches[i]]][["B"]][[cond]][[rep]] <- NA
      
      if(!is.atomic(fit) && fit@Fit@converged){ # if exists and optimizer says converged
        ListStats[[approaches[i]]][["conv"]][[cond]][[rep]] <- fit@Fit@converged
        ListStats[[approaches[i]]][["time"]][[cond]][[rep]] <- unname(fit@timing[["total"]]) # in s
        if ("WFcovreg" %in% approaches[i]){
          ListStats[[approaches[i]]][["time"]][[cond]][[rep]] <- round(unname(fit@timing[["total"]]) + fitList$WFcovreg$covreg_time, 2) # in s
        }
        if ("LF" %in% approaches[i]){
          
          ## parameter estimates
          ListStats[[approaches[i]]][["W"]][[cond]][[rep]] <- fit@Fit@x[1:pc]
          ListStats[[approaches[i]]][["B"]][[cond]][[rep]] <- fit@Fit@x[(pc + 1):(2 * pc)] # the following elements would be the between level intercepts
          
          ## coverage
          logParams_W <- c(rep(TRUE, pc), # within variances and covariances
                           rep(FALSE, p), # within means and (between) factor leadings
                           rep(FALSE, pc), # between variances and covariances
                           rep(FALSE, p)) # between means
          ListStats[[approaches[i]]][["coverage_W"]][[cond]][[rep]] <- popParams_W %()% as.matrix(parameterEstimates(fit)[logParams_W,c("ci.lower", "ci.upper")])
          logParams_B <- c(rep(FALSE, pc), # within variances and covariances
                           rep(FALSE, p), # within means and (between) factor leadings
                           rep(TRUE, pc), # between variances and covariances
                           rep(FALSE, p)) # between means
          ListStats[[approaches[i]]][["coverage_B"]][[cond]][[rep]] <- popParams_B %()% as.matrix(parameterEstimates(fit)[logParams_B, c("ci.lower", "ci.upper")])
          
          # ICC_S 
          # == lavInspect(fit, "h1")$groups$cov /  (lavInspect(fit, "h1")$groups$cov + lavInspect(fit, "h1")$within$cov)
          ListStats[[approaches[i]]][["ICC_S"]][[cond]][[rep]] <- #round(mean(
            lavInspect(fit, "icc")#), 2)
          
          # ICC_Sigma
          # == round(diag(covMulti(data_T[,1:p], data_T[,p+1])$Sigma_B / (covMulti(data_T[,1:p], data_T[,p+1])$Sigma_B + covMulti(data_T[,1:p], data_T[,p+1])$S_pw)), 3)
          ListStats[[approaches[i]]][["ICC_Sigma"]][[cond]][[rep]] <- #round(mean(
            diag(fit@SampleStats@YLp[[1]][[2]][["Sigma.B"]] / (fit@SampleStats@YLp[[1]][[2]][["Sigma.B"]] + fit@SampleStats@YLp[[1]][[2]][["Sigma.W"]]))#), 3)
          
        } else {
          
          ## parameter estimates
          # within: n*p times
          idxW_var <- c() # all equals (n) consecutively
          for (j in 1:p){
            idxW_var <- append(idxW_var, (j-1)*n+1)
          }
          endVarW <- p*n
          idxW_cov <- (endVarW+1):(endVarW+c) # all uniques, then equal (n) again
          idxW_all <- c(idxW_var, idxW_cov) 
          # between: p times
          startB_var <- tail(idxW_all,1)+1+c*(n-1) 
          idxB_var <- startB_var:(startB_var+(p-1))
          startB_cov <- tail(idxB_var,1)+1
          idxB_cov <- startB_cov:(startB_cov+(c-1))
          idxB_all <- c(idxB_var, idxB_cov)
          # the following elements would be the between level (factor) intercepts (means)
          
          ListStats[[approaches[i]]][["W"]][[cond]][[rep]] <- fit@Fit@x[idxW_all]  
          ListStats[[approaches[i]]][["B"]][[cond]][[rep]] <- fit@Fit@x[idxB_all] 
          
          ## coverage
          wvc <- rep(FALSE, p*n+c*n)
          wvc[c(seq(1, p*n, n), # variances (duplicate e.g 111222)
                seq(p*n+1, p*n+c, 1)) # covariances (in row e.g. 121212)
              ] <- TRUE # select only first of duplicate elements (because of equality constraints)
          logParams_W <- c(wvc, # within variances and covariances
                         rep(FALSE, p*n+p*n), # within means and (between) factor leadings
                         rep(FALSE, pc), # between variances and covariances
                         rep(FALSE, p)) # between means
          ListStats[[approaches[i]]][["coverage_W"]][[cond]][[rep]] <- popParams_W %()% as.matrix(parameterEstimates(fit)[logParams_W,c("ci.lower", "ci.upper")])
          
          logParams_B <- c(rep(FALSE, p*n+c*n), # within variances and covariances
                           rep(FALSE, p*n+p*n), # within means and (between) factor leadings
                           rep(TRUE, pc), # between variances and covariances
                           rep(FALSE, p)) # between means
          ListStats[[approaches[i]]][["coverage_B"]][[cond]][[rep]] <- popParams_B %()% as.matrix(parameterEstimates(fit)[logParams_B, c("ci.lower", "ci.upper")])
        }
      } else if (!is.null(attr(fit, "condition"))) { # save error note (mainly not positive definite S in WF?)
        errorMessages[[approaches[i]]][[cond]][[rep]] <- attr(fit, "condition")[[1]] # only error, not function call
      }
    }
    gc()
  }
  print(rep)
}

saveRDS(ListStats, file = paste0(pathO, "ListStats.rds"))



### ListStats --> TableStats  ######################################################################################################################################################
# summarize info from list into table

matrx <- c("ns", "allZero", "pd", "psd", "nd", "nsd", "id", "kappa", "Infkappa")
estim <- c("MAE", "MSE", "RMSE", "Bias", "Var")
critsTable <- c(paste(matrx,  "W", sep="_"),
                paste(matrx,  "B", sep="_"),
  
                "conv", "time", "timesd",
                
                ## within
                # variance
                paste(estim, "var", "W", sep="_"),
                paste(paste0("rel", estim), "var", "W", sep="_"),
                # covariance
                paste(estim, "cov", "W", sep="_"),
                paste(paste0("rel", estim), "cov", "W", sep="_"),
                # overall
                paste(estim, "W", sep="_"),
                paste(paste0("rel", estim), "W", sep="_"),
                # negative variances
                "negVar_W",
                
                ## between
                # variance
                paste(estim, "var", "B", sep="_"),
                paste(paste0("rel", estim), "var", "B", sep="_"),
                # covariance
                paste(estim, "cov", "B", sep="_"),
                paste(paste0("rel", estim), "cov", "B", sep="_"),
                # overall
                paste(estim, "B", sep="_"),
                paste(paste0("rel", estim), "B", sep="_"),
                
                # negative variances
                "negVar_B",
                
                # coverage
                "coverage_W",
                "coverage_B",
                
                # Intraclass Correlation
                "ICC_S",
                "ICC_S_neg", # percentage of negative ICCs (which are set to zero for the mean ICC)
                "ICC_Sigma",
                "ICC_Sigma_neg",
                "ICC_Sigma_theta",
                "ICC_Sigma_theta_neg"
                )

colTable <- c()
for (i in 1:napproach){
  colTable <- append(colTable, paste(critsTable, approaches[i], sep="_"))
}
tmp <- as.data.frame(matrix(NA, ncol=length(colTable), nrow=ncond))
colnames(tmp) <- colTable
TableStats <- TableConds
TableStats$df_LF <- 0
TableStats <- dplyr::mutate(TableStats, df_WF= ((p*n)*(p*n+1)/2 + (p*n))  - ((p)*(p+1) + p), .keep="all") # (p)*(p+1) weil within und between vars+covs; p weil nur between means geschätzt
TableStats <- cbind(TableStats, tmp)

# note that the relative parameters (e.g., relRMSE) are not computed in the order according to the formula, 
# but the order is changed on mathematical basis to allow for more efficient code

for (cond in 1:ncond){
  
  # required for stat props
  p <- TableConds$p[cond]
  n <- TableConds$n[cond]
  popVar_B <- TableConds$ICC[cond]
  popVar_W <- (1 - TableConds$ICC[cond])
  popVar_B <- TableConds$ICC[cond]
  popVar_W <- (1 - TableConds$ICC[cond])
  if (cB){
    popCov_B <- TableConds$cor_B[cond] * popVar_B
  } else{
    popCov_B <- cor_B * popVar_B
  }
  if (cW){
    popCov_W <- TableConds$cor_W[cond] * popVar_W
  } else{
    popCov_W <- cor_W * popVar_W
  }
  pc <- p * (p + 1) / 2
  c <- pc - p
  # sorting of params in LF and WF: varW (p), covW (c), varB (p), covB (c)
  
  for (i in 1:napproach){

    # matrix properties
    if ("LF" %in% approaches[i]){
      
      tmpList <- unname(unlist(ListStats[[approaches[i]]][["Sigma_W_ns"]][[cond]]))
      tmpList <- factor(tmpList, levels=c("TRUE", "FALSE"))
      tmp <- table(tmpList)
      TableStats[cond, paste("ns_W", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
      
      tmpList <- unname(unlist(ListStats[[approaches[i]]][["Sigma_W_allZero"]][[cond]]))
      tmpList <- factor(tmpList, levels=c("TRUE", "FALSE"))
      tmp <- table(tmpList)
      TableStats[cond, paste("allZero_W", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
      
      tmpList <- unname(unlist(ListStats[[approaches[i]]][["Sigma_W_def"]][[cond]]))
      tmpList <- factor(tmpList, levels=c("pd", "psd", "nd", "nsd", "id"))
      tmp <- table(tmpList)
      TableStats[cond, paste("pd_W", approaches[i], sep="_")] <- tmp[1] / length(tmpList) * 100
      TableStats[cond, paste("psd_W", approaches[i], sep="_")] <- (tmp[1] + tmp[2]) / length(tmpList) * 100
      TableStats[cond, paste("nd_W", approaches[i], sep="_")] <- tmp[3] / length(tmpList) * 100
      TableStats[cond, paste("nsd_W", approaches[i], sep="_")] <- (tmp[3] + tmp[4]) / length(tmpList) * 100
      TableStats[cond, paste("id_W", approaches[i], sep="_")] <- tmp[5] / length(tmpList) * 100
      
      tmpList <- unname(unlist(ListStats[[approaches[i]]][["Sigma_W_kappa"]][[cond]]))
      tmpInf <- is.infinite(tmpList)
      if (any(tmpInf)){
        tmpList[tmpInf] <- NA
      }
      TableStats[cond, paste("kappa_W", approaches[i], sep="_")] <- round(mean(tmpList, na.rm=TRUE), 3)
      tmpInf <- factor(tmpInf, levels=c("TRUE", "FALSE"))
      tmp <- table(tmpInf)
      TableStats[cond, paste("Infkappa_W", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
      
      # ICC_S
      ICC_S <- unname(unlist(ListStats[[approaches[i]]][["ICC_S"]][[cond]]))
      tmp <- sum(ICC_S < 0)
      if (tmp > 0){
        TableStats[cond, paste("ICC_S_neg", approaches[i], sep="_")] <- round(tmp/length(ICC_S) * 100, 2)
        ICC_S[ICC_S < 0] <- 0
      }
      TableStats[cond, paste("ICC_S", approaches[i], sep="_")] <- round(mean(ICC_S, na.rm=TRUE), 2)
      
      # ICC_Sigma
      ICC_Sigma <- unname(unlist(ListStats[[approaches[i]]][["ICC_Sigma"]][[cond]]))
      tmp <- sum(ICC_Sigma < 0)
      if (tmp > 0){
        TableStats[cond, paste("ICC_Sigma_neg", approaches[i], sep="_")] <- round(tmp/length(ICC_Sigma) * 100, 2)
        ICC_Sigma[ICC_Sigma < 0] <- 0
      }
      TableStats[cond, paste("ICC_Sigma", approaches[i], sep="_")] <- round(mean(ICC_Sigma, na.rm=TRUE), 2)
    
    } 
      
      tmpList <- unname(unlist(ListStats[[approaches[i]]][["Sigma_B_ns"]][[cond]]))
      tmpList <- factor(tmpList, levels=c("TRUE", "FALSE"))
      tmp <- table(tmpList)
      TableStats[cond, paste("ns_B", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
      
      tmpList <- unname(unlist(ListStats[[approaches[i]]][["Sigma_B_allZero"]][[cond]]))
      tmpList <- factor(tmpList, levels=c("TRUE", "FALSE"))
      tmp <- table(tmpList)
      TableStats[cond, paste("allZero_B", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
      
      tmpList <- unname(unlist(ListStats[[approaches[i]]][["Sigma_B_def"]][[cond]]))
      tmpList <- factor(tmpList, levels=c("pd", "psd", "nd", "nsd", "id"))
      tmp <- table(tmpList)
      TableStats[cond, paste("pd_B", approaches[i], sep = "_")] <- tmp[1] / length(tmpList) * 100
      TableStats[cond, paste("psd_B", approaches[i], sep = "_")] <- (tmp[1] + tmp[2]) / length(tmpList) * 100
      TableStats[cond, paste("nd_B", approaches[i], sep = "_")] <- tmp[3] / length(tmpList) * 100
      TableStats[cond, paste("nsd_B", approaches[i], sep = "_")] <- (tmp[3] + tmp[4]) / length(tmpList) * 100
      TableStats[cond, paste("id_B", approaches[i], sep = "_")] <- tmp[5] / length(tmpList) * 100
      
      tmpList <-
        unname(unlist(ListStats[[approaches[i]]][["Sigma_B_kappa"]][[cond]]))
      tmpInf <- is.infinite(tmpList)
      if (any(tmpInf)) {
        tmpList[tmpInf] <- NA
      }
      TableStats[cond, paste("kappa_B", approaches[i], sep = "_")] <-
        round(mean(tmpList, na.rm = TRUE), 3)
      tmpInf <- factor(tmpInf, levels = c("TRUE", "FALSE"))
      tmp <- table(tmpInf)
      TableStats[cond, paste("Infkappa_B", approaches[i], sep = "_")] <-
        tmp[1] / (tmp[1] + tmp[2]) * 100
      

    # model properties
    
    # convergence
    tmpList <- unname(unlist(ListStats[[approaches[i]]][["conv"]][[cond]]))
    tmpList <- factor(tmpList, levels=c("TRUE", "FALSE"))
    tmp <- table(tmpList)
    TableStats[cond, paste("conv", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
    
    if (TableStats[cond, paste("conv", approaches[i], sep="_")] > 0){
      
      # computation time
      tmp <- unname(unlist(ListStats[[approaches[i]]][["time"]][[cond]])) 
      TableStats[cond, paste("time", approaches[i], sep="_")] <- round(mean(tmp, na.rm=TRUE), 2) 
      TableStats[cond, paste("timesd", approaches[i], sep="_")] <- round(sd(tmp, na.rm=TRUE), 2) 
      
      ## Parameter Estimates
      
      ## within
      W <- na.omit(as.data.frame(do.call(rbind, ListStats[[approaches[i]]][["W"]][[cond]])))
      
      ## variance
      var_W <- W[,1:p]
      
      # negative variances 
      tv <- c()
      for (rep in 1:nrep){
        tv[rep] <- any(var_W[rep,] < 0) 
      }
      tv <- factor(tv, levels=c("TRUE", "FALSE"))
      tmp <- table(tv)
      TableStats[cond, paste("negVar_W", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
        
      MAE_var_W <- mean(colMeans(abs(var_W - popVar_W))) 
      TableStats[cond, paste("MAE_var_W", approaches[i], sep="_")] <- round(MAE_var_W, 2) 
      TableStats[cond, paste("relMAE_var_W", approaches[i], sep="_")] <- round( ((MAE_var_W / popVar_W)* 100), 2) 
        
      MSE_var_W <- mean(colMeans((var_W - popVar_W)^2)) 
      TableStats[cond, paste("MSE_var_W", approaches[i], sep="_")] <- round(MSE_var_W, 2) 
      TableStats[cond, paste("relMSE_var_W", approaches[i], sep="_")] <- round( ((MSE_var_W / popVar_W)* 100), 2)
        
      RMSE_var_W <- mean(sqrt(colMeans((var_W - popVar_W)^2)))
      TableStats[cond, paste("RMSE_var_W", approaches[i], sep="_")] <- round(RMSE_var_W, 2) 
      TableStats[cond, paste("relRMSE_var_W", approaches[i], sep="_")] <- round( ((RMSE_var_W / popVar_W)* 100), 2) 
        
      Bias_var_W <- mean(colMeans((var_W - popVar_W)))
      TableStats[cond, paste("Bias_var_W", approaches[i], sep="_")] <- round(Bias_var_W, 2)
      TableStats[cond, paste("relBias_var_W", approaches[i], sep="_")] <- round( ((Bias_var_W / popVar_W)* 100), 2)
      
      Var_var_W <- mean(colMeans( ((var_W - mean(colMeans((var_W))))^2) ))
      TableStats[cond, paste("Var_var_W", approaches[i], sep="_")] <- round(Var_var_W, 2)
      TableStats[cond, paste("relVar_var_W", approaches[i], sep="_")] <- round( ((Var_var_W / popVar_W)* 100), 2)
        
      ## covariance
      cov_W <- as.data.frame(W[,(p+1):(p+c)])
        
      MAE_cov_W <- mean(colMeans(abs(cov_W - popCov_W))) 
      TableStats[cond, paste("MAE_cov_W", approaches[i], sep="_")] <- round(MAE_cov_W, 2) 
      TableStats[cond, paste("relMAE_cov_W", approaches[i], sep="_")] <- round( ((MAE_cov_W / popCov_W)* 100), 2)
        
      MSE_cov_W <- mean(colMeans((cov_W - popCov_W)^2)) 
      TableStats[cond, paste("MSE_cov_W", approaches[i], sep="_")] <- round(MSE_cov_W, 2) 
      TableStats[cond, paste("relMSE_cov_W", approaches[i], sep="_")] <- round( ((MSE_cov_W / popCov_W)* 100), 2)
        
      RMSE_cov_W <- mean(sqrt(colMeans((cov_W - popCov_W)^2)))
      TableStats[cond, paste("RMSE_cov_W", approaches[i], sep="_")] <- round(RMSE_cov_W, 2) 
      TableStats[cond, paste("relRMSE_cov_W", approaches[i], sep="_")] <- round( ((RMSE_cov_W / popCov_W)* 100), 2)
        
      Bias_cov_W <- mean(colMeans((cov_W - popCov_W)))
      TableStats[cond, paste("Bias_cov_W", approaches[i], sep="_")] <- round(Bias_cov_W, 2)
      TableStats[cond, paste("relBias_cov_W", approaches[i], sep="_")] <- round( ((Bias_cov_W / popCov_W)* 100), 2)
        
      Var_cov_W <- mean(colMeans( ((cov_W - mean(colMeans((cov_W))))^2) ))
      TableStats[cond, paste("Var_cov_W", approaches[i], sep="_")] <- round(Var_cov_W, 2)
      TableStats[cond, paste("relVar_cov_W", approaches[i], sep="_")] <- round( ((Var_cov_W / popCov_W)* 100), 2)
        
      ## overall (variance and covariance)
        
      TableStats[cond, paste("MAE_W", approaches[i], sep="_")] <- round( mean( c(MAE_var_W, MAE_cov_W) ), 2)
      TableStats[cond, paste("MSE_W", approaches[i], sep="_")] <- round( mean( c(MSE_var_W, MSE_cov_W) ), 2)
      TableStats[cond, paste("RMSE_W", approaches[i], sep="_")] <- round( mean( c(RMSE_var_W, RMSE_cov_W) ), 2)
      TableStats[cond, paste("Bias_W", approaches[i], sep="_")] <- round( mean( c(Bias_var_W, Bias_cov_W) ), 2)
      TableStats[cond, paste("Var_W", approaches[i], sep="_")] <- round( mean( c(Var_var_W, Var_cov_W) ), 2)
        
      TableStats[cond, paste("relMAE_W", approaches[i], sep="_")] <- round( mean( c(((MAE_var_W / popVar_W)* 100), ((MAE_cov_W / popCov_W)* 100)) ), 2)
      TableStats[cond, paste("relMSE_W", approaches[i], sep="_")] <- round( mean( c(((MSE_var_W / popVar_W)* 100), ((MSE_cov_W / popCov_W)* 100)) ), 2)
      TableStats[cond, paste("relRMSE_W", approaches[i], sep="_")] <- round( mean( c(((RMSE_var_W / popVar_W)* 100), ((RMSE_cov_W / popCov_W)* 100)) ), 2)
      TableStats[cond, paste("relBias_W", approaches[i], sep="_")] <- round( mean( c(((Bias_var_W / popVar_W)* 100), ((Bias_cov_W / popCov_W)* 100)) ), 2)
      TableStats[cond, paste("relVar_W", approaches[i], sep="_")] <- round( mean( c(((Var_var_W / popVar_W)* 100), ((MSE_cov_W / popCov_W)* 100)) ), 2)
        
      
      ### between
      B <- na.omit(as.data.frame(do.call(rbind, ListStats[[approaches[i]]][["B"]][[cond]])))
      
      # variances
      var_B <- B[,1:p]
      
      # negative variances 
      tv <- c()
      for (rep in 1:nrep){
        tv[rep] <- any(var_B[rep,] < 0) 
      }
      tv <- factor(tv, levels=c("TRUE", "FALSE"))
      tmp <- table(tv)
      TableStats[cond, paste("negVar_B", approaches[i], sep="_")] <- tmp[1] / ( tmp[1] + tmp[2]) * 100
      
      MAE_var_B <- mean(colMeans(abs(var_B - popVar_B))) 
      TableStats[cond, paste("MAE_var_B", approaches[i], sep="_")] <- round(MAE_var_B, 2) 
      TableStats[cond, paste("relMAE_var_B", approaches[i], sep="_")] <- round( ((MAE_var_B / popVar_B)* 100), 2) 
      
      MSE_var_B <- mean(colMeans((var_B - popVar_B)^2)) 
      TableStats[cond, paste("MSE_var_B", approaches[i], sep="_")] <- round(MSE_var_B, 2) 
      TableStats[cond, paste("relMSE_var_B", approaches[i], sep="_")] <- round( ((MSE_var_B / popVar_B)* 100), 2)
      
      RMSE_var_B <- mean(sqrt(colMeans((var_B - popVar_B)^2)))
      TableStats[cond, paste("RMSE_var_B", approaches[i], sep="_")] <- round(RMSE_var_B, 2) 
      TableStats[cond, paste("relRMSE_var_B", approaches[i], sep="_")] <- round( ((RMSE_var_B / popVar_B)* 100), 2) 
      
      Bias_var_B <- mean(colMeans((var_B - popVar_B)))
      TableStats[cond, paste("Bias_var_B", approaches[i], sep="_")] <- round(Bias_var_B, 2)
      TableStats[cond, paste("relBias_var_B", approaches[i], sep="_")] <- round( ((Bias_var_B / popVar_B)* 100), 2)
      
      Var_var_B <- mean(colMeans( ((var_B - mean(colMeans((var_B))))^2) ))
      TableStats[cond, paste("Var_var_B", approaches[i], sep="_")] <- round(Var_var_B, 2)
      TableStats[cond, paste("relVar_var_B", approaches[i], sep="_")] <- round( ((Var_var_B / popVar_B)* 100), 2)
      
      ## covariance
      cov_B <- as.data.frame(B[,(p+1):(p+c)])
      
      MAE_cov_B <- mean(colMeans(abs(cov_B - popCov_B))) 
      TableStats[cond, paste("MAE_cov_B", approaches[i], sep="_")] <- round(MAE_cov_B, 2) 
      TableStats[cond, paste("relMAE_cov_B", approaches[i], sep="_")] <- round( ((MAE_cov_B / popCov_B)* 100), 2)
      
      MSE_cov_B <- mean(colMeans((cov_B - popCov_B)^2)) 
      TableStats[cond, paste("MSE_cov_B", approaches[i], sep="_")] <- round(MSE_cov_B, 2) 
      TableStats[cond, paste("relMSE_cov_B", approaches[i], sep="_")] <- round( ((MSE_cov_B / popCov_B)* 100), 2)
      
      RMSE_cov_B <- mean(sqrt(colMeans((cov_B - popCov_B)^2)))
      TableStats[cond, paste("RMSE_cov_B", approaches[i], sep="_")] <- round(RMSE_cov_B, 2) 
      TableStats[cond, paste("relRMSE_cov_B", approaches[i], sep="_")] <- round( ((RMSE_cov_B / popCov_B)* 100), 2)
      
      Bias_cov_B <- mean(colMeans((cov_B - popCov_B)))
      TableStats[cond, paste("Bias_cov_B", approaches[i], sep="_")] <- round(Bias_cov_B, 2)
      TableStats[cond, paste("relBias_cov_B", approaches[i], sep="_")] <- round( ((Bias_cov_B / popCov_B)* 100), 2)
      
      Var_cov_B <- mean(colMeans( ((cov_B - mean(colMeans((cov_B))))^2) ))
      TableStats[cond, paste("Var_cov_B", approaches[i], sep="_")] <- round(Var_cov_B, 2)
      TableStats[cond, paste("relVar_cov_B", approaches[i], sep="_")] <- round( ((Var_cov_B / popCov_B)* 100), 2)
      
      ## overall (variance and covariance)
      
      TableStats[cond, paste("MAE_B", approaches[i], sep="_")] <- round( mean( c(MAE_var_B, MAE_cov_B) ), 2)
      TableStats[cond, paste("MSE_B", approaches[i], sep="_")] <- round( mean( c(MSE_var_B, MSE_cov_B) ), 2)
      TableStats[cond, paste("RMSE_B", approaches[i], sep="_")] <- round( mean( c(RMSE_var_B, RMSE_cov_B) ), 2)
      TableStats[cond, paste("Bias_B", approaches[i], sep="_")] <- round( mean( c(Bias_var_B, Bias_cov_B) ), 2)
      TableStats[cond, paste("Var_B", approaches[i], sep="_")] <- round( mean( c(Var_var_B, Var_cov_B) ), 2)
      
      TableStats[cond, paste("relMAE_B", approaches[i], sep="_")] <- round( mean( c(((MAE_var_B / popVar_B)* 100), ((MAE_cov_B / popCov_B)* 100)) ), 2)
      TableStats[cond, paste("relMSE_B", approaches[i], sep="_")] <- round( mean( c(((MSE_var_B / popVar_B)* 100), ((MSE_cov_B / popCov_B)* 100)) ), 2)
      TableStats[cond, paste("relRMSE_B", approaches[i], sep="_")] <- round( mean( c(((RMSE_var_B / popVar_B)* 100), ((RMSE_cov_B / popCov_B)* 100)) ), 2)
      TableStats[cond, paste("relBias_B", approaches[i], sep="_")] <- round( mean( c(((Bias_var_B / popVar_B)* 100), ((Bias_cov_B / popCov_B)* 100)) ), 2)
      TableStats[cond, paste("relVar_B", approaches[i], sep="_")] <- round( mean( c(((Var_var_B / popVar_B)* 100), ((MSE_cov_B / popCov_B)* 100)) ), 2)
       
      # ICC_Sigma_theta (based on model parameters of unstructured model)
      ICC_Sigma_theta <- unlist((var_B/(var_W+var_B)))
      tmp <- sum(ICC_Sigma_theta < 0)
      if (tmp > 0){
        TableStats[cond, paste("ICC_Sigma_theta_neg", approaches[i], sep="_")] <- round(tmp/length(ICC_Sigma_theta) * 100, 2)
        ICC_Sigma_theta[ICC_Sigma_theta < 0] <- 0
      }
      TableStats[cond, paste("ICC_Sigma_theta", approaches[i], sep="_")] <- round(mean(ICC_Sigma_theta, na.rm=TRUE), 2)
      
      # coverage
      coverage_W <- unname(unlist(ListStats[[approaches[i]]][["coverage_W"]][[cond]]))
      TableStats[cond, paste("coverage_W", approaches[i], sep="_")] <- round(sum(coverage_W)/length(coverage_W) * 100, 2)
      coverage_B <- unname(unlist(ListStats[[approaches[i]]][["coverage_B"]][[cond]]))
      TableStats[cond, paste("coverage_B", approaches[i], sep="_")] <- round(sum(coverage_B)/length(coverage_B) * 100, 2)
      
    }  
  }
}

# aggregate over within and between

for (i in 1:napproach){
  TableStats[paste("MAE", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("MAE_W", "MAE_B"), approaches[i], sep="_"))), 2)
  TableStats[paste("MSE", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("MSE_W", "MSE_B"), approaches[i], sep="_"))), 2)
  TableStats[paste("RMSE", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("RMSE_W", "RMSE_B"), approaches[i], sep="_"))), 2)
  TableStats[paste("Bias", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("Bias_W", "Bias_B"), approaches[i], sep="_"))), 2)
  TableStats[paste("Var", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("Var_W", "Var_B"), approaches[i], sep="_"))), 2)
  
  TableStats[paste("relMAE", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("relMAE_W", "relMAE_B"), approaches[i], sep="_"))), 2)
  TableStats[paste("relMSE", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("relMSE_W", "relMSE_B"), approaches[i], sep="_"))), 2)
  TableStats[paste("relRMSE", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("relRMSE_W", "relRMSE_B"), approaches[i], sep="_"))), 2)
  TableStats[paste("relBias", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("relBias_W", "relBias_B"), approaches[i], sep="_"))), 2)
  TableStats[paste("relVar", approaches[i], sep="_")] <- round(rowMeans(select(TableStats, paste(c("relVar_W", "relVar_B"), approaches[i], sep="_"))), 2)
}

TableStats <- dplyr::mutate(TableStats, S_LF= ((p)*(p+1)/2 + (p)), .keep="all") # unique elements of S_LF
TableStats <- dplyr::mutate(TableStats, S_WF= ((p*n)*(p*n+1)/2 + (p*n)), .keep="all") # unique elements of S_WF
TableStats <- dplyr::mutate(TableStats, ratio_LF= (p)/g, .keep="all") # cols:rows LF-B
TableStats <- dplyr::mutate(TableStats, ratio_WF= (p*n)/g, .keep="all") # cols:rows WF-T

write.xlsx( TableStats, file = paste0(pathO, "TableStats.xlsx"))



### TableStats --> dat  ######################################################################################################################################################
# rearrange sumamrized table into data frame that can be used for figures

df <- c()
for (i in 1:napproach){
  if (grepl("LF", approaches[i])){
    df <- append(df, TableStats$df_LF)
  } else if (grepl("WF", approaches[i])){
    df <- append(df, TableStats$df_WF)
  }
}

ratio <- c()
for (i in 1:length(approaches)){
  if (grepl("LF", approaches[i])){
    ratio <- append(ratio, TableStats$ratio_LF)
  } else if (grepl("WF", approaches[i])){
    ratio <- append(ratio, TableStats$ratio_WF)
  }
}

dat <- data.frame(cond = rep(1:ncond, napproach),
                  n = rep(TableStats$n, napproach),
                  g = rep(TableStats$g, napproach),
                  p = rep(TableStats$p, napproach),
                  ratio  = ratio,
                  ICC = rep(TableStats$ICC, napproach),
                  # cor_W = rep(TableStats$cor_W, napproach),
                  # cor_B = rep(TableStats$cor_B, napproach),
                  approach = factor(rep(approaches, each=ncond), level=approaches, ordered=TRUE),
                  df = df,
                  # matrix properties
                  ns_B = unlist(TableStats[, grepl("ns_B", colnames(TableStats))], use.names=FALSE),
                  pd_B = unlist(TableStats[, grepl("pd_B", colnames(TableStats))], use.names=FALSE),
                  psd_B = unlist(TableStats[, grepl("psd_B", colnames(TableStats))], use.names=FALSE),
                  nd_B = unlist(TableStats[, grepl("nd_B", colnames(TableStats))], use.names=FALSE),
                  nsd_B = unlist(TableStats[, grepl("nsd_B", colnames(TableStats))], use.names=FALSE),
                  id_B = unlist(TableStats[, grepl("id_B", colnames(TableStats))], use.names=FALSE),
                  kappa_B = unlist(TableStats[, grepl("^kappa_B", colnames(TableStats))], use.names=FALSE),
                  Infkappa_B = unlist(TableStats[, grepl("Infkappa_B", colnames(TableStats))], use.names=FALSE),
                  # model properties
                  conv = unlist(TableStats[, grepl("conv", colnames(TableStats))], use.names=FALSE),
                  time = c(unlist(TableStats[, grepl("time_", colnames(TableStats))], use.names=FALSE)),
                  ## type
                  # variance
                  MAE_var_B = unlist(TableStats[, paste0("MAE_var_B_", approaches)], use.names=FALSE),
                  MAE_var_W = unlist(TableStats[, paste0("MAE_var_W_", approaches)], use.names=FALSE),
                  MSE_var_B = unlist(TableStats[, paste0("MSE_var_B_", approaches)], use.names=FALSE),
                  MSE_var_W = unlist(TableStats[, paste0("MSE_var_W_", approaches)], use.names=FALSE),
                  RMSE_var_B = unlist(TableStats[, paste0("RMSE_var_B_", approaches)], use.names=FALSE),
                  RMSE_var_W = unlist(TableStats[, paste0("RMSE_var_W_", approaches)], use.names=FALSE),
                  Bias_var_B = unlist(TableStats[, paste0("Bias_var_B_", approaches)], use.names=FALSE),
                  Bias_var_W = unlist(TableStats[, paste0("Bias_var_W_", approaches)], use.names=FALSE),
                  Var_var_B = unlist(TableStats[, paste0("Var_var_B_", approaches)], use.names=FALSE),
                  Var_var_W = unlist(TableStats[, paste0("Var_var_W_", approaches)], use.names=FALSE),
                  relMAE_var_B = unlist(TableStats[, paste0("relMAE_var_B_", approaches)], use.names=FALSE),
                  relMAE_var_W = unlist(TableStats[, paste0("relMAE_var_W_", approaches)], use.names=FALSE),
                  relMSE_var_B = unlist(TableStats[, paste0("relMSE_var_B_", approaches)], use.names=FALSE),
                  relMSE_var_W = unlist(TableStats[, paste0("relMSE_var_W_", approaches)], use.names=FALSE),
                  relRMSE_var_B = unlist(TableStats[, paste0("relRMSE_var_B_", approaches)], use.names=FALSE),
                  relRMSE_var_W = unlist(TableStats[, paste0("relRMSE_var_W_", approaches)], use.names=FALSE),
                  relBias_var_B = unlist(TableStats[, paste0("relBias_var_B_", approaches)], use.names=FALSE),
                  relBias_var_W = unlist(TableStats[, paste0("relBias_var_W_", approaches)], use.names=FALSE),
                  relVar_var_B = unlist(TableStats[, paste0("relVar_var_B_", approaches)], use.names=FALSE),
                  relVar_var_W = unlist(TableStats[, paste0("relVar_var_W_", approaches)], use.names=FALSE),
                  negVar_B = unlist(TableStats[, paste0("negVar_B_", approaches)], use.names=FALSE),
                  negVar_W = unlist(TableStats[, paste0("negVar_W_", approaches)], use.names=FALSE),
                  # covariance
                  MAE_cov_B = unlist(TableStats[, paste0("MAE_cov_B_", approaches)], use.names=FALSE),
                  MAE_cov_W = unlist(TableStats[, paste0("MAE_cov_W_", approaches)], use.names=FALSE),
                  MSE_cov_B = unlist(TableStats[, paste0("MSE_cov_B_", approaches)], use.names=FALSE),
                  MSE_cov_W = unlist(TableStats[, paste0("MSE_cov_W_", approaches)], use.names=FALSE),
                  RMSE_cov_B = unlist(TableStats[, paste0("RMSE_cov_B_", approaches)], use.names=FALSE),
                  RMSE_cov_W = unlist(TableStats[, paste0("RMSE_cov_W_", approaches)], use.names=FALSE),
                  Bias_cov_B = unlist(TableStats[, paste0("Bias_cov_B_", approaches)], use.names=FALSE),
                  Bias_cov_W = unlist(TableStats[, paste0("Bias_cov_W_", approaches)], use.names=FALSE),
                  Var_cov_B = unlist(TableStats[, paste0("Var_cov_B_", approaches)], use.names=FALSE),
                  Var_cov_W = unlist(TableStats[, paste0("Var_cov_W_", approaches)], use.names=FALSE),
                  relMAE_cov_B = unlist(TableStats[, paste0("relMAE_cov_B_", approaches)], use.names=FALSE),
                  relMAE_cov_W = unlist(TableStats[, paste0("relMAE_cov_W_", approaches)], use.names=FALSE),
                  relMSE_cov_B = unlist(TableStats[, paste0("relMSE_cov_B_", approaches)], use.names=FALSE),
                  relMSE_cov_W = unlist(TableStats[, paste0("relMSE_cov_W_", approaches)], use.names=FALSE),
                  relRMSE_cov_B = unlist(TableStats[, paste0("relRMSE_cov_B_", approaches)], use.names=FALSE),
                  relRMSE_cov_W = unlist(TableStats[, paste0("relRMSE_cov_W_", approaches)], use.names=FALSE),
                  relBias_cov_B = unlist(TableStats[, paste0("relBias_cov_B_", approaches)], use.names=FALSE),
                  relBias_cov_W = unlist(TableStats[, paste0("relBias_cov_W_", approaches)], use.names=FALSE),
                  relVar_cov_B = unlist(TableStats[, paste0("relVar_cov_B_", approaches)], use.names=FALSE),
                  relVar_cov_W = unlist(TableStats[, paste0("relVar_cov_W_", approaches)], use.names=FALSE),
                  # level
                  MAE_B = unlist(TableStats[, paste0("MAE_B_", approaches)], use.names=FALSE),
                  MAE_W = unlist(TableStats[, paste0("MAE_W_", approaches)], use.names=FALSE),
                  MSE_B = unlist(TableStats[, paste0("MSE_B_", approaches)], use.names=FALSE),
                  MSE_W = unlist(TableStats[, paste0("MSE_W_", approaches)], use.names=FALSE),
                  RMSE_B = unlist(TableStats[, paste0("RMSE_B_", approaches)], use.names=FALSE),
                  RMSE_W = unlist(TableStats[, paste0("RMSE_W_", approaches)], use.names=FALSE),
                  Bias_B = unlist(TableStats[, paste0("Bias_B_", approaches)], use.names=FALSE),
                  Bias_W = unlist(TableStats[, paste0("Bias_W_", approaches)], use.names=FALSE),
                  Var_B = unlist(TableStats[, paste0("Var_B_", approaches)], use.names=FALSE),
                  Var_W = unlist(TableStats[, paste0("Var_W_", approaches)], use.names=FALSE),
                  relMAE_B = unlist(TableStats[, paste0("relMAE_B_", approaches)], use.names=FALSE),
                  relMAE_W = unlist(TableStats[, paste0("relMAE_W_", approaches)], use.names=FALSE),
                  relMSE_B = unlist(TableStats[, paste0("relMSE_B_", approaches)], use.names=FALSE),
                  relMSE_W = unlist(TableStats[, paste0("relMSE_W_", approaches)], use.names=FALSE),
                  relRMSE_B = unlist(TableStats[, paste0("relRMSE_B_", approaches)], use.names=FALSE),
                  relRMSE_W = unlist(TableStats[, paste0("relRMSE_W_", approaches)], use.names=FALSE),
                  relBias_B = unlist(TableStats[, paste0("relBias_B_", approaches)], use.names=FALSE),
                  relBias_W = unlist(TableStats[, paste0("relBias_W_", approaches)], use.names=FALSE),
                  relVar_B = unlist(TableStats[, paste0("relVar_B_", approaches)], use.names=FALSE),
                  relVar_W = unlist(TableStats[, paste0("relVar_W_", approaches)], use.names=FALSE),
                  # aggregated
                  MAE = unlist(TableStats[, paste0("MAE_", approaches)], use.names=FALSE),
                  MSE = unlist(TableStats[, paste0("MSE_", approaches)], use.names=FALSE),
                  RMSE = unlist(TableStats[, paste0("RMSE_", approaches)], use.names=FALSE),
                  Bias = unlist(TableStats[, paste0("Bias_", approaches)], use.names=FALSE),
                  Var = unlist(TableStats[, paste0("Var_", approaches)], use.names=FALSE),
                  relMAE = unlist(TableStats[, paste0("relMAE_", approaches)], use.names=FALSE),
                  relMSE = unlist(TableStats[, paste0("relMSE_", approaches)], use.names=FALSE),
                  relRMSE = unlist(TableStats[, paste0("relRMSE_", approaches)], use.names=FALSE),
                  relBias = unlist(TableStats[, paste0("relBias_", approaches)], use.names=FALSE),
                  relVar = unlist(TableStats[, paste0("relVar_", approaches)], use.names=FALSE),
                  # coverage
                  coverage_W = unlist(TableStats[, paste0("coverage_W_", approaches)], use.names=FALSE),
                  coverage_B = unlist(TableStats[, paste0("coverage_B_", approaches)], use.names=FALSE),
                  # ICC
                  ICC_S = unlist(TableStats[, paste0("ICC_S_", approaches)], use.names=FALSE),
                  ICC_Sigma = unlist(TableStats[, paste0("ICC_Sigma_", approaches)], use.names=FALSE),
                  ICC_Sigma_theta = unlist(TableStats[, paste0("ICC_Sigma_theta_", approaches)], use.names=FALSE)
)


saveRDS(dat, paste0(pathO, "dat.rds"))
