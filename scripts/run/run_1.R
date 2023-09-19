# multiple scripts in individual sessions are run for using multiple cores
# copy "run_1.R" times the number of cores you want to use and change the number in the name and in the respective script (it's used for the random seed below)
script <- 1

# set paths
pathS <- "scripts/functions/" # path for local functions 
pathD <- "data/raw/" # path for data 
pathF <- "data/fit/" # path for fit 

# source local functions
source(paste0(pathS, "simulateData.R"))
source(paste0(pathS, "analyseData.R"))
source(paste0(pathS, "projectInfo.R"))

# (install and) load packages
packages <- c("lavaan", 
              "tidyr"
)
newPackages <-
  packages[!(packages %in% installed.packages()[, "Package"])]
if (length(newPackages)){install.packages(newPackages)}
lapply(packages, require, character.only = TRUE)

# get simulation design
pI <- projectInfo() # when the function is executed the first time, SimConds.rds is created
TableConds <- readRDS("objects/SimConds.rds")
nconds <- nrow(TableConds)

# use multiple cores through initializing multiple R sessions
ncores <- 20
nrep <- pI$nrep
scriptReps <- seq(script, nrep, ncores)
set.seed(as.numeric(paste0(script, nrep, ncores, collapse="")))

# Do you want to center the data in LF and WF?
center <- FALSE

for (rep in scriptReps){
  for (cond in 1:nconds){
    
    # generate data
    n <- TableConds$n[cond]
    p <- TableConds$p[cond]
    g <- TableConds$g[cond]
    data <- simulateData(N=TableConds$N[cond], n=n, g=g, p=p, var_B=TableConds$ICC[cond], var_W=(1-TableConds$ICC[cond]))
    saveRDS(data, file = paste0(pathD, "data_C", cond, "_R", rep, ".rds"))

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

    # analyse data
    fit <- analyseData(data_LF=data_LF,
                       data_WF=data_WF,
                       p=p, n=n, g=g, center=center)
    saveRDS(fit, file = paste0(pathF, "fit_C", cond, "_R", rep, ".rds"))
    
  }
  gc()
  print(rep)
}
