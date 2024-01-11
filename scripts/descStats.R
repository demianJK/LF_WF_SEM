## Table 1: Descriptive Statistics of the Model Performance Criteria

# load packages
library(dplyr)
library(huxtable)

# set paths and load data
pathO <- "objects/"
dat <- readRDS(paste0(pathO, "dat.rds"))

# select conditions in which WF (and thus, LF as well) converged 100% (for descriptives estimation accuracy)
which <- filter(dat, approach == "WF", conv == 100)$cond
datWFconv <- filter(dat, cond %in% which)

# compute descriptive statistics for all performance measures
desStat <- c() 
desStat <- rbind(desStat, psych::describe(dat[,"conv"]))
desStat <- rbind(desStat, psych::describe(datWFconv[,"relRMSE"]))
desStat <- rbind(desStat, psych::describe(datWFconv[,"relBias"]))
desStat <- rbind(desStat, psych::describe(datWFconv[,"relVar"]))

# using only simulation conditions that resulted in 100% convergence in both approaches does not alter the descriptive stats substantially
# check with:
# desStat <- rbind(desStat, psych::describe(dat[,"relRMSE"]))
# desStat <- rbind(desStat, psych::describe(dat[,"relBias"]))
# desStat <- rbind(desStat, psych::describe(dat[,"relVar"]))

# set nb of decimals in table
ndec <- 2 # number of decimals

desStat <- select(desStat, mean, sd, median, min, max)
desc <- data.frame(crit = c("Convergence Rate", "Relative RMSE", "Relative Bias", "Relative Variance"))
desc <- cbind(desc, desStat)
desc <- as_hux(desc)
number_format(desc) <- ndec
desc <- map_align(desc, by_cols(".")) # align at decimal
desc <- set_align(desc, 1:nrow(desc), 1, "left") # first col left-aligned
desc <- set_align(desc, 1, 1:ncol(desc), "center") # all headers centered
desc <- set_contents(desc, 1, 1:ncol(desc), c("Criterion", "M", "SD", "Median", "Min", "Max")) # header names
desc <- set_italic(desc, 1, 2:6, value = TRUE)
desc <- set_header_rows(desc, 1, TRUE)
desc <- theme_article(desc)
desc <- style_headers(desc, bold=FALSE)
print_screen(desc)
quick_latex(desc, file="objects/tab1.tex")
