projectInfo <- function(){
  
  out <- list(
    
    approaches = c("LF", "WF"),
    nrep = 1000,
    
    # population covariance matrix
    p = c(2, 5, 10, 20),
    ICC = c(0.05, 0.25, 0.5),
    
    # sample characteristics
    N = c(20, 50, 100, 200, 1000),
    n = c(2, 5, 10, 20, 100)
  )
  

  if ( !file.exists("objects/SimConds.rds") ){ # create Sim Table if it does not exist

    # cross population covariance matrix and sample characteristics
    final <- expand.grid(p=out$p, ICC=out$ICC, N=out$N, n=out$n) 
    final$g <- final$N/final$n
    
    ## sort out impossible conditions 
    # for sample characteristics: N=n, N>n, g not an integer, g<2
    # for population x sample characteristics: p>N
    ind <- c()
    for (i in 1:nrow(final)){
      if (final$n[i] == final$N[i] || final$n[i] > final$N[i] || final$g[i]%%1!=0 || final$g[i]<2 || final$p[i] > final$N[i]){ 
        ind <- c(ind, i)
      }
    }
    final <- final[-c(ind),]
    
    final <- final[order(final$p),] # sort with increasing p
    cond <- 1:nrow(final) # give simulation conditions nb
    table <- cbind(cond, final)
    xlsx::write.xlsx( table, file="objects/SimConds.xlsx" )
    saveRDS(table, file = "objects/SimConds.rds" ) 
  }
  
  return(out)
  
}


