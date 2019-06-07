mean_rm.na <- 
  ## Standard mean function with na.rm bound to TRUE.
  function(X){
    return(mean(X, na.rm=TRUE))
  }

sd_rm.na <- function(X){
  ## Standard sd function with na.rm bound to TRUE.
  return(sd(X, na.rm=TRUE))
}
