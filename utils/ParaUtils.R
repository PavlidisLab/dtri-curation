constructParaUtils <- function() {
  
  public <- list()
  
  public$lapply <- function(x, func, nCores = 10) {
    cl <- makeCluster(10, type = "FORK")
    result <- parLapply(cl, x, func)
    stopCluster(cl)
    return(result)
  }
  
  public$sapply <- function(x, func, nCores = 10) {
    cl <- makeCluster(10, type = "FORK")
    result <- parSapply(cl, x, func)
    stopCluster(cl)
    return(result)
  }
  
  return(public)
}
