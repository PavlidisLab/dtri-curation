constructBedsDataLoader <- function() {
  
  public <- list()
  private <- list()
  
  # Transcriptomic datasets
  public$datasetIds$sun_2017 <- "GSE76315"
  # =====
  
  private$bedsDataDir <- paste0(globalSession$static$HOME_DIR, "beds_data/clean/")
  
  # Data loader function
  public$loadData <- function(datasetId) {
    return(readRDS(paste0(private$bedsDataDir, datasetId, ".rds")))
  }
  
  return(public)
}
