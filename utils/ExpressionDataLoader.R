constructExpressionDataLoader <- function() {
  
  public <- list()
  private <- list()
  
  # Transcriptomic datasets
  public$datasetIds$brainspan_2011 <- "GSE25219"
  public$datasetIds$braincloud_2011 <- "GSE30272"
  public$datasetIds$sun_2017 <- "GSE76315"
  public$datasetIds$swisa_2016 <- "GSE87530"
  # =====
  
  private$expressionDataDir <- paste0(globalSession$static$HOME_DIR, "expression_data/clean/")
  
  # Data loader function
  public$loadData <- function(datasetId) {
    return(readRDS(paste0(private$expressionDataDir, datasetId, ".rds")))
  }
  
  return(public)
}
