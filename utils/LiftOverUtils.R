library(liftOver)

constructLiftOverUtils <- function(session) {
  
  public <- list()
  private <- list()
  
  private$localSource <- paste0(session$UTILS_DIR, "LiftOverUtils_DATA/")
  
  public$hg19ToHg38 <- function(gRanges) {
    chain <- import.chain(paste0(private$localSource, "hg19ToHg38.over.chain"))
    result <- liftOver(gRanges, chain)
    result %>% unlist() 
  }
  
  return(public)
}
