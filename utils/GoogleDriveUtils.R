constructGoogleDrive <- function(remote, local) {
  
  public <- list()
  
  public$remote <- remote
  public$local <- local
  
  public$files <- names(public$remote) %>% lapply(function(currName) { currName })
  names(public$files) <- public$files
  
  public$loadFileData <- function(file) {
    localFile <- public$local[[file]]
    return(session$dataWrangler$readTsv(localFile))
  }
  
  public$pullFromDrive <- function(file) {
    localFile <- public$local[[file]]
    remoteFile <- public$remote[[file]]
    if (file.exists(localFile)) {
      file.copy(localFile, paste0(localFile, "_OLD"), overwrite = TRUE)
    }
    drive_download(remoteFile, localFile, overwrite = TRUE)
    public$loadFileData(file)
  }
  
  public$reversePull <- function(file) {
    localFile <- public$local[[file]]
    oldLocalFile <- paste0(localFile, "_OLD")
    file.copy(oldLocalFile, localFile, overwrite = TRUE)
    file.remove(oldLocalFile)
    public$loadFileData(file)
  }
  
  public$pushToDrive <- function(tibble, file) {
    localFile <- public$local[[file]]
    remoteFile <- public$remote[[file]]
    if (file.exists(localFile)) {
      file.copy(localFile, paste0(localFile, "_OLD"), overwrite = TRUE)
    }
    tibble %>% write.csv(localFile, row.names = FALSE)
    drive_update(remoteFile, localFile)
    public$loadFileData(file)
  }
  
  public$reversePush <- function(file) {
    localFile <- public$local[[file]]
    oldLocalFile <- paste0(localFile, "_OLD")
    remoteFile <- public$remote[[file]]
    file.copy(oldLocalFile, localFile, overwrite = TRUE)
    file.remove(oldLocalFile)
    drive_update(remoteFile, localFile)
    public$loadFileData(file)
  }
  
  public$updateFileDates <- function(oldVersion) {
    allFiles <- paste0(oldVersion, "_", public$remote %>% names(), ".csv")
    newFiles <- paste0(public$currentVersion, "_", public$remote %>% names(), ".csv")
    for (i in 1:length(allFiles)) {
      currFilePath <- paste0(public$syncDir, allFiles[i])
      newFilePath <-  paste0(public$syncDir, newFiles[i])
      file.rename(currFilePath, newFilePath)
    }
  }
  
  return(public)
}