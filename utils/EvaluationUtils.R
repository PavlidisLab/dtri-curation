constructEvaluationUtils <- function() {
  
  public <- list()
  
  private <- list()
  
  private$computeRecall <- function(p, n, t) {
    tp <- p %>% intersect(t) %>% length()
    fn <- t %>% intersect(n) %>% length()
    return(tp / (tp + fn))
  }
  
  private$computeFPR <- function(p, n, t) {
    fp <- p %>% setdiff(t) %>% length()
    tn <- n %>% setdiff(t) %>% length()
    return(fp / (fp + tn))
  }
  
  private$computePrecision <- function(p, n, t) {
    tp <- p %>% intersect(t) %>% length()
    fp <- p %>% setdiff(t) %>% length()
    return(tp /(tp + fp))
  }
  
  public$roc <- function(ranking, trueSet) {
    t <- trueSet %>% intersect(ranking)
    roc <- do.call(rbind, 1:length(ranking) %>% session$collectionUtils$lapply(function(i) {
      p <- ranking[1:i]
      n <- ranking %>% setdiff(p)
      tibble(recall = private$computeRecall(p, n, t), 
             false_positive_rate = private$computeFPR(p, n, t))
    }, reporter = "ROC"))
    return(roc)
  }
  
  public$precisionRecall <- function(ranking, trueSet) {
    t <- trueSet %>% intersect(ranking)
    precisionRecall <- do.call(rbind, 1:length(ranking) %>% session$collectionUtils$lapply(function(i) {
      p <- ranking[1:i]
      n <- ranking %>% setdiff(p)
      tibble(recall = private$computeRecall(p, n, t), 
             precision = private$computePrecision(p, n, t))
    }, reporter = "PR"))
    return(precisionRecall)
  }
  
  public$plotAuc <- function(x, y) {
    data <- tibble(x = x, y = y)
    data <- data %>% 
      mutate(x_diff = diff(c(0, x))) %>% 
      mutate(unit_area = x_diff * y) 
    data$unit_area %>% sum()
  }
  
  public$auroc <- function(ranking, trueSet) {
    setStatus <- tibble(item = ranking) %>% mutate(in_set = item %in% trueSet) %>% session$dataWrangler$extractColumn("in_set")
    auc <- colAUC(1:length(ranking), setStatus)[1,1]
    names(auc) <- NULL
    return(auc)
  }
  
  public$aurocCI <- function(ranking, trueSet, iteration) {
    f <- ranking %>% setdiff(trueSet)
    t <- trueSet
    sampleAuc <- 1:iteration %>% sapply(function(i) {
      currF <- f %>% sample(length(f), replace = TRUE) %>% unique()
      currT <- t %>% sample(length(t), replace = TRUE) %>% unique()
      currRanking <- ranking[which(ranking %in% c(currF, currT))]
      public$auroc(currRanking, currT)
    })
    absCI <- (sd(sampleAuc)*1.96)
    mean <- public$auroc(ranking, trueSet)
    return(tibble(lower = mean - absCI, mean = mean, upper = mean + absCI))
  }
  
  
  public$precision <- function(ranking, trueSet, recall) {
    
    recallN <- round((length(trueSet) * recall), 0)
    
    trueRanks <- tibble(item = ranking, rank = 1:length(ranking)) %>% 
      mutate(in_set = item %in% trueSet) %>% 
      filter(in_set)
    
    recallRank <- trueRanks[recallN,]$rank
    
    p <- ranking[1:recallRank]
    n <- ranking %>% setdiff(p)
    t <- trueSet
    
    private$computePrecision(p, n, t)
  }
  
  
  public$precisionCI <- function(ranking, trueSet, recall, iteration) {
    f <- ranking %>% setdiff(trueSet)
    t <- trueSet
    samplePrecision <- 1:iteration %>% sapply(function(i) {
      currF <- f %>% sample(length(f), replace = TRUE) %>% unique()
      currT <- t %>% sample(length(t), replace = TRUE) %>% unique()
      currRanking <- ranking[which(ranking %in% c(currF, currT))]
      public$precision(currRanking, currT, recall)
    })
    absCI <- (sd(samplePrecision)*1.96)
    mean <- public$precision(ranking, trueSet, recall)
    return(tibble(lower = mean - absCI, mean = mean, upper = mean + absCI))
  }
  
  
  public$computeOra <- function(hits, trueSet, background) {
    
    allPositives <- trueSet %>% intersect(background) %>% na.omit()
    allHits <- hits %>% intersect(background) %>% na.omit()
    
    truePositives <- allHits %>% intersect(allPositives)
    
    hits_frac = length(hits) / length(background)
    tp_frac = length(truePositives) / length(allPositives)
    
    pvalue = phyper(length(truePositives), length(allPositives), length(background) - length(allPositives), length(hits), lower.tail = FALSE)
    
    tibble(n_background = length(background), 
           n_sample = length(hits), 
           n_true = length(allPositives), 
           n_tp = length(truePositives), 
           sample_frac = hits_frac, 
           tp_frac = tp_frac,
           p_value = pvalue)
  }
  
  
  return(public)
}












  