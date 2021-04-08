constructOntology <- function(session) {
  
  
  private <- list()
  
  private$remoteSource <- "/space/grp/nlim/CronGemmaDump/Ontology/"
  private$localSource <- paste0(session$UTILS_DIR, "OntologyUtils_DATA/")
  
  private$allFiles <- system(paste0("ls ", private$remoteSource), intern = TRUE)
  
  private$rewnewData <- function() {
    remoteFiles <- paste0(private$remoteSource, private$allFiles)
    localFiles <- paste0(private$localSource, private$allFiles)
    for (i in 1:length(private$allFiles)) {
      system(paste0("cp ", remoteFiles[i], " ", localFiles[i]))
    }
  }
  
  public <- list()
  
  public$ontologies <- str_extract(grep("_DEF.TSV$", private$allFiles, value = TRUE), "(?<=Ontology_Dump_).*(?=_DEF.TSV)") %>% unique()
  names(public$ontologies) <- public$ontologies
  public$ontologies <- public$ontologies %>% as.list()
  
  
  public$definitions <- public$ontologies %>% session$collectionUtils$lapply(function(currOntology) {
    session$dataWrangler$readTsv(paste0(private$localSource, "Ontology_Dump_", currOntology, "_DEF.TSV")) %>% 
      dplyr::select(node = Node, label = Definition, scope = OntologyScope)
  }, reporter = "OntologyUtils::definitions")
  
  
  public$nodeLabels <- public$definitions$MERGED %>% 
    dplyr::select(node, label) %>% 
    unique() %>% 
    na.omit() %>% 
    group_by(node) %>% 
    summarize(label = dplyr::first(label)) %>% 
    unique()
  
  public$relationships <-  public$ontologies %>% session$collectionUtils$lapply(function(currOntology) {
    session$dataWrangler$readTsv(paste0(private$localSource, "Ontology_Dump_", currOntology, ".TSV")) %>% 
      dplyr::select(parent = ParentNode, child = ChildNode, relationship = RelationType, scope = OntologyScope)
  }, reporter = "OntologyUtils::relationships")
  
  
  public$getChildNodeRecursive <- function(nodes, scope) {
    oldNodes <- nodes %>% unique()
    currNodes <- public$relationships[[scope]] %>% 
      filter(parent %in% oldNodes) %>% 
      session$dataWrangler$extractColumn("child") %>% 
      unique()
    newNodes <- currNodes %>% setdiff(oldNodes)
    if (length(newNodes) > 0) {
      return(public$getChildNodeRecursive(c(oldNodes, newNodes), scope))
    } else {
      return(oldNodes)
    }
  }
  
  public$getChildNodes <- function(node, scope) {
    public$relationships[[scope]] %>% 
      filter(parent %in% node) %>% 
      left_join(public$definitions[[scope]] %>% dplyr::select(parent = node, parent_label = label), by = "parent") %>% 
      left_join(public$definitions[[scope]] %>% dplyr::select(child = node, child_label = label), by = "child")
  }
  
  public$getParentNodes <- function(node, scope) {
    public$relationships[[scope]] %>% 
      filter(child %in% node) %>% 
      left_join(public$definitions[[scope]] %>% dplyr::select(parent = node, parent_label = label), by = "parent") %>% 
      left_join(public$definitions[[scope]] %>% dplyr::select(child = node, child_label = label), by = "child")
  }
  
  public$getLabels <- function(nodes, scope) {
    public$definitions[[scope]] %>% filter(node %in% nodes)
  }
  
  public$estimateRootDistance <- function(node, scope) {
    parents <- public$relationships[[scope]] %>% 
      filter(child %in% node) %>% 
      session$dataWrangler$extractColumn("parent") %>% 
      unique()
    if (length(parents) <= 0) {
      return(0)
    } else {
      # distanceRecursive <- parents %>% sapply(function(currParent) { public$estimateRootDistance(currParent, scope) }) %>% min() # assumption: no circularity
      distanceRecursive <- public$estimateRootDistance(parents[1], scope)
      return(1 + distanceRecursive)
    }
  }
  
  return(public)
}
