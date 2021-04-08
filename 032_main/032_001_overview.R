# run 032_SETUP.R first

ws032$m001 <- list()
ws032$m001$figurePath <- paste0(ws032$figuresDir, "f01_")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - overall curation workflow

pData <- list()

pData$candidatePool <- ws032$data$s01.tsv

pData$curatedExperiments <- ws032$data$s03.tsv

pData$tfs <- ws032$data$s02.tsv

# number of papers from external sources
pData$candidatePool %>% filter(Database != "PubMedQuery") %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

# number of TFs from external sources
pData$candidatePool %>% filter(Database != "PubMedQuery") %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique() %>% length()

# number of papers from our own query
pData$candidatePool %>% filter(Database == "PubMedQuery") %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

# number of papers from our own query
fData <- list()
fData$queryOnly <- pData$candidatePool %>% filter(Database == "PubMedQuery") %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% 
  setdiff((pData$candidatePool %>% filter(Database != "PubMedQuery") %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique()))
fData$queryOnly %>% length()
pData$candidatePool %>% 
  filter(PubMed_ID %in% fData$queryOnly) %>% 
  session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique() %>% length()

# total number of papers
pData$candidatePool %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

# total number of TFs
pData$candidatePool %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique() %>% length()

# total number of TFs with no paper
pData$tfs %>% filter((TF_Entrez_ID %in% pData$candidatePool$TF_Entrez_ID))
pData$tfs %>% filter(!(TF_Entrez_ID %in% pData$candidatePool$TF_Entrez_ID)) %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique() %>% length()


# total number of papers curated or discarded
pData$candidatePool %>%
  filter(Status %in% c("curated", "discarded")) %>% 
  session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

# total number of TFs whose papers were curated or discarded
pData$candidatePool %>%
  filter(Status %in% c("curated", "discarded")) %>% 
  session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique() %>% length()

# number of papers curated
pData$candidatePool %>%
  filter(Status %in% c("curated")) %>% 
  session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

# number of experiments curated
pData$curatedExperiments %>%
  session$dataWrangler$extractColumn("Experiment_ID") %>% unique() %>% length()

# number of interactions curated
pData$curatedExperiments %>%
  session$dataWrangler$extractColumn("DTRI_ID") %>% unique() %>% length()

# number of TFs curated
pData$curatedExperiments %>%
  session$dataWrangler$extractColumn("TF_Entrez_ID_Human") %>% unique() %>% length()



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - show off the network


# interactions (edges)
# 1. multiple exp vs single exp
# 2. cns vs. other

# genes (nodes)
# 1. is_tf
# 2. is_nd_tf
# 3. n_target


pData <- list()

pData$outputFile <- paste0(ws032$m001$figurePath, "01")


pData$edges <- ws032$data$s03.tsv %>% 
  dplyr::select(Experiment_ID, DTRI_ID, TF_Entrez_ID_Human, Target_Entrez_ID_Human, 
                TF_Symbol_Human, Target_Symbol_Human, Experiment_Type, Cell_Type)

pData$cnsTerms <- session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$UBERON) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$CL)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$UBERON)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$CL)) %>%
  unique() # include both "central nervous system" and "future central nervous system"

pData$edgeCns <- pData$edges %>% 
  mutate(cns = Cell_Type %in% pData$cnsTerms) %>% 
  group_by(DTRI_ID) %>% 
  summarize(cns = any(cns)) %>% 
  ungroup()

pData$edgeExpType <- pData$edges %>% 
  dplyr::select(DTRI_ID, Experiment_Type) %>% 
  unique() %>% 
  group_by(DTRI_ID) %>% 
  summarize(n_exp_type = n()) %>% 
  ungroup()

pData$edges <- pData$edges %>% 
  dplyr::select(DTRI_ID, TF_Entrez_ID_Human, Target_Entrez_ID_Human, TF_Symbol_Human, Target_Symbol_Human) %>% 
  unique() %>% 
  left_join(pData$edgeCns, by = "DTRI_ID") %>% 
  left_join(pData$edgeExpType, by = "DTRI_ID") %>% 
  arrange(desc(cns), desc(n_exp_type))

pData$nodes <- tibble(Entrez_ID = ws032$data$s03.tsv$TF_Entrez_ID_Human %>% c(ws032$data$s03.tsv$Target_Entrez_ID) %>% unique())

pData$nodesNTarget <- ws032$data$s03.tsv %>% 
  dplyr::select(TF_Entrez_ID_Human, Target_Entrez_ID) %>% 
  unique() %>% 
  group_by(TF_Entrez_ID_Human) %>% 
  summarize(N_Target = n()) %>% 
  dplyr::select(Entrez_ID = TF_Entrez_ID_Human, everything()) %>% 
  arrange(desc(N_Target))

pData$nodesTf <- ws032$data$s03.tsv %>% 
  dplyr::select(TF_Entrez_ID_Human) %>% 
  unique() %>% 
  mutate(Is_TF = TRUE) %>% 
  dplyr::select(Entrez_ID = TF_Entrez_ID_Human, everything())

pData$nodesNdTf <- ws032$data$s03.tsv %>% 
  dplyr::select(TF_Entrez_ID_Human) %>% 
  unique() %>% 
  left_join(ws032$data$s02.tsv %>% dplyr::select(TF_Entrez_ID_Human = TF_Entrez_ID, Neurodev_TF) %>% unique(), by = "TF_Entrez_ID_Human") %>% 
  dplyr::select(Entrez_ID = TF_Entrez_ID_Human, everything())

pData$nodes <- pData$nodes %>% 
  left_join(pData$nodesNTarget, by = "Entrez_ID") %>% 
  session$dataWrangler$fillNa("n_target", 0) %>% 
  left_join(pData$nodesTf, by = "Entrez_ID") %>% 
  left_join(pData$nodesNdTf, by = "Entrez_ID") %>% 
  session$dataWrangler$fillNa(c("is_tf", "nd_tf"), FALSE)

pData$nodes <- pData$nodes %>% 
  left_join(ws032$genesLabs %>% dplyr::select(Entrez_ID = entrez, name = gene), by = "Entrez_ID")


print(paste0("output: ", paste0(pData$outputFile, "_edges.tsv")))
pData$edges %>% write_tsv(paste0(pData$outputFile, "_edges.tsv"))

print(paste0("output: ", paste0(pData$outputFile, "_nodes.tsv")))
pData$nodes %>% write_tsv(paste0(pData$outputFile, "_nodes.tsv"))


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# compute stats about the largest component of the curated network

pData <- list()

pData$network <- ws032$data$s04.tsv %>% 
  dplyr::select(source = TF_Entrez_ID_Human, target = Target_Entrez_ID_Human) %>% 
  unique() %>% 
  filter(source != target)

pData$graph <- pData$network %>% 
  mutate(edge = paste0(source, ",", target)) %>% 
  session$dataWrangler$extractColumn("edge") %>% 
  paste0(collapse = ",") %>% 
  strsplit(",") %>% 
  unlist() %>% 
  make_graph()

pData$genesMainComponent <- components(pData$graph)$membership %>% 
  as.data.frame() %>% 
  session$dataWrangler$setRownameAsColumn("gene") %>% 
  filter(. == 1) %>% 
  session$dataWrangler$extractColumn("gene")

pData$networkMainComponent <- ws032$data$s04.tsv %>% 
  dplyr::select(source = TF_Entrez_ID_Human, target = Target_Entrez_ID_Human) %>% 
  unique() %>% 
  filter(source != target) %>% 
  filter((source %in% pData$genesMainComponent) & (target %in% pData$genesMainComponent))

pData$networkMainComponent %>% nrow() # number of edges 

pData$networkMainComponent$source %>% c(pData$networkMainComponent$target) %>% unique() %>% length() # number of nodes


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# main text - overview

pData <- list()

pData$candidatePool <- ws032$data$s01.tsv
pData$curatedExperiments <- ws032$data$s03.tsv
pData$curatedDtris <- ws032$data$s04.tsv

# number of TFs
pData$curatedDtris %>% 
  session$dataWrangler$extractColumn("TF_Symbol_Human") %>% unique() %>% length()

# number of targets
pData$curatedDtris %>% 
  session$dataWrangler$extractColumn("Target_Symbol_Human") %>% unique() %>% length()

# number of papers curated
pData$curatedExperiments %>% 
  session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()
pData$candidatePool %>% 
  filter(Status == "curated") %>% 
  session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()


