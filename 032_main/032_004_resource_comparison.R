# run 032_SETUP.R first

ws032$m004 <- list()

ws032$m004$figurePath <- paste0(ws032$figuresDir, "f04_")

ws032$m004$dtris <- ws032$data$s09.tsv %>% dplyr::select(DTRI_ID, Database) %>% 
  rbind(ws032$data$s04.tsv %>% dplyr::select(DTRI_ID) %>% mutate(Database = "Current")) %>% 
  unique()

ws032$m004$neuroTfs <- ws032$data$s02.tsv %>% filter(Neurodev_TF == "TRUE") %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique()

ws032$m004$s03 <- ws032$data$s03.tsv %>% 
  mutate(species_str = paste0(TF_Species, ",", Target_Species)) %>% 
  mutate(Species = species_str %>% sapply(function(currSpecies) {
    mouse <- grepl("mouse", currSpecies)
    human <- grepl("human", currSpecies)
    if (mouse & human) { "mixed" }
    else if (mouse) { "mouse" }
    else if (human) { "human" }
    else { currSpecies }
  })) %>% dplyr::select(-TF_Species, -Target_Species, -species_str) %>% 
  unique()

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - flowchart of previously curated papers

pData <- list()

pData$candidates <- ws032$data$s01.tsv
pData$candidatesExternal <- pData$candidates %>% filter(Database != "PubMedQuery")
pData$curatedExperiments <- ws032$data$s03.tsv

# total number of pmids from external Databases
pData$candidates %>% filter(Database != "PubMedQuery") %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

pData$candidatesExternal %>% 
  filter(Status != "unexamined") %>% 
  session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

pData$candidatesExternal %>% 
  filter(Status != "unexamined") %>% 
  session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique() %>% length()

pData$candidatesExternal%>% 
  filter(Status == "curated") %>% 
  session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

pData$candidatesExternal %>% 
  filter(Status == "curated") %>% 
  session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique() %>% length()


fData <- list()
fData$curatedExternal <- pData$candidatesExternal %>% 
  filter(Status == "curated") %>% 
  session$dataWrangler$extractColumn("PubMed_ID") %>% unique() 
pData$curatedExperiments %>% 
  filter(PubMed_ID %in% fData$curatedExternal) %>% 
  session$dataWrangler$extractColumn("DTRI_ID") %>% 
  unique() %>% length()


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - distribution of DTRIs by number of Databases

pData <- list()

pData$outputFile <- paste0(ws032$m004$figurePath, "01.png")

pData$records <- ws032$m004$dtris %>%
  dplyr::select(tri_id = DTRI_ID, database = Database) %>%
  unique() 

pData$nDb <- pData$records %>% 
  group_by(tri_id) %>% 
  summarize(n_database = n()) %>% 
  arrange(desc(n_database))

pData$nDb <- pData$nDb %>% 
  group_by(n_database) %>% 
  summarize(n_interaction = n())

pData$nDb %>% 
  session$graphingUtils$ggplot(aes(x = n_database, y = n_interaction)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n_interaction), vjust = -1, size = 6.5) +
  ylim(0, 21000) +
  xlab("Number of Databases") +
  ylab("Number of DTRIs") +
  theme(legend.position = "none", axis.text = element_text(size = 15))


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile,
       width = 4.5, height = 6)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - venn diagram of number of DTRIs recorded in external databases vs by current

pData <- list()

pData$outputFile <- paste0(ws032$m004$figurePath, "02.png")

pData$interactions <- ws032$m004$dtris %>% 
  dplyr::select(tri_id = DTRI_ID, database = Database) %>%
  unique() %>% 
  mutate(database = database %>% sapply(function(currDb) {
    if (currDb == "Current") { "Current" }
    else { "DBs" }
  })) %>% 
  unique()

pData$sets <- pData$interactions$database %>% unique()
names(pData$sets) <- pData$sets
pData$sets <- pData$sets %>% lapply(function(currDb) {
  pData$interactions %>% 
    filter(database == currDb) %>% 
    session$dataWrangler$extractColumn("tri_id")
}) 

print(paste0("printing: ", pData$outputFile))
venn.diagram (
  x = pData$sets,
  filename = pData$outputFile,
  output = TRUE, 
  width = 3300, 
  height = 3000, 
  cex = 3.5, 
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0,
  lwd = 2,
  lty = 'blank',
  fill = c("#DCF0E7", "#FDE8D9")
)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - venn diagram of number of DTRIs recorded in external databases vs by current, filter for ND TF

pData <- list()

pData$outputFile <- paste0(ws032$m004$figurePath, "03.png")

pData$interactions <- ws032$m004$dtris %>% 
  left_join(ws032$data$s08.tsv %>% dplyr::select(DTRI_ID, TF_Entrez_ID_Human), by = "DTRI_ID") %>% 
  filter(TF_Entrez_ID_Human %in% ws032$m004$neuroTfs) %>% 
  dplyr::select(tri_id = DTRI_ID, database = Database) %>%
  unique() %>% 
  mutate(database = database %>% sapply(function(currDb) {
    if (currDb == "Current") { "Current" }
    else { "DBs" }
  })) %>% 
  unique()

pData$sets <- pData$interactions$database %>% unique()
names(pData$sets) <- pData$sets
pData$sets <- pData$sets %>% lapply(function(currDb) {
  pData$interactions %>% 
    filter(database == currDb) %>% 
    session$dataWrangler$extractColumn("tri_id")
}) 


print(paste0("printing: ", pData$outputFile))
venn.diagram (
  x = pData$sets,
  filename = pData$outputFile,
  output = TRUE, 
  width = 3300, 
  height = 3000, 
  cex = 3, 
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = 0,
  lwd = 2,
  lty = 'blank',
  fill = c("#DCF0E7", "#FDE8D9")
)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# f04_05 - number DTRIs (without considering overlaps between sources), simple bar, fill = nd_tf

pData <- list()

pData$outputFile <- paste0(ws032$m004$figurePath, "04.png")

pData$records <- ws032$m004$dtris %>%
  dplyr::select(tri_id = DTRI_ID, database = Database) %>% 
  left_join(ws032$data$s08.tsv %>% dplyr::select(tri_id = DTRI_ID, tf_gene_h = TF_Entrez_ID_Human), by = "tri_id") %>% 
  dplyr::select(database, tri_id, tf_gene_h) %>% 
  unique() %>% 
  mutate(neurodev = tf_gene_h %>% sapply(function(currTf) { currTf %in% ws032$m004$neuroTfs})) %>% 
  dplyr::select(database, tri_id, neurodev) %>% 
  unique()

pData$bar <- pData$records  %>% 
  group_by(database, neurodev) %>% 
  summarize(n_dtri = n()) %>% 
  ungroup() %>% 
  arrange(desc(n_dtri))

pData$sum <- pData$bar %>% 
  group_by(database) %>% 
  summarize(n_dtri = sum(n_dtri)) %>% 
  arrange(desc(n_dtri))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = database, y = n_dtri)) +
  geom_bar(aes(fill = neurodev), stat = "identity") +
  scale_fill_manual(values = c("grey80", "#3182BC")) +
  geom_text(data = pData$sum, aes(label = n_dtri), hjust = -0.2, size = 8) +
  ylim(0, 18000) + 
  coord_flip() + 
  scale_x_discrete(limits = rev(c("Current", pData$sum$database[!grepl("Current", pData$sum$database)])), 
                   labels = rev(c("Current", "TRRUST", "InnateDB", "TFactS", "HTRIdb_LC", "TFe", "CytReg", "ENdb", "ORegAnno_LC"))) +
  ylab("Number of DTRIs") +
  xlab("") +
  theme(legend.position = "none", axis.text = element_text(size = 20))

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 10, height = 6)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# CONSOLIDATING DTRI COMBINATION OF EXPERIMENTAL APPROACHES

pData <- list()

pData$perturbation <- ws032$data$s05.tsv %>% 
  left_join(ws032$m004$s03, by = "Experiment_ID") %>% 
  dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, species = Species, context_type = Context_Type, effect = Effect, induced = Induced) %>% 
  unique()

pData$perturbation <- pData$perturbation %>% 
  group_by(tri_id) %>% 
  summarize(effect_str = paste0(effect, collapse = ","), 
            induced_str = paste0(induced, collapse = ",")) %>% 
  mutate(perturbation = 1) %>% 
  dplyr::select(-effect_str, -induced_str)

pData$binding <- ws032$data$s06.tsv %>% 
  left_join(ws032$m004$s03, by = "Experiment_ID") %>% 
  dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, species = Species, context_type = Context_Type, 
                method = Method, tf_source_type = TF_Source_Type) %>% 
  unique()

pData$binding <- pData$binding %>% 
  group_by(tri_id) %>% 
  summarize(method_str = paste0(method, collapse = ","), 
            tf_source_type_str = paste0(tf_source_type, collapse = ",")) %>% 
  mutate(binding = 1) %>% 
  dplyr::select(-method_str, -tf_source_type_str)

pData$reporter <- ws032$data$s07.tsv %>% 
  left_join(ws032$m004$s03, by = "Experiment_ID")  %>% 
  dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, species = Species, context_type = Context_Type, 
                mutated = Mutated, binding_verified = Binding_Verified) %>% 
  unique()

pData$reporter <- pData$reporter %>% 
  group_by(tri_id) %>% 
  summarize(mutated_str = paste0(mutated, collapse = ",")) %>% 
  mutate(reporter = 1) %>% 
  dplyr::select(-mutated_str)

pData$overall <- ws032$m004$s03 %>% 
  dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, species = Species, context_type = Context_Type, cell_type = Cell_Type) %>% 
  unique()

pData$cnsTerms <- session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$UBERON) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$CL)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$UBERON)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$CL)) %>% 
  unique()

pData$overall <- pData$overall %>% 
  mutate(cns = (cell_type %in% pData$cnsTerms)) %>% 
  group_by(tri_id) %>% 
  summarize(context_type_str = paste0(context_type, collapse = ","), 
            species_str = paste0(species, collapse = ","), 
            cns = any(cns)) %>% 
  mutate(overall.cns = as.numeric(cns), 
         overall.primary_tissue = as.numeric(grepl("primary_tissue", context_type_str)),
         overall.primary_cell = as.numeric(grepl("primary_cell", context_type_str)), 
         overall.primary = as.numeric(grepl("primary", context_type_str)), 
         overall.cell_line = as.numeric(grepl("cell_line", context_type_str)), 
         overall.species = species_str %>% sapply(function(currStr) {
           human <- grepl("human", currStr)
           mouse <- grepl("mouse", currStr)
           mixed <- grepl("mixed", currStr)
           if ((human & mouse) | mixed) { "both" }
           else if (human) { "human" }
           else if (mouse) { "mouse" }
           else { currStr }
         })) %>% 
  dplyr::select(-context_type_str, -species_str, -cns)


pData$sets <- c("perturbation", 
                "binding", 
                "reporter",
                "overall.primary_tissue", 
                "overall.primary_cell",
                "overall.primary", 
                "overall.cell_line", 
                "overall.cns")

pData$master <- pData$perturbation %>% 
  full_join(pData$binding, by = "tri_id") %>% 
  full_join(pData$reporter, by = "tri_id") %>% 
  full_join(pData$overall, by = "tri_id")

pData$master <- pData$master %>% 
  session$dataWrangler$fillNa(pData$sets, 0)


# ------------------------------------------------------------------------------
# FIGURE - Number of DTRIs done for different evidence combinations uniquely curated DTRIs in our db

pData$outputFile <- paste0(ws032$m004$figurePath, "05.pdf")

pData$prevDtris <- ws032$data$s09.tsv %>% 
  session$dataWrangler$extractColumn("DTRI_ID") %>% 
  unique()

pData$data <- pData$master %>% 
  filter(!(tri_id %in% pData$prevDtris)) %>% 
  dplyr::select(perturbation, binding, reporter, overall.primary, overall.cns)

pData$plot <- pData$data %>% 
  mutate(overall.primary = overall.primary %>% sapply(function(currNum) { as.character(as.logical(currNum)) }), 
         overall.cns = overall.cns %>% sapply(function(currNum) { as.character(as.logical(currNum)) })) %>% 
  mutate(context = paste0(overall.primary, overall.cns)) %>% 
  as.data.frame() %>% 
  upset(sets = rev(c("perturbation", "binding", "reporter")),
        queries = list(list(query = elements, params = list("context", c("FALSEFALSE", "TRUEFALSE", "TRUETRUE")), color = "grey85", active = TRUE), 
                       list(query = elements, params = list("context", c("TRUEFALSE", "TRUETRUE")), color = "#95B2D5", active = TRUE), 
                       list(query = elements, params = list("context", c("TRUETRUE")), color = "#5382B3", active = TRUE)),
        order.by = "degree",
        text.scale = c(5, 3, 5, 3, 5, 5), 
        keep.order = TRUE,
        mainbar.y.label = "", 
        sets.x.label = "", 
        point.size = 5, line.size = 1)
pData$plot 

print(paste0("printing: ", pData$outputFile))
pdf(file = pData$outputFile, 
    width = 15, height = 10) # or other device
pData$plot
dev.off()


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# TISSUE DISTRIBUTION 


ws032$m004$tissueExperiments$primaryTissue <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "primary_tissue.tsv"))

ws032$m004$tissueExperiments$primaryCell <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "primary_cell.tsv"))

ws032$m004$tissueExperiments$primaryCombined <- ws032$m004$tissueExperiments$primaryTissue %>%
  rbind(ws032$m004$tissueExperiments$primaryCell) %>%
  dplyr::select(node, label, record_id) %>%
  unique() %>%
  group_by(node) %>%
  mutate(n_experiment = n()) %>%
  ungroup() %>%
  dplyr::select(node, label, n_experiment, record_id)

ws032$m004$tissueExperiments$primaryCombined <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "primary_combined.tsv"))

# map "future central nervous system" -UBERON_0016879 to "central nervous system" -UBERON_0001017
ws032$m004$tissueExperiments$primaryCombined <- ws032$m004$tissueExperiments$primaryCombined %>% 
  mutate(node = node %>% sapply(function(currNode) { if (currNode == "UBERON_0016879") { "UBERON_0001017" } else { currNode } }), 
         label = label %>% sapply(function(currLabel) { if (currLabel == "future central nervous system") { "central nervous system" } else { currLabel } }))

ws032$m004$tissueExperiments$cellLine <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "cell_line.tsv"))

ws032$m004$terms$anatomicalSystems <- tibble(node = c("UBERON_0001017", # CNS
                                                      "UBERON_0001032", # sensory
                                                      "CL_0000988",     # hematopoietic
                                                      "UBERON_0002405", # immune
                                                      "UBERON_0002204", # muscuskeletal
                                                      "UBERON_0000026", # appendage
                                                      "UBERON_0002384", # connective tissue
                                                      "UBERON_0004535", # cardiovascular
                                                      "UBERON_0001004", # respiratory
                                                      "UBERON_0001007", # digestive
                                                      "UBERON_0002416", # integumental
                                                      "UBERON_0001008", # renal 
                                                      "UBERON_0000990" # reproductive
))


ws032$m004$terms$anatomicalSystems <- ws032$m004$terms$anatomicalSystems %>% 
  left_join(session$ontologyUtils$getLabels(ws032$m004$terms$anatomicalSystems$node, scope = session$ontologyUtils$ontologies$CL), by = "node")
ws032$m004$terms$anatomicalSystems <- ws032$m004$terms$anatomicalSystems %>% 
  rbind(tibble(node = "other", label = "other", scope = "CL"))

ws032$m004$terms$cellLineCategories <- tibble(
  node = c("CLO_0000220",
           "CLO_0000166",
           "CLO_0000225",
           "CLO_0000235",
           "CLO_0000257",
           "CLO_0000219",
           "CLO_0000236",
           "CLO_0000212",
           "CLO_0000248",
           "CLO_0000239",
           "CLO_0000253",
           "CLO_0000238",
           "CLO_0000228",
           "CLO_0000328",
           "CLO_0000197",
           "CLO_0000249")) 

ws032$m004$terms$cellLineCategories <- ws032$m004$terms$cellLineCategories %>% 
  left_join(session$ontologyUtils$getLabels(ws032$m004$terms$cellLineCategories$node, scope = session$ontologyUtils$ontologies$CLO), by = "node")
ws032$m004$terms$cellLineCategories <- ws032$m004$terms$cellLineCategories %>% 
  rbind(tibble(node = "other", label = "other", scope = "CLO"))

ws032$m004$terms$cellLineCategories <- ws032$m004$terms$cellLineCategories %>% 
  dplyr::select(node, scope, label.y = label) %>% 
  mutate(label.x = label.y %>% str_extract("(?<=immortal).*(?=cell line cell)") %>% trimws()) %>% 
  session$dataWrangler$mergeColumnsXy(colName = "label")

ws032$m004$terms$cellLineCategories <- ws032$m004$terms$cellLineCategories %>% 
  dplyr::select(node, scope, label.y = label) %>% 
  mutate(label.x = label.y %>% str_extract(".*(?=-derived)") %>% trimws()) %>% 
  session$dataWrangler$mergeColumnsXy(colName = "label")


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - primary combined, n experiment across anatomical systems, fill = nd_tf, simple bar, filter = embryonic

pData <- list()

pData$outputFile <- paste0(ws032$m004$figurePath, "06.png")

# first address coverage

pData$targetExperiments <- ws032$m004$tissueExperiments$primaryCombined %>% filter(node %in% ws032$m004$terms$anatomicalSystems$node) %>% 
  dplyr::select(record_id, node, label) %>% 
  unique()

pData$otherExperiments <- ws032$m004$tissueExperiments$primaryCombined %>% filter(!(record_id %in% pData$targetExperiments$record_id)) %>% unique()
pData$otherExperiments %>% group_by(node) %>% mutate(n_experiment = n()) %>% arrange(n_experiment) %>% ungroup()

pData$otherExperiments <- pData$otherExperiments %>% 
  dplyr::select(record_id) %>% 
  unique() %>% 
  mutate(node = "other", label = "other")

# attach "others"
pData$targetExperiments <- pData$targetExperiments %>% 
  rbind(pData$otherExperiments) %>% 
  unique()

pData$targetExperiments <- pData$targetExperiments %>% 
  left_join(ws032$data$s03.tsv %>% dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, tf_gene_h = TF_Entrez_ID_Human, age = Age), by = "record_id") %>% 
  dplyr::select(record_id, tri_id, tf_gene_h, label, node, age) %>% 
  mutate(nd_tf = (tf_gene_h %in% ws032$m004$neuroTfs)) %>% 
  unique()

# filter for stuff found in other dbs as well ***
pData$targetExperiments <- pData$targetExperiments %>% 
  filter(tri_id %in% (ws032$data$s09.tsv %>% session$dataWrangler$extractColumn("DTRI_ID"))) %>% 
  dplyr::select(-tri_id) %>% 
  unique()
# ===


pData$bar <- pData$targetExperiments %>% 
  group_by(node, nd_tf) %>%
  summarize(label = dplyr::first(label), n_experiment = n()) %>% 
  arrange(desc(n_experiment)) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(node) %>% 
  summarize(n_experiment = sum(n_experiment))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = node, y = n_experiment)) +
  geom_bar(stat = "identity", aes(fill = nd_tf)) +
  geom_text(data = pData$sums, aes(label = n_experiment), hjust = -0.5, size = 6) +
  coord_flip() +
  scale_x_discrete(breaks = ws032$m004$terms$anatomicalSystems$node, 
                   labels = ws032$m004$terms$anatomicalSystems$label,
                   limits = rev(ws032$m004$terms$anatomicalSystems$node)) +
  scale_fill_manual(values = c("#CCCCCC", "#5382B2")) +
  ylim(0, 140) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), legend.position = "none")


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 5, height = 5.5)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - n_experiment across different cell lines, fill = nd_tf, simple bar

pData <- list()

pData$outputFile <- paste0(ws032$m004$figurePath, "07.png")

# first address coverage

pData$targetExperiments <- ws032$m004$tissueExperiments$cellLine %>% filter(node %in% ws032$m004$terms$cellLineCategories$node) %>% unique() %>% 
  dplyr::select(record_id, node, label)

pData$otherExperiments <- ws032$m004$tissueExperiments$cellLine %>% filter(!(record_id %in% pData$targetExperiments$record_id)) %>% unique()
# pData$otherExperiments %>% group_by(node) %>% summarize(label = first(label), n_experiment = n()) %>% ungroup() %>% arrange(n_experiment) 

pData$otherExperiments <- pData$otherExperiments %>% 
  dplyr::select(record_id) %>% 
  unique() %>% 
  mutate(node = "other", label = "other")

# attach "others"
pData$targetExperiments <- pData$targetExperiments %>% 
  rbind(pData$otherExperiments)

pData$targetExperiments <- pData$targetExperiments %>% 
  left_join(ws032$data$s03.tsv %>% dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, tf_gene_h = TF_Entrez_ID_Human, age = Age), by = "record_id") %>% 
  dplyr::select(record_id, tri_id, tf_gene_h, label, node, age) %>% 
  mutate(nd_tf = (tf_gene_h %in% ws032$m004$neuroTfs)) %>% 
  unique()

# filter for stuff found in other dbs as well ***
pData$targetExperiments <- pData$targetExperiments %>% 
  filter(tri_id %in% (ws032$data$s09.tsv %>% session$dataWrangler$extractColumn("DTRI_ID"))) %>% 
  dplyr::select(-tri_id) %>% 
  unique()
# ===


pData$a <- pData$targetExperiments %>% 
  group_by(node, nd_tf) %>%
  summarize(label = dplyr::first(label), n_experiment = n()) %>% 
  arrange(desc(n_experiment)) %>% 
  ungroup()

pData$b <- pData$a %>% 
  group_by(node) %>% 
  summarize(n_experiment = sum(n_experiment))

pData$a %>% 
  session$graphingUtils$ggplot(aes(x = node, y = n_experiment)) +
  geom_bar(stat = "identity", aes(fill = nd_tf)) +
  geom_text(data = pData$b, aes(label = n_experiment), hjust = -0.5, size = 6) +
  coord_flip() +
  scale_x_discrete(breaks = ws032$m004$terms$cellLineCategories$node, 
                   labels = ws032$m004$terms$cellLineCategories$label,
                   limits = rev(ws032$m004$terms$cellLineCategories$node)) +
  scale_fill_manual(values = c("#CCCCCC", "#5382B2")) +
  ylim(0, 250) + 
  xlab("") +
  ylab("")  +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), legend.position = "none")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6, height = 5.5)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - bubble chart, pairwise between all pairs of dbs, for pmids

pData <- list()

pData$outputFile <- paste0(ws032$m004$figurePath, "08.png")

pData$a <- ws032$m004$dtris %>%
  dplyr::select(database = Database, tri_id = DTRI_ID) %>%
  unique()

pData$b <- pData$a$database %>% unique()
names(pData$b) <- pData$b

pData$c <- do.call(rbind, pData$b %>% session$collectionUtils$lapply(function(currDb) {
  do.call(rbind, pData$b %>% session$collectionUtils$lapply(function(targetDb) {
    if (currDb == targetDb) {
      return(tibble(curr_db = currDb, target_db = targetDb, n_intersect = NA, frac_intersect = NA))
    }
    currSet <- pData$a %>% filter(database == currDb) %>% session$dataWrangler$extractColumn("tri_id") %>% unique()
    targetSet <- pData$a %>% filter(database == targetDb) %>% session$dataWrangler$extractColumn("tri_id") %>% unique()
    nIntersect <- currSet %>% intersect(targetSet) %>% length()
    fracTarget <- nIntersect / length(targetSet)
    return(tibble(curr_db = currDb, target_db = targetDb, n_intersect = nIntersect, frac_intersect = round(fracTarget, 2)))
  }, reporter = pData$outputFile))
}, reporter = pData$outputFile)) 

pData$d <- pData$a %>% 
  group_by(database) %>% summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  session$dataWrangler$extractColumn("database")

pData$d <- c("Current", pData$d[!grepl("Current", pData$d)])

pData$c %>% 
  session$graphingUtils$ggplot(aes(x = target_db, y = curr_db)) +
  geom_point(aes(size = frac_intersect)) +
  geom_text(data = pData$c %>% filter(frac_intersect >= 0.05), aes(label = frac_intersect), vjust = -1, hjust = 0.2, size = 5) +
  scale_x_discrete(limits = pData$d)  +
  scale_y_discrete(limits = pData$d)  +
  session$graphingUtils$tiltX(angle = 90) +
  ylab(label = "Base") + 
  xlab(label = "Target") +
  theme(legend.position = "none", axis.text = element_text(size = 18))


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 8, height = 6)












