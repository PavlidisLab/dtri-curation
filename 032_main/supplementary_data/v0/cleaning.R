# this is the script to clean the supplementary data & generate the final set of datasets
# worry about the numbers & figures that go into the paper later; the important thing here is making sure that everything is consistent

# curated stuff - s03, s04, s05, s06, s07
# external records - s08
# external stuff + curated stuff > dtris - s09
# curated stuff + external stuff > candidate papers - s01
# candidate papers + lambert > tfs - s02
# 
# no mesh counts > s10
# dea > s11
# negatives > s12

# candidate papers are from external papers + pubmed search + curated papers; if a paper is curated for at least one TF, then it cannot be unexamined
# tfs come from TFs + genes with least one candidate paper


ws032$cleaningV0 <- list()

ws032$cleaningV0$compilation <- readRDS(paste0(ws032$workspaceDir, "misc/", "COMPILATION_2020_0915.rds"))

ws032$cleaningV0$inputDir <- paste0(ws032$suppDataDir, "v0/")

ws032$cleaningV0$outputDir <- paste0(ws032$suppDataDir, "vfinal/")

ws032$cleaningV0$inputFiles <- ws032$cleaningV0$inputDir %>% list.files() %>% session$dataWrangler$attachNames()
ws032$cleaningV0$inputFiles <- ws032$cleaningV0$inputFiles[grepl("*.tsv", ws032$cleaningV0$inputFiles)]
ws032$cleaningV0$inputFiles <- ws032$cleaningV0$inputFiles %>% 
  session$collectionUtils$lapply(function(currTsv) { session$dataWrangler$readTsv(paste0(ws032$cleaningV0$inputDir, currTsv)) })

# put together supplementary tables one by one

ws032$cleaningV0$output <- list()


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 3 - curated records, by experiment

pData <- list()

pData$currTibble <-  ws032$cleaningV0$inputFiles$s03.tsv

pData$currTibble <- pData$currTibble %>% 
  dplyr::select(-DTRI_ID, -TF_Symbol_Human, -Target_Symbol_Human) %>% 
  dplyr::select(Experiment_ID = EXPERIMENT_ID, everything()) %>% 
  unique()

# make the table for DTRIs, then map it back

pData$dtris <- pData$currTibble %>% 
  dplyr::select(TF_Entrez_ID, Target_Entrez_ID) %>% 
  unique() %>% 
  left_join(ws032$geneOrthologs %>% dplyr::select(TF_Entrez_ID = entrez, TF_Entrez_ID_Human = entrez_h) %>% unique(), by = "TF_Entrez_ID") %>% 
  left_join(ws032$geneOrthologs %>% dplyr::select(Target_Entrez_ID = entrez, Target_Entrez_ID_Human = entrez_h) %>% unique(), by = "Target_Entrez_ID") %>% 
  unique() 

pData$dtriHuman <- pData$dtris %>% 
  dplyr::select(TF_Entrez_ID_Human, Target_Entrez_ID_Human) %>% 
  unique()

pData$dtriHuman <- pData$dtriHuman %>% 
  arrange(TF_Entrez_ID_Human, Target_Entrez_ID_Human) %>% 
  mutate(DTRI_ID = paste0("tri_", formatC(1:nrow(pData$dtriHuman), digits = 3, flag = "0"))) %>% 
  dplyr::select(DTRI_ID, everything())

pData$dtriHuman <- pData$dtriHuman %>%
  left_join(ws032$genes %>% dplyr::select(TF_Entrez_ID_Human = entrez, TF_Symbol_Human = gene), by = "TF_Entrez_ID_Human") %>% 
  left_join(ws032$genes %>% dplyr::select(Target_Entrez_ID_Human = entrez, Target_Symbol_Human = gene), by = "Target_Entrez_ID_Human") %>% 
  unique() 

pData$dtris <- pData$dtris %>% left_join(pData$dtriHuman)


# now map it back over to the experiments

pData$currTibble <- pData$currTibble %>% left_join(pData$dtris) %>% unique() 

pData$currTibble <- pData$currTibble %>% 
  dplyr::select(DTRI_ID, 
                Experiment_ID,
                TF_Symbol_Human, 
                Target_Symbol_Human, 
                TF_Entrez_ID_Human, 
                Target_Entrez_ID_Human, 
                TF_Entrez_ID, 
                Target_Entrez_ID, 
                everything())


# COMMIT IT!
ws032$cleaningV0$output$s03.tsv <- pData$currTibble


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 4 - curated records, by dtris

pData <- list()

pData$currTibble <- ws032$cleaningV0$inputFiles$s04.tsv

pData$experiments <- ws032$cleaningV0$output$s03.tsv 

pData$main <- pData$experiments %>% 
  mutate(species_str = paste0(TF_Species, ",", Target_Species)) %>% 
  mutate(Species = species_str %>% sapply(function(currSpecies) {
    mouse <- grepl("mouse", currSpecies)
    human <- grepl("human", currSpecies)
    if (mouse & human) { "mixed" }
    else if (mouse) { "mouse" }
    else if (human) { "human" }
    else { currSpecies }
  })) %>% dplyr::select(-TF_Species, -Target_Species, -species_str)

pData$cnsTerms <- session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$UBERON) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$CL)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$UBERON)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$CL)) %>% 
  unique()

pData$main <- pData$main %>% 
  mutate(CNS = Cell_Type %in% pData$cnsTerms) 

pData$attributes <- pData$main %>% 
  group_by(DTRI_ID) %>% 
  summarize(Species = paste0(Species, collapse = ","), 
            Experiment_Type = paste0(Experiment_Type, collapse = ","), 
            Context_Type = paste0(Context_Type, collapse = ","), 
            CNS = any(CNS),
            Age = paste0(Age, collapse = ","), 
            TFBS = paste0(TFBS, collapse = ","), 
            Mode = paste0(Mode, collapse = ",")) %>% 
  mutate(
    Species = Species %>% sapply(function(currSpecies) {
      human = grepl("human", currSpecies) 
      mouse = grepl("mouse", currSpecies)
      mixed = grepl("mixed", currSpecies)
      if (mixed) { "both" }
      else if (human & !mouse) { "human" }
      else if (mouse & !human) { "mouse" }
      else { "both" }
    }), 
    Binding = grepl("binding", Experiment_Type), 
    Perturbation = grepl("perturbation", Experiment_Type), 
    Reporter =  grepl("reporter", Experiment_Type), 
    Primary =  grepl("primary", Context_Type), 
    Cell_Line =  grepl("cell_line", Context_Type), 
    Prenatal = grepl("E", Age), 
    Postnatal = grepl("P", Age)
  ) %>% 
  mutate(Mode = Mode %>% sapply(function(currMode) {
    hasActivation = grepl("activation", currMode)
    hasRepression = grepl("repression", currMode)
    if (hasActivation & hasRepression) { "both" }
    else if (hasActivation) { "activation" }
    else if (hasRepression) { "repression" }
    else { "unknown" }
  })) %>% 
  mutate(TFBS = TFBS %>% sapply(function(currPosition) {
    hasProximal = grepl("proximal", currPosition)
    hasDistal = grepl("distal", currPosition)
    if (hasProximal & hasDistal) { "both" }
    else if (hasProximal) { "proximal" }
    else if (hasDistal) { "distal" }
    else { "unknown" }
  })) %>% 
  dplyr::select(-Experiment_Type, -Context_Type, -Age)

pData$main <- pData$main %>% 
  group_by(DTRI_ID) %>% 
  summarize(
    TF_Symbol_Human = dplyr::first(TF_Symbol_Human), 
    Target_Symbol_Human = dplyr::first(Target_Symbol_Human), 
    TF_Entrez_ID_Human = dplyr::first(TF_Entrez_ID_Human),
    Target_Entrez_ID_Human = dplyr::first(Target_Entrez_ID_Human),
    Experiments = paste0(unique(Experiment_ID), collapse = ", ")
  ) %>% 
  left_join(pData$attributes, by = "DTRI_ID")

pData$main <- pData$main %>% 
  group_by(TF_Symbol_Human) %>% 
  mutate(x = n(), y = nchar(Experiments)) %>% 
  arrange(desc(x), desc(y), TF_Symbol_Human, Target_Symbol_Human) %>% 
  ungroup() %>% 
  dplyr::select(DTRI_ID, 
                TF_Symbol_Human, Target_Symbol_Human, 
                TF_Entrez_ID_Human, Target_Entrez_ID_Human, 
                Species, Binding, Perturbation, Reporter, Primary, CNS, Cell_Line, Prenatal, Postnatal, TFBS, Mode, Experiments)


# COMMIT 
ws032$cleaningV0$output$s04.tsv <- pData$main



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 5 - perturbation experiment details

pData <- list()

pData$currTibble <- ws032$cleaningV0$inputFiles$s05.tsv

pData$experiments <- ws032$cleaningV0$output$s03.tsv 

pData$currTibble <- pData$currTibble %>% 
  dplyr::select(Experiment_ID, Experiment_Type, Effect, Induced) %>% 
  unique() %>% 
  arrange(Experiment_ID)

# COMMIT 
ws032$cleaningV0$output$s05.tsv <- pData$currTibble


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 6 - binding experiment details

pData <- list()

pData$currTibble <- ws032$cleaningV0$inputFiles$s06.tsv

pData$experiments <- ws032$cleaningV0$output$s03.tsv 

pData$currTibble <- pData$currTibble %>% 
  dplyr::select(Experiment_ID, Experiment_Type, Method, TF_Source_Type) %>% 
  unique() %>% 
  arrange(Experiment_ID)

# COMMIT 
ws032$cleaningV0$output$s06.tsv <- pData$currTibble


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 7 - reporter experiment details

pData <- list()

pData$currTibble <- ws032$cleaningV0$inputFiles$s07.tsv

pData$experiments <- ws032$cleaningV0$output$s03.tsv 

pData$currTibble <- pData$currTibble %>% 
  dplyr::select(Experiment_ID, Experiment_Type, Mutated, Binding_Verified) %>% 
  unique() %>% 
  arrange(Experiment_ID)

# COMMIT 
ws032$cleaningV0$output$s07.tsv <- pData$currTibble



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 8 - imported regulatory interaction records
# this one is derived from ws032$cleaningV0$compilation$records

pData <- list()

pData$currTibble <- ws032$cleaningV0$compilation$records %>% 
  dplyr::select(TF_Entrez_ID = tf_entrez, 
                Target_Entrez_ID = target_entrez, 
                TF_Species = tf_species, 
                Target_Species = target_species,
                PubMed_ID = pmid, 
                Mode = mode, 
                Database = database) %>% 
  unique()

pData$currTibble <- pData$currTibble %>% 
  mutate(Database = Database %>% sapply(function(currStr) {
    if (currStr == "cytreg") { "CytReg" }
    else if (currStr == "endb") { "ENdb" }
    else if (currStr == "htridb") { "HTRIdb_LC" }
    else if (currStr == "innatedb") { "InnateDB" }
    else if (currStr == "oreganno_lt") { "ORegAnno_LC" }
    else if (currStr == "trrust") { "TRRUST" }
    else if (currStr == "TFe") { "TFe" }
    else if (currStr == "tfacts") { "TFactS" }
    else if (currStr == "current") { "Current" }
    else { currStr }
  })) %>% 
  filter(Database != "Current") %>% 
  unique() %>% 
  na.omit()

# first deduplicate

# one record == one combination of 1. TF_Entrez_ID, 2. Target_Entrez_ID, 3. PubMed_ID
pData$currTibble <- pData$currTibble %>% 
  dplyr::select(TF_Entrez_ID, Target_Entrez_ID, PubMed_ID, Database, Mode) %>% 
  unique() %>% 
  group_by(TF_Entrez_ID, Target_Entrez_ID, PubMed_ID, Database) %>% 
  summarize(Mode = paste0(Mode, collapse = ",")) %>% 
  mutate(Mode = Mode %>% sapply(function(currMode) {
    hasActivation <- grepl("activation", currMode)
    hasRepression <- grepl("repression", currMode)
    if (hasActivation & hasRepression) { "both" }
    else if (hasActivation) { "activation" }
    else if (hasRepression) { "repression" }
    else { "unknown" }
  })) %>% 
  ungroup()


pData$currTibble <- pData$currTibble %>% 
  left_join(ws032$genes %>% dplyr::select(TF_Entrez_ID = entrez, tf_species = species) %>% unique()) %>% 
  left_join(ws032$genes %>% dplyr::select(Target_Entrez_ID = entrez, target_species = species) %>% unique())

# make the table for DTRIs, then map it back

pData$dtris <- pData$currTibble %>% 
  dplyr::select(TF_Entrez_ID, Target_Entrez_ID) %>% 
  unique() %>% 
  left_join(ws032$geneOrthologs %>% dplyr::select(TF_Entrez_ID = entrez, TF_Entrez_ID_Human = entrez_h) %>% unique(), by = "TF_Entrez_ID") %>% 
  left_join(ws032$geneOrthologs %>% dplyr::select(Target_Entrez_ID = entrez, Target_Entrez_ID_Human = entrez_h) %>% unique(), by = "Target_Entrez_ID") %>% 
  unique() 

pData$dtriHuman <- pData$dtris %>% 
  dplyr::select(TF_Entrez_ID_Human, Target_Entrez_ID_Human) %>% 
  unique()

pData$dtriHuman <- pData$dtriHuman %>% 
  left_join(ws032$cleaningV0$output$s04.tsv %>% 
              dplyr::select(DTRI_ID.x = DTRI_ID, TF_Entrez_ID_Human, Target_Entrez_ID_Human)) %>% 
  mutate(DTRI_ID.y = paste0("etri_", formatC(1:nrow(pData$dtriHuman), digits = 4, flag = "0"))) %>% 
  session$dataWrangler$mergeColumnsXy(colName = "DTRI_ID")

pData$dtriHuman <- pData$dtriHuman %>%
  left_join(ws032$genes %>% dplyr::select(TF_Entrez_ID_Human = entrez, TF_Symbol_Human = gene), by = "TF_Entrez_ID_Human") %>% 
  left_join(ws032$genes %>% dplyr::select(Target_Entrez_ID_Human = entrez, Target_Symbol_Human = gene), by = "Target_Entrez_ID_Human") %>% 
  unique() 

pData$dtris <- pData$dtris %>% left_join(pData$dtriHuman)


# now map it back over to the experiments

pData$currTibble <- pData$currTibble %>% left_join(pData$dtris) %>% unique() 

pData$currTibble <- pData$currTibble %>% 
  dplyr::select(DTRI_ID, 
                TF_Symbol_Human, 
                Target_Symbol_Human, 
                TF_Entrez_ID_Human, 
                Target_Entrez_ID_Human, 
                TF_Entrez_ID, 
                Target_Entrez_ID, 
                everything())

# COMMIT
ws032$cleaningV0$output$s09.tsv <- pData$currTibble


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 8 - master DTRI table, including current

pData <- list()

pData$external <- ws032$cleaningV0$output$s09.tsv

pData$curated <- ws032$cleaningV0$output$s04.tsv

pData$external <- pData$external %>% 
  dplyr::select(DTRI_ID, TF_Symbol_Human, Target_Symbol_Human, TF_Entrez_ID_Human, Target_Entrez_ID_Human, Database) %>% 
  unique()

pData$curated <- pData$curated %>% 
  dplyr::select(DTRI_ID, TF_Symbol_Human, Target_Symbol_Human, TF_Entrez_ID_Human, Target_Entrez_ID_Human) %>% 
  mutate(Database = "Current") %>% 
  unique()

pData$currTibble <- pData$curated %>% rbind(pData$external)

pData$currTibble <- pData$currTibble %>% 
  group_by(DTRI_ID, TF_Symbol_Human, Target_Symbol_Human, TF_Entrez_ID_Human, Target_Entrez_ID_Human) %>% 
  summarize(Databases = paste0(sort(unique(Database)), collapse = ", ")) %>% 
  ungroup() %>% 
  arrange(Databases, TF_Symbol_Human, Target_Symbol_Human) 

# COMMIT
ws032$cleaningV0$output$s08.tsv <- pData$currTibble 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 2 - TFs (lambert + anything curated)

pData <- list()

pData$dtris <- ws032$cleaningV0$output$s08.tsv

pData$curatedTfs <- pData$dtris$TF_Entrez_ID_Human %>% unique()

pData$lambertTfs <- ws032$genes %>% filter(lambert_tf == "TRUE") %>% session$dataWrangler$extractColumn("entrez_h") %>% unique()

pData$tfs <- pData$curatedTfs %>% c(pData$lambertTfs) %>% unique()

pData$currTibble <- ws032$genes %>% 
  dplyr::select(entrez, gene, species, lambert_tf, sfari_score, sfari_syndromic, go_neurodev) %>% 
  filter(entrez %in% pData$tfs) %>% 
  unique()

pData$currTibble <- pData$currTibble %>%
  dplyr::select(TF_Symbol = gene, TF_Entrez_ID = entrez, Species = species, Lambert = lambert_tf, 
                SFARI_Score = sfari_score, SFARI_Syndomic = sfari_syndromic, GO_Neurodev = go_neurodev)

pData$currTibble <- pData$currTibble %>% 
  mutate(GO_Neurodev = as.logical(GO_Neurodev), 
         SFARI_Score = as.numeric(SFARI_Score), 
         SFARI_Syndomic = as.numeric(SFARI_Syndomic)) %>% 
  mutate(Neurodev_TF = (SFARI_Score > 0 | SFARI_Syndomic > 0 | GO_Neurodev))

# COMMIT - run this supplementary before running the rest in this section...
ws032$cleaningV0$output$s02.tsv <- pData$currTibble 


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 1 - candidate papers

pData <- list()

pData$original <- ws032$cleaningV0$inputFiles$s01.tsv %>% unique()

pData$discarded <- pData$original %>% 
  filter(Status == "discarded") %>% 
  filter(TF_Symbol != "not_applicable") %>% 
  dplyr::select(TF_Symbol, PubMed_ID) %>% 
  unique() %>% 
  left_join(ws032$genes %>% dplyr::select(TF_Symbol = gene, TF_Entrez_ID = entrez_h) %>% unique(), by = "TF_Symbol") %>% 
  filter(TF_Entrez_ID != "100131938") %>% # remove expired entrez id for PBX1
  dplyr::select(TF_Entrez_ID, PubMed_ID) %>% 
  mutate(Status = "discarded")

pData$curated <- ws032$cleaningV0$output$s03.tsv %>% 
  dplyr::select(TF_Entrez_ID = TF_Entrez_ID_Human, PubMed_ID) %>% 
  mutate(Status = "curated") %>% 
  unique()

pData$pubmedQueried <- pData$original %>% 
  filter(Source == "PubMedQuery") %>% 
  filter(TF_Symbol != "not_applicable") %>% 
  dplyr::select(TF_Symbol, PubMed_ID) %>% 
  unique() %>% 
  left_join(ws032$genes %>% dplyr::select(TF_Symbol = gene, TF_Entrez_ID = entrez_h) %>% unique(), by = "TF_Symbol") %>% 
  filter(TF_Entrez_ID != "100131938") %>% # remove expired entrez id for PBX1
  dplyr::select(TF_Entrez_ID, PubMed_ID) %>% 
  mutate(Database = "PubMedQuery")

pData$externals <- ws032$cleaningV0$output$s09.tsv %>% 
  dplyr::select(TF_Entrez_ID = TF_Entrez_ID_Human, PubMed_ID, Database) %>% 
  filter(TF_Entrez_ID %in% ws032$cleaningV0$output$s02.tsv$TF_Entrez_ID) %>% 
  unique()

pData$currTibble <- pData$externals %>% 
  rbind(pData$pubmedQueried) %>% 
  full_join(pData$curated, by = c("TF_Entrez_ID", "PubMed_ID")) %>% 
  full_join(pData$discarded, by = c("TF_Entrez_ID", "PubMed_ID")) %>% 
  session$dataWrangler$mergeColumnsXy("Status") %>% 
  session$dataWrangler$fillNa("Database", "PubMedQuery") %>% 
  session$dataWrangler$fillNa("Status", "unexamined") 

pData$currTibble <- pData$currTibble %>% 
  left_join(ws032$genes %>% dplyr::select(TF_Entrez_ID = entrez, TF_Symbol = gene) %>% unique()) %>% 
  dplyr::select(TF_Symbol, everything()) %>% 
  arrange((Status), TF_Symbol, Database) %>% 
  unique()

# COMMIT
ws032$cleaningV0$output$s01.tsv <- pData$currTibble 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 2 - modify s02!!!

# COMMIT - run this supplementary before running the rest in this section...

pData <- list()

pData$currTibble <- ws032$cleaningV0$output$s02.tsv

# use this set of TFs to come up with the list of candidate papers, and then attach the number of candidate papers back here!

pData$candidates <- ws032$cleaningV0$output$s01.tsv # this is available once the next section is run (s01 and s02 are inter-dependent)
pData$candidates <- pData$candidates %>% 
  dplyr::select(TF_Entrez_ID, PubMed_ID) %>% 
  unique() %>% 
  group_by(TF_Entrez_ID) %>% 
  summarize(N_Candidate_Paper = n())

pData$currTibble <- pData$currTibble %>% 
  left_join(pData$candidates, by = "TF_Entrez_ID") %>% 
  session$dataWrangler$fillNa("N_Candidate_Paper", 0) %>% 
  arrange(dplyr::desc(Neurodev_TF), dplyr::desc(N_Candidate_Paper), TF_Symbol)


# final step: add in N PubMed Paper from query with no mesh terms

pData$allPapers <- ws032$cleaningV0$compilation$pmids %>% 
  filter(source == "nomesh") %>% 
  dplyr::select(tf_gene_h, pmid) %>% 
  unique() %>% 
  group_by(tf_gene_h) %>% 
  summarize(N_Paper_Total = n()) %>% 
  arrange(desc(N_Paper_Total)) %>% 
  filter(N_Paper_Total < 10000) %>% # filter out the outliers as the queries contain too many false positives (e.g. "REST" includes many papers irrelevant to the TF)
  dplyr::select(TF_Symbol = tf_gene_h, everything())

pData$currTibble <- pData$currTibble %>%
  left_join(pData$allPapers, by = "TF_Symbol")


# COMMIT AGAIN
ws032$cleaningV0$output$s02.tsv <- pData$currTibble 


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 10 - walcher et al. differential expression

pData <- list()

pData$expressionDataFile <- paste0(ws032$workspaceDir, "misc/", "walcher_2013/", "6128_GSE35260_expmat.unfilt.data.txt")
pData$designFile <- paste0(ws032$workspaceDir, "misc/", "walcher_2013/", "6128_GSE35260_expdesign.data.txt")

pData$expressions <- read.table(pData$expressionDataFile, 
                                header = TRUE,
                                skip = 6, 
                                sep = "\t",
                                stringsAsFactors = FALSE, 
                                quote = "") %>% as_tibble()

pData$probes <- pData$expressions %>% dplyr::select(probe = Probe, gene = GeneSymbol, entrez_id = NCBIid)

pData$expressions <- pData$expressions %>% dplyr::select(-Sequence, -GeneSymbol, -GeneName, -GemmaId, -NCBIid)
pData$expressions <- pData$expressions %>% dplyr::select(probe = Probe, everything()) %>% 
  session$dataWrangler$setColAsRownames("probe")
names(pData$expressions) <- str_extract(names(pData$expressions), "GSE35260_Biomat_[0-9]*")

pData$design <- read.table(pData$designFile, 
                           header = TRUE,
                           skip = 9, 
                           sep = "\t",
                           stringsAsFactors = FALSE, 
                           quote = "") %>% as_tibble()

pData$design <- pData$design %>% filter(grepl("Sey", Bioassay)) %>% arrange(genotype)

pData$design <- pData$design %>% mutate(genotype = genotype %>% sapply(function(currGenotype) {
  if (grepl("wild", currGenotype)) { "wildtype" } 
  else { "pax6_sey" }
}))

pData$design <- pData$design %>% 
  mutate(bioassay = str_extract(Bioassay, "GSE35260_Biomat_[0-9]*")) %>% 
  dplyr::select(bioassay, external_id = ExternalID, genotype, batch, strain)

pData$expressions <- pData$expressions[pData$design$bioassay] # making sure the order is consistent


# DIFFERENTIAL EXPRESSION ANALYSIS CALL

pData$design$genotype <- pData$design$genotype %>% factor(levels = c("wildtype", "pax6_sey"))

pData$modelMatrix <- model.matrix(~genotype, pData$design)

pData$expressions %>% names()

pData$lmFit <- lmFit(pData$expressions, pData$modelMatrix)
pData$lmFitEb <- eBayes(pData$lmFit)

pData$result <- topTable(pData$lmFitEb, number = Inf) %>% session$dataWrangler$setRownameAsColumn("probe")

pData$result <- pData$result %>% 
  left_join(pData$probes %>% dplyr::select(probe, entrez = entrez_id) %>% unique(), by = "probe") %>% 
  dplyr::select(probe, entrez, everything()) %>% 
  na.omit()

pData$result <- pData$result %>% 
  left_join(ws032$geneOrthologs %>% dplyr::select(entrez, entrez_h), by = "entrez") %>% 
  dplyr::select(-entrez) %>% 
  dplyr::select(probe, entrez = entrez_h, everything()) %>% 
  left_join(ws032$genes %>% dplyr::select(entrez, gene) %>% unique(), by = "entrez") %>% 
  dplyr::select(-probe) %>% 
  dplyr::select(entrez, gene, everything())

pData$result <- pData$result %>% 
  filter(grepl("[0-9]*", entrez)) %>% 
  group_by(entrez) %>% 
  filter(P.Value == min(P.Value)) %>% 
  ungroup() %>% 
  unique() %>% 
  arrange(adj.P.Val)

# COMMIT
ws032$cleaningV0$output$s10.tsv <- pData$result %>% na.omit()



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# SUPPLEMENTARY TABLE 12 - negative cases 
# derived manually







# write latest version to disk

ws032$cleaningV0$output %>% session$collectionUtils$lapplyWithName(function(currName, currTsv) {
  currTsv %>% write_tsv(paste0(ws032$cleaningV0$outputDir, currName))
})




















