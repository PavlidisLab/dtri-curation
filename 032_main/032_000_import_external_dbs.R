
# in this script... 
# 1. generate the paper batches from previous curation: TFe, TRRUST, OReganno, TFactS, TRRD, TRED-LC, HTRD-LC
# ---- for each source, 3 tibbles: records, interactions, papers
# 2. generate the paper batches from pubmed


source("GlobalSession.R"); session <- initGlobalSession()

ws032 <- list()

ws032$workspaceDir <- paste0(session$WORKSPACES, "032_FIRST_PAPER_PUB_READY/") 

ws032$litSourcesDir <- paste0(ws032$workspaceDir, "misc/", "literature_sources/")

ws032$genes <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "genes_master.tsv"))

# ---- 1. TFe
# **** DEPRECATED ****
# use the API to fetch the target list for every transcription factor in TFe 

# ws012$an001$preCurationData$tfe$tfCodes <- (GET("http://www.cisreg.ca/cgi-bin/tfe/api.pl?code=all-tfids") %>% as.character() %>% strsplit("\n"))[[1]]
# ws012$an001$preCurationData$tfe$tfSpecies <- ws012$an001$preCurationData$tfe$tfCodes %>% session$collectionUtils$lapply(function(currTfCode) {
#   (GET(paste0("http://cisreg.cmmt.ubc.ca/cgi-bin/tfe/api.pl?tfid=", currTfCode, "&code=species")) %>% as.character() %>% strsplit("\t\n"))[[1]]
# })
# ws012$an001$preCurationData$tfe$tfSymbols <- ws012$an001$preCurationData$tfe$tfCodes %>% session$collectionUtils$lapply(function(currTfCode) {
#   (GET(paste0("http://cisreg.cmmt.ubc.ca/cgi-bin/tfe/api.pl?tfid=", currTfCode, "&code=symbol")) %>% as.character() %>% strsplit("\t\n"))[[1]]
# })
# 
# ws012$an001$preCurationData$tfe$tfsTibble <- tibble(regulator_gene = unlist(ws012$an001$preCurationData$tfe$tfSymbols), 
#                                                     tfe_id = ws012$an001$preCurationData$tfe$tfCodes,
#                                                     species = unlist(ws012$an001$preCurationData$tfe$tfSpecies))
# 
# 
# ws012$an001$preCurationData$tfe$regulationsTibble <- do.call(rbind, ws012$an001$preCurationData$tfe$tfsTibble$tfe_id %>% session$collectionUtils$lapply(function(currId) {
#   response <- GET(paste0("http://cisreg.cmmt.ubc.ca/cgi-bin/tfe/api.pl?tfid=", currId, "&code=targets"))
#   do.call(rbind, response %>% as.character() %>% strsplit("\n") %>% unlist() %>% 
#             session$collectionUtils$lapply(function(currLine) {
#               currLine %>% strsplit("\t") %>% unlist() %>% as.data.frame() %>% t()
#             })) %>% as_tibble() %>% mutate(tfe_id = currId)
# }))
# 
# ws012$an001$preCurationData$tfe$regulationsTibble <- ws012$an001$preCurationData$tfe$regulationsTibble %>% 
#   left_join(ws012$an001$preCurationData$tfe$tfsTibble, by = "tfe_id") %>% 
#   dplyr::select(regulator_gene, target_gene = V2, pmid = V5, species) %>% 
#   mutate(paper_source = "TFe")
# 
# ws012$an001$preCurationData$tfe$regulationsTibble <- ws012$an001$preCurationData$tfe$regulationsTibble %>% 
#   mutate(species = species %>% sapply(function(currSpecies) {
#     if (currSpecies == "Homo sapiens") { "human" } 
#     else if (currSpecies == "Mus musculus") { "mouse" }
#     else if (currSpecies == "Rattus norvegicus") { "rat" }
#     else { "unknown" }
#   }))
# 
# # write TFe to regulationsTibble object
# ws012$an001$preCurationData$regulationsTibble$tfe <- ws012$an001$preCurationData$tfe$regulationsTibble

# import from previous source


# records table
# regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, pmid, tf_species, target_species, paper_source, throughput
# throughput == how many targets per tf per pmid

ws032$m001$records$tfe <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/raw_tfe.tsv"))

ws032$m001$records$tfe <- ws032$m001$records$tfe %>% filter(species %in% c("human", "mouse"))

ws032$m001$records$tfe <- ws032$m001$records$tfe %>% dplyr::select(regulator_symbol = regulator_gene, target_symbol = target_gene, everything())

ws032$m001$records$tfe <- ws032$m001$records$tfe %>% 
  dplyr::select(regulator_symbol, target_symbol, pmid, tf_species = species, paper_source) %>% 
  left_join(ws032$genes %>% dplyr::select(target_symbol = gene, target_species = species) %>% unique(), by = "target_symbol") %>% 
  dplyr::select(regulator_symbol, target_symbol, pmid, tf_species, target_species, database = paper_source) %>% 
  unique()

ws032$m001$records$tfe <- ws032$m001$records$tfe %>% 
  mutate(regulator_gene = regulator_symbol, target_gene = target_symbol) %>% 
  left_join(ws032$genes %>% dplyr::select(tf_species = species, regulator_gene = gene, regulator_entrez = entrez, regulator_entrez_h = entrez_h, regulator_gene_h = gene_h), by = c("tf_species", "regulator_gene")) %>% 
  left_join(ws032$genes %>% dplyr::select(target_species = species, target_gene = gene, target_entrez = entrez, target_entrez_h = entrez_h, target_gene_h = gene_h), by = c("target_species", "target_gene")) %>% 
  dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, 
                pmid, tf_species, target_species, database, 
                regulator_symbol_ = regulator_symbol, target_symbol_ = target_symbol)


ws032$m001$records$tfe <- ws032$m001$records$tfe %>% unique()

# ws032$m001$records$tfe %>% write_tsv(paste0(ws032$litSourcesDir, "records_tfe.tsv"))
ws032$m001$records$tfe <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "records_tfe.tsv"))


# ---- 2. Oreganno
# download from http://www.oreganno.org/dump/


ws032$m001$raw$oreganno <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/oreganno_combined.tsv"))

ws032$m001$raw$oreganno <- ws032$m001$raw$oreganno %>% 
  mutate(species = Species %>% sapply(function(currSpecies) {
    if (currSpecies == "Homo sapiens") { "human" } 
    else if (currSpecies == "Mus musculus") { "mouse" }
    else { "unknown" }
  })) %>% 
  filter(species %in% c("human", "mouse"))

ws032$m001$records$oreganno <- ws032$m001$raw$oreganno

ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% 
  dplyr::select(regulator_symbol = Regulatory_Element_Symbol, regulator_id = Regulatory_Element_ID, target_symbol = Gene_Symbol, target_id = Gene_ID, pmid = PMID, tf_species = species, database = Dataset, type = Type) %>% 
  mutate(target_species = tf_species)

ws032$m001$records$oreganno %>% group_by(database) %>% summarize(n = n())

# get rid of miRNA stuff 
ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% filter(!grepl("miRNA", type))

ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% filter(regulator_symbol != "N/A" | regulator_id != "N/A")
ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% filter(target_symbol != "N/A" | target_id != "N/A")
ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% unique()

# ws032$m001$records$oreganno %>% saveRDS(paste0(ws032$litSourcesDir, "raw/oreganno_TEMPORARY.rds"))
# ws032$m001$records$oreganno <- readRDS(paste0(ws032$litSourcesDir, "raw/oreganno_TEMPORARY.rds")) ##### INTERMEDIATE DATA STATE

ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% 
  left_join(ws032$genes %>% dplyr::select(regulator_id = ensembl, regulator_entrez = entrez), by = "regulator_id") %>% 
  left_join(ws032$genes %>% dplyr::select(regulator_id = entrez) %>% mutate(regulator_entrez = regulator_id), by = "regulator_id") %>% 
  session$dataWrangler$mergeColumnsXy("regulator_entrez") %>% 
  left_join(ws032$genes %>% dplyr::select(regulator_symbol = gene,  regulator_entrez = entrez), by = "regulator_symbol") %>% 
  session$dataWrangler$mergeColumnsXy("regulator_entrez")

ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>%
  left_join(ws032$genes %>% dplyr::select(target_id = ensembl, target_entrez = entrez), by = "target_id") %>% 
  left_join(ws032$genes %>% dplyr::select(target_id = entrez) %>% mutate(target_entrez = target_id), by = "target_id") %>% 
  session$dataWrangler$mergeColumnsXy("target_entrez") %>% 
  left_join(ws032$genes %>% dplyr::select(target_symbol = gene, target_entrez = entrez), by = "target_symbol") %>% 
  session$dataWrangler$mergeColumnsXy("target_entrez")


# organize everything into "pazar", "nfiregulome", "oreganno_lt", "oreganno_ht" for database
ws032$m001$records$oreganno %>% group_by(database) %>% summarize(n = n())
ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% mutate(database = database %>% sapply(function(currDatabase) {
  if (currDatabase == "PAZAR") { "pazar" }
  else if (currDatabase == "NFIRegulomeDB") { "nfi_regulome_db" }
  else if (currDatabase == "NRSF/REST ChIPSeq sites") { "oreganno_ht" } 
  else { "oreganno_lt" }
}))

# add in human ortholog annotations
ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% 
  left_join(ws032$genes %>% dplyr::select(regulator_entrez = entrez, regulator_gene_h = gene_h, regulator_entrez_h = entrez_h), by = "regulator_entrez") %>% 
  left_join(ws032$genes %>% dplyr::select(target_entrez = entrez, target_gene_h = gene_h, target_entrez_h = entrez_h), by = "target_entrez") 


ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% 
  dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, 
                pmid, tf_species, target_species, database, 
                regulator_symbol_ = regulator_symbol, target_symbol_ = target_symbol, regulator_id_ = regulator_id, target_id_ = target_id, type_ = type)

ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% unique()

ws032$m001$records$oreganno <- ws032$m001$records$oreganno %>% mutate(pmid = pmid %>% sapply(function(currPmid) { if (currPmid == "N/A") NA else currPmid })) 

# ws032$m001$records$oreganno %>% write_tsv(paste0(ws032$litSourcesDir, "records_oreganno.tsv"))
ws032$m001$records$oreganno <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "records_oreganno.tsv"))

# ---- 3. TRRUST
# download from https://www.grnpedia.org/trrust/downloadnetwork.php 


ws032$m001$raw$trrust <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/trrust_rawdata.human.tsv")) %>% mutate(species = "human") %>% 
  rbind(session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/trrust_rawdata.mouse.tsv")) %>% mutate(species = "mouse")) %>% 
  mutate(database = "trrust") 

ws032$m001$records$trrust <- ws032$m001$raw$trrust

ws032$m001$records$trrust <- ws032$m001$records$trrust %>% 
  dplyr::select(tf_species = species, everything()) %>% 
  mutate(target_species = tf_species) %>% 
  dplyr::select(regulator_symbol = regulator_gene, target_symbol = target_gene, everything())

ws032$m001$records$trrust <- ws032$m001$records$trrust %>% 
  left_join(ws032$genes %>% dplyr::select(regulator_symbol = gene, tf_species = species, regulator_entrez = entrez), by = c("regulator_symbol", "tf_species")) %>% 
  left_join(ws032$genes %>% dplyr::select(target_symbol = gene, target_species = species, target_entrez = entrez), by = c("target_symbol", "target_species"))

ws032$m001$records$trrust <- ws032$m001$records$trrust %>% 
  mutate(pmid = strsplit(pmid, ";")) %>% 
  unnest(cols = "pmid")

ws032$m001$records$trrust <- ws032$m001$records$trrust %>% mutate(mode = tolower(mode))


# now attach the human ortholog annotations

ws032$m001$records$trrust <- ws032$m001$records$trrust %>% 
  left_join(ws032$genes %>% dplyr::select(regulator_entrez = entrez, regulator_gene_h = gene_h, regulator_entrez_h = entrez_h), by = "regulator_entrez") %>% 
  left_join(ws032$genes %>% dplyr::select(target_entrez = entrez, target_gene_h = gene_h, target_entrez_h = entrez_h), by = "target_entrez")


# now map to standard format

ws032$m001$records$trrust <- ws032$m001$records$trrust %>% 
  dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, 
                pmid, tf_species, target_species, mode, database, 
                regulator_symbol_ = regulator_symbol, target_symbol_ = target_symbol)

ws032$m001$records$trrust <- ws032$m001$records$trrust %>% unique()


# ws032$m001$records$trrust %>% write_tsv(paste0(ws032$litSourcesDir, "records_trrust.tsv"))
ws032$m001$records$trrust <-  session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "records_trrust.tsv"))


# ---- 4. TFACTS
# download from http://www.tfacts.org/TFactS-new/TFactS-v2/index1.html 

ws032$m001$raw$tfacts <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/tfacts_signless.txt")) %>% 
  left_join(session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/tfacts_signed.txt")) %>% dplyr::select(regulator_gene, target_gene, mode), 
            by = c("regulator_gene", "target_gene"))

ws032$m001$raw$tfacts <- ws032$m001$raw$tfacts %>% session$dataWrangler$fillNa("mode", "unknown")

ws032$m001$records$tfacts <- ws032$m001$raw$tfacts

ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% 
  mutate(species = species %>% sapply(function(currSpecies) {
    hasHuman <- grepl("Homo sapien", currSpecies)
    hasMouse <- grepl("Mus musculus", currSpecies)
    if (hasHuman & hasMouse) { "human,mouse" }
    else if (hasHuman) { "human" }
    else if (hasMouse) { "mouse" }
    else { "" }
  }))

ws032$m001$records$tfacts %>% group_by(species) %>% summarize(n = n())

ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% mutate(species = strsplit(species, ",")) %>% unnest(cols = c(species))


# keep pubmed(tfacts), tred, trrd

ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% 
  mutate(database = database %>% sapply(function(currDatabase) {
    strs <- c("pubmed", "tred", "trrd")
    hasPubmed <- grepl("pubmed", tolower(currDatabase))
    hasTred <- grepl("tred", tolower(currDatabase))
    hasTrrd <- grepl("trrd", tolower(currDatabase))
    strs[c(hasPubmed, hasTred, hasTrrd)]
  })) %>% unnest(cols = c(database))

ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% 
  mutate(accession = pmid %>% sapply(function(currPmid) {
    ids <- currPmid %>% trimws() %>% strsplit(";") %>% unlist()
    ids <- ids[grepl(".+", ids)] %>% trimws()
    if (length(ids) >= 1) {
      return(ids)
    } else {
      return("unknown")
    }
  })) %>% 
  unnest() %>% 
  dplyr::select(-pmid)


ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% 
  mutate(pmid = accession %>% str_extract("^[0-9].*$")) 

ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% mutate(mode = tolower(mode))

ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% dplyr::select(regulator_symbol = regulator_gene, target_symbol = target_gene, everything())


ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% 
  left_join(ws032$genes %>% dplyr::select(regulator_symbol = gene_h, regulator_entrez = entrez, species), by = c("regulator_symbol", "species")) %>% 
  left_join(ws032$genes %>% dplyr::select(target_symbol = gene_h, target_entrez = entrez, species), by = c("target_symbol", "species"))

# attach the human annotations

ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>%
  left_join(ws032$genes %>% dplyr::select(regulator_entrez = entrez, regulator_gene_h = gene_h, regulator_entrez_h = entrez_h), by = "regulator_entrez") %>% 
  left_join(ws032$genes %>% dplyr::select(target_entrez = entrez, target_gene_h = gene_h, target_entrez_h = entrez_h), by = "target_entrez")


# standard form

ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% 
  mutate(target_species = species) %>% 
  dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, 
                pmid, tf_species = species, target_species, mode, database, 
                regulator_symbol_ = regulator_symbol, target_symbol_ = target_symbol, accession_ = accession)


ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% unique()

ws032$m001$records$tfacts <- ws032$m001$records$tfacts %>% mutate(database = database %>% sapply(function(currDatabase) { if (currDatabase == "pubmed") { "tfacts" } else { currDatabase } }))

# ws032$m001$records$tfacts %>% write_tsv(paste0(ws032$litSourcesDir, "records_tfacts.tsv"))
ws032$m001$records$tfacts <-  session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "records_tfacts.tsv"))



# ---- 5. HTRIDB
# download from http://www.lbbc.ibb.unesp.br/htri

ws032$m001$raw$htridb <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/HTRIdb_data.txt"))

ws032$m001$records$htridb <- ws032$m001$raw$htridb %>%
  dplyr::select(regulator_entrez = GENEID_TF, target_entrez = GENEID_TG, technique = TECHNIQUE, pmid = PUBMED_ID) %>% 
  mutate(tf_species = "human") %>% 
  mutate(target_species = tf_species) %>% 
  mutate(database = "htridb")

ws032$m001$records$htridb %>% group_by(technique) %>% summarize(n = n())

ws032$m001$records$htridb %>% 
  dplyr::select(regulator_entrez, target_entrez, pmid, technique) %>% unique() %>% 
  group_by(pmid, regulator_entrez) %>% mutate(throughput = n()) %>% ungroup() %>% 
  arrange(desc(throughput)) %>% 
  session$dataWrangler$extractColumn("throughput") %>% unique()

ws032$m001$records$htridb <- ws032$m001$records$htridb %>% 
  filter(technique %in% c("Chromatin Immunoprecipitation", "Electrophoretic Mobility Shift Assay"))


ws032$m001$records$htridb <- ws032$m001$records$htridb %>%
  left_join(ws032$genes %>% dplyr::select(regulator_entrez = entrez, regulator_gene_h = gene_h, regulator_entrez_h = entrez_h), by = "regulator_entrez") %>% 
  left_join(ws032$genes %>% dplyr::select(target_entrez = entrez, target_gene_h = gene_h, target_entrez_h = entrez_h), by = "target_entrez")

ws032$m001$records$htridb <- ws032$m001$records$htridb %>% 
  dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, 
                pmid, tf_species, target_species, database, 
                technique_ = technique)

ws032$m001$records$htridb %>% write_tsv(paste0(ws032$litSourcesDir, "records_htridb.tsv"))
# ws032$m001$records$htridb <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "records_htridb.tsv"))

# ---- 6. TRED (RegNetwork)
# download from http://www.regnetworkweb.org/
# join these with the TRED from TFACTS, there might be some duplicates... 


ws032$m001$raw$tred  <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/tred_regnetwork_human.csv")) %>% mutate(species = "human") %>% 
  rbind(session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/tred_regnetwork_mouse.csv")) %>% mutate(species = "mouse")) 

ws032$m001$records$tred <- ws032$m001$raw$tred %>% 
  dplyr::select(regulator_entrez = regulator_id, target_entrez = target_id, everything()) %>% 
  mutate(database = "tred_regnet")

ws032$m001$records$tred <- ws032$m001$records$tred %>%
  left_join(ws032$genes %>% dplyr::select(regulator_entrez = entrez, regulator_gene_h = gene_h, regulator_entrez_h = entrez_h), by = "regulator_entrez") %>% 
  left_join(ws032$genes %>% dplyr::select(target_entrez = entrez, target_gene_h = gene_h, target_entrez_h = entrez_h), by = "target_entrez")

ws032$m001$records$tred <- ws032$m001$records$tred %>% 
  mutate(tf_species = species, pmid = NA) %>% 
  dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, 
                pmid, tf_species, target_species = tf_species, database, 
                regulator_symbol_ = regulator_symbol, target_symbol_ = target_symbol, evidence_ = evidence, confidence_ = confidence)

ws032$m001$records$tred <- ws032$m001$records$tred %>% unique()


# ws032$m001$records$tred %>% write_tsv(paste0(ws032$litSourcesDir, "records_tred.tsv"))
ws032$m001$records$tred <-  session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "records_tred.tsv"))


# ---- 7. CytReg
# download from their paper supplementary


ws032$m001$raw$cytreg  <- session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "raw/cytreg_raw.tsv"))

ws032$m001$records$cytreg <- ws032$m001$raw$cytreg %>% 
  dplyr::select(mode = Activation.Repression, everything()) %>% 
  mutate(species = tolower(species), 
         mode = tolower(mode))

ws032$m001$records$cytreg <- ws032$m001$records$cytreg %>% 
  filter(species %in% c("human", "mouse"))

ws032$m001$records$cytreg <- ws032$m001$records$cytreg %>% 
  dplyr::select(regulator_symbol = TF, target_symbol = Cytokine, everything()) %>% 
  mutate(database = "cytreg")

ws032$m001$records$cytreg <- ws032$m001$records$cytreg %>%
  left_join(ws032$genes %>% dplyr::select(regulator_symbol = gene, regulator_entrez = entrez, regulator_gene_h = gene_h, regulator_entrez_h = entrez_h), by = "regulator_symbol") %>% 
  left_join(ws032$genes %>% dplyr::select(target_symbol = gene, target_entrez = entrez, target_gene_h = gene_h, target_entrez_h = entrez_h), by = "target_symbol")

ws032$m001$records$cytreg <- ws032$m001$records$cytreg %>% 
  mutate(tf_species = species) %>% 
  dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, 
                pmid = PMIDs, tf_species, target_species = species, mode, database, 
                regulator_symbol_ = regulator_symbol, target_symbol_ = target_symbol, assay_type_ = Assay.type, year_of_publication_ = Year.of.publication) # - fill "unknown" in place of empty string for mode... 

ws032$m001$records$cytreg <- ws032$m001$records$cytreg %>% session$dataWrangler$fillEmpty("mode", "unknown") 

ws032$m001$records$cytreg <- ws032$m001$records$cytreg %>% unique()


# ws032$m001$records$cytreg %>% write_tsv(paste0(ws032$litSourcesDir, "records_cytreg.tsv"))
ws032$m001$records$cytreg <-  session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "records_cytreg.tsv"))



# ---- 8. InnateDB
# download from their website: https://www.innatedb.com/redirect.do?go=downloadCurated


ws032$m001$raw$innatedb <- read.table(paste0(ws032$litSourcesDir, "raw/innatedb_all.txt"), sep = "\t", quote = "", colClasses = "character") %>% as_tibble() 

ws032$m001$records$innatedb <- ws032$m001$raw$innatedb %>% 
  dplyr::select(ensembl_1 = V3, ensembl_2 = V4, pmid = V9, species = V29, type_1 = V21, type_2 = V22)

ws032$m001$records$innatedb <- ws032$m001$records$innatedb %>% 
  mutate(ensembl_1 = ensembl_1 %>% str_extract("(?<=ensembl:).*"), 
         ensembl_2 = ensembl_2 %>% str_extract("(?<=ensembl:).*"), 
         pmid = pmid %>% str_extract("(?<=pubmed:).*"), 
         type_1 = type_1 %>% str_extract("MI:[0-9]*"), 
         type_2 = type_2 %>% str_extract("MI:[0-9]*")) %>% 
  mutate(species = species %>% sapply(function(currSpecies) {
    if (grepl("9606", currSpecies)) { "human" }
    else if (grepl("10090", currSpecies)) { "mouse" }
    else { currSpecies }
  }))


ws032$m001$records$innatedb <- ws032$m001$records$innatedb %>% filter(species %in% c("human", "mouse"))

ws032$m001$records$innatedb <- ws032$m001$records$innatedb <- ws032$m001$records$innatedb %>% 
  mutate(has_protein = grepl("0326", type_1) | grepl("0326", type_2), 
         has_dna = grepl("0319", type_1) | grepl("0319", type_2)) %>% 
  filter(has_protein, has_dna)


# assign temporary id
ws032$m001$records$innatedb <- ws032$m001$records$innatedb %>% 
  mutate(temp_id = paste0("temp_", formatC(1:n(), digits = floor(log10(n())), flag = "0"))) 

ws032$m001$tempData$innatedb$geneResolve <- ws032$m001$records$innatedb %>% 
  dplyr::select(temp_id, type = type_1, ensembl = ensembl_1) %>% 
  rbind(ws032$m001$records$innatedb %>% 
          dplyr::select(temp_id, type = type_2, ensembl = ensembl_2)) %>% 
  spread(key = type, value = ensembl) %>% 
  dplyr::select(temp_id, regulator_ensembl = `MI:0326`, target_ensembl = `MI:0319`)

ws032$m001$records$innatedb <- ws032$m001$records$innatedb %>% 
  left_join(ws032$m001$tempData$innatedb$geneResolve, by = "temp_id") %>% 
  dplyr::select(regulator_ensembl, target_ensembl, pmid, species) %>% 
  mutate(database = "innatedb")

##### MAP BY ENSEMBL DIRECTLY HERE. 

ws032$m001$records$innatedb <- ws032$m001$records$innatedb %>%
  left_join(ws032$genes %>% dplyr::select(regulator_ensembl = ensembl, regulator_entrez = entrez, regulator_gene_h = gene_h, regulator_entrez_h = entrez_h), by = "regulator_ensembl") %>% 
  left_join(ws032$genes %>% dplyr::select(target_ensembl = ensembl, target_entrez = entrez, target_gene_h = gene_h, target_entrez_h = entrez_h), by = "target_ensembl")

ws032$m001$records$innatedb <- ws032$m001$records$innatedb %>% mutate(mode = "unknown")

ws032$m001$records$innatedb <- ws032$m001$records$innatedb %>% 
  mutate(tf_species = species) %>% 
  dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, 
                pmid, tf_species, target_species = species, database)


ws032$m001$records$innatedb <- ws032$m001$records$innatedb %>% unique()


# ws032$m001$records$innatedb %>% write_tsv(paste0(ws032$litSourcesDir, "records_innatedb.tsv"))
ws032$m001$records$innatedb <-  session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "records_innatedb.tsv"))


# ---- 9. ENdb
# download from their website: http://www.licpathway.net/ENdb/Download.php

ws032$m001$raw$endb  <- read_tsv(paste0(ws032$litSourcesDir, "raw/endb_enhancer.txt"))

ws032$m001$records$endb <- ws032$m001$raw$endb %>% 
  dplyr::select(regulator_symbol = TF_name, target_symbol = Target_gene, pmid = Pubmed_ID, species = Species) %>% 
  mutate(species = tolower(species)) %>%
  mutate(pmid = as.character(pmid)) 

ws032$m001$records$endb <- ws032$m001$records$endb %>% 
  filter(regulator_symbol != "--", target_symbol != "--", pmid != "--") # keep only valid records

ws032$m001$records$endb <- ws032$m001$records$endb %>% 
  mutate(database = "endb")

ws032$m001$records$endb <- ws032$m001$records$endb %>% # break lists of TFs & targets into individual records
  mutate(regulator_symbol = strsplit(regulator_symbol, ",")) %>% 
  unnest(cols = "regulator_symbol") %>% 
  mutate(target_symbol = strsplit(target_symbol, ",")) %>% 
  unnest(cols = "target_symbol") 

ws032$m001$records$endb <- ws032$m001$records$endb %>% # now annotate the human orthologs
  left_join(ws032$genes %>% dplyr::select(regulator_symbol = gene, regulator_entrez = entrez, regulator_gene_h = gene_h, regulator_entrez_h = entrez_h), by = "regulator_symbol") %>% 
  left_join(ws032$genes %>% dplyr::select(target_symbol = gene, target_entrez = entrez, target_gene_h = gene_h, target_entrez_h = entrez_h), by = "target_symbol")

ws032$m001$records$endb <- ws032$m001$records$endb %>% 
  mutate(tf_species = species, mode = "unknown") %>% 
  dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, 
                pmid, tf_species, target_species = species, mode, database) # - fill "unknown" in place of empty string for mode... 

ws032$m001$records$endb <- ws032$m001$records$endb %>% session$dataWrangler$fillEmpty("mode", "unknown") 

ws032$m001$records$endb <- ws032$m001$records$endb %>% unique()


# ws032$m001$records$endb %>% write_tsv(paste0(ws032$litSourcesDir, "records_endb.tsv"))
ws032$m001$records$endb <-  session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, "records_endb.tsv"))








############################
######## RECORDS ###########


# load in all the records here

ws032$m001$records <- list()

ws032$m001$recordFiles <- list.files(ws032$litSourcesDir)
ws032$m001$recordFiles <- ws032$m001$recordFiles[grepl("records_", ws032$m001$recordFiles)]
ws032$m001$recordFiles <- ws032$m001$recordFiles[!grepl("master", ws032$m001$recordFiles)]
names(ws032$m001$recordFiles) <- ws032$m001$recordFiles %>% str_extract("(?<=records_).*(?=.tsv)")

ws032$m001$records <- ws032$m001$recordFiles %>% session$collectionUtils$lapply(function(currFile) {
  session$dataWrangler$readTsv(paste0(ws032$litSourcesDir, currFile))
}, reporter = "reading in cache")


# regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, pmid, mode, tf_species, target_species, database

ws032$m001$records %>% names()
ws032$m001$records %>% names() %>% length()

ws032$m001$master <- do.call(rbind, ws032$m001$records %>% lapply(function(currTibble) {
  cols <- c("regulator_gene_h", "target_gene_h", "regulator_entrez_h", "target_entrez_h", "regulator_entrez", "target_entrez", 
            "pmid", "mode", "tf_species", "target_species", "database", 
            "regulator_symbol_", "target_symbol_")
  missingCols <- cols %>% setdiff(names(currTibble))
  clone <- currTibble
  for (currCol in missingCols) {
    clone[[currCol]] <- "unknown"
  }
  return(clone[cols])
})) %>% unique()



# map the final gene aliases

# ALIAS MAPPINGS

ws032$genesAliases <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/genes/", "ncbi_genes_human_200707.txt")) %>% 
  dplyr::select(gene = Symbol, entrez = GeneID, alias = Aliases) %>% mutate(species = "human") %>% 
  rbind(session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/genes/", "ncbi_genes_mouse_200707.txt")) %>% 
          dplyr::select(gene = Symbol, entrez = GeneID, alias = Aliases) %>% mutate(species = "mouse"))

ws032$genesAliases <- ws032$genesAliases %>%
  mutate(alias = strsplit(alias, ",")) %>% 
  unnest(cols = "alias") %>% 
  mutate(alias = trimws(alias)) %>% 
  dplyr::select(alias, species, entrez)



# regulator aliases

ws032$m001$master <- ws032$m001$master %>% 
  filter(!session$dataWrangler$isInvalid(regulator_entrez)) %>% 
  rbind(ws032$m001$master %>% 
          filter(session$dataWrangler$isInvalid(regulator_entrez)) %>% 
          dplyr::select(-regulator_entrez, -regulator_entrez_h, -regulator_gene_h) %>% 
          left_join(ws032$genesAliases %>% dplyr::select(regulator_symbol_ = alias, tf_species = species, regulator_entrez = entrez), by = c("regulator_symbol_", "tf_species")) %>% 
          left_join(ws032$genes %>% dplyr::select(regulator_entrez = entrez, regulator_gene_h = gene_h, regulator_entrez_h = entrez_h), by = "regulator_entrez") %>% 
          session$dataWrangler$removeInvalidRows("regulator_entrez") %>% 
          dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, pmid, mode, tf_species, target_species, database, regulator_symbol_, target_symbol_))


ws032$m001$master %>% 
  filter(session$dataWrangler$isInvalid(regulator_entrez))


# target aliases

ws032$m001$master <- ws032$m001$master %>% 
  filter(!session$dataWrangler$isInvalid(target_entrez)) %>% 
  rbind(ws032$m001$master %>% 
          filter(session$dataWrangler$isInvalid(target_entrez)) %>% 
          dplyr::select(-target_entrez, -target_entrez_h, -target_gene_h) %>% 
          left_join(ws032$genesAliases %>% dplyr::select(target_symbol_ = alias, target_species = species, target_entrez = entrez), by = c("target_symbol_", "target_species")) %>% 
          left_join(ws032$genes %>% dplyr::select(target_entrez = entrez, target_gene_h = gene_h, target_entrez_h = entrez_h), by = "target_entrez") %>% 
          session$dataWrangler$removeInvalidRows("target_entrez") %>% 
          dplyr::select(regulator_gene_h, target_gene_h, regulator_entrez_h, target_entrez_h, regulator_entrez, target_entrez, pmid, mode, tf_species, target_species, database, regulator_symbol_, target_symbol_))


ws032$m001$master %>% 
  filter(session$dataWrangler$isInvalid(target_entrez))



# final step changes ++++ 

ws032$m001$master <- ws032$m001$master %>% 
  session$dataWrangler$removeInvalidRows(c("regulator_entrez", "target_entrez", "pmid")) %>% 
  dplyr::select(tf_gene_h = regulator_gene_h, target_gene_h, tf_entrez_h = regulator_entrez_h, target_entrez_h, tf_entrez = regulator_entrez, target_entrez, 
                pmid, mode, tf_species, target_species, database)


ws032$m001$master <- ws032$m001$master %>% 
  mutate(mode = mode %>% sapply(function(currMode) { 
    if (currMode %in% c("activation", "up")) { "activation" }
    else if (currMode %in% c("repression", "down")) { "repression" }
    else { currMode }
  }))


# assign ids now


ws032$m001$master <- ws032$m001$master %>% 
  group_by(database) %>% 
  mutate(record_id = paste0(database, "_", formatC(1:n(), digits = floor(log10(n())), flag = "0"))) %>% 
  ungroup()


ws032$m001$master %>% write_tsv(paste0(ws032$litSourcesDir, "MASTER.tsv"))


