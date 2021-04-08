source("GlobalSession.R"); session <- initGlobalSession()

ws032 <- list()

ws032$workspaceDir <- paste0(session$WORKSPACES, "032_FIRST_PAPER_PUB_READY/") 
ws032$suppDataDir <- paste0(ws032$workspaceDir, "supplementary_data/") 

ws032$dataDir <- paste0(ws032$suppDataDir, "vfinal/")

ws032$figuresDir <- paste0(ws032$workspaceDir, "figures/")

ws032$genes <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "genes_master.tsv"))

ws032$geneOrthologs <- ws032$genes %>% 
  dplyr::select(entrez, entrez_h) %>% 
  unique() %>% 
  filter(!(entrez %in% c("13047", "17268", "14955", "21674", "21414", "12577"))) %>% # manual mappings
  rbind(tibble(entrez = c("13047", "17268", "14955", "21674", "21414", "12577"), 
               entrez_h = c("1523", "4211", "283120", "21674", "6932", "1028"))) %>% 
  unique()

ws032$geneOrthologs <- ws032$geneOrthologs %>% 
  rbind(tibble(entrez = c("4513", "7682", "360171", "441026", "79424", "319115", "6957"), 
               entrez_h = c("4513", "7682", "360171", "441026", "79424", "319115", "6957")))

ws032$geneOrthologs <- ws032$geneOrthologs %>% unique()

# for anything that doesn't have a human ortholog, just use the mouse
ws032$geneOrthologs$entrez_h[which(ws032$geneOrthologs$entrez_h == "not_applicable")] <- ws032$geneOrthologs$entrez[which(ws032$geneOrthologs$entrez_h == "not_applicable")]

ws032$genesLabs <- ws032$genes %>% dplyr::select(entrez, gene) %>% unique()



# LOAD THE SUPPLEMENTARY DATA FOR DOWN STREAM ANALYSIS

ws032$data <- ws032$dataDir %>% 
  list.files() %>% 
  session$dataWrangler$attachNames() %>% 
  lapply(function(currFile) { session$dataWrangler$readTsv(paste0(ws032$dataDir, currFile)) })

