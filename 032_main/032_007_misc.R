# run 032_SETUP.R first

ws032$m007 <- list()

ws032$m007$figurePath <- paste0(ws032$figuresDir, "f04_")

ws032$m007$dtris <- ws032$data$s09.tsv %>% dplyr::select(DTRI_ID, Database) %>% 
  rbind(ws032$data$s04.tsv %>% dplyr::select(DTRI_ID) %>% mutate(Database = "Current")) %>% 
  unique()

ws032$m007$neuroTfs <- ws032$data$s02.tsv %>% filter(Neurodev_TF == "TRUE") %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique()


# This file contains scripts for any ad hoc calculations for the main text

ws032$data$s04.tsv$TF_Entrez_ID_Human %>% unique() %>% length()

ws032$data$s04.tsv$Target_Entrez_ID_Human %>% unique() %>% length()


# TRRUST recovered # of pubmed in other resources? 
pData <- list()
pData$trrustPmids <- ws032$data$s01.tsv %>% filter(Database == "TRRUST") %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique()
pData$notTrrust <- ws032$data$s01.tsv %>% filter(Database != "TRRUST") %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique()

(pData$notTrrust %>% intersect(pData$trrustPmids) %>% length()) / length(pData$notTrrust)

pData$notTrrust %>% setdiff(pData$trrustPmids) %>% length()


# correlation between n candidate paper & n total paper

as.numeric(ws032$data$s02.tsv$N_Candidate_Paper) %>% 
  cor(as.numeric(ws032$data$s02.tsv$N_Paper_Total), use = "complete.obs", method = "spearman")



# breakdown of DTRIs by species


ws032$data$s04.tsv %>% 
  group_by(Species) %>% 
  summarize(n = n())


# number of TFs with more than ten targets

pData <- list()
pData$main <-  ws032$data$s04.tsv %>% 
  dplyr::select(TF_Entrez_ID_Human, Target_Entrez_ID_Human) %>% 
  unique() %>% 
  group_by(TF_Entrez_ID_Human) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n))
pData$main %>% 
  filter(n >= 10) %>% 
  session$dataWrangler$extractColumn("n") %>% 
  sum()
pData$main %>% 
  filter(n < 10) %>% 
  session$dataWrangler$extractColumn("n") %>% 
  sum()

pData$main %>% 
  filter(n >= 10) %>% 
  mutate(neurodev = (TF_Entrez_ID_Human %in% ws032$m007$neuroTfs)) %>% 
  arrange(desc(neurodev)) %>% View


# number of targets with more than ten TFs

pData <- list()
pData$main <-  ws032$data$s04.tsv %>% 
  dplyr::select(Target_Entrez_ID_Human, TF_Entrez_ID_Human) %>% 
  unique() %>% 
  group_by(Target_Entrez_ID_Human) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n))
pData$main %>% 
  filter(n >= 10) %>% 
  session$dataWrangler$extractColumn("n") %>% 
  sum()
pData$main %>% 
  filter(n < 10) %>% 
  session$dataWrangler$extractColumn("n") %>% 
  sum()

pData$main %>% 
  filter(n >= 10) %>% 
  mutate(neurodev = (Target_Entrez_ID_Human %in% ws032$m007$neuroTfs)) %>% 
  arrange(desc(neurodev)) %>% 
  left_join(ws032$genesLabs %>% dplyr::select(Target_Entrez_ID_Human = entrez, gene))



# number of experiments in cell line vs primary tissue or cells

ws032$data$s03.tsv %>% 
  dplyr::select(Experiment_ID, Context_Type) %>% 
  unique() %>% 
  group_by(Context_Type) %>% 
  summarize(n = n())


# number of DTRIs with primary context experiment

ws032$data$s04.tsv %>% 
  group_by(Primary) %>% 
  summarize(n = n())




# how many candidate papers hasn't been examined

ws032$data$s01.tsv %>% 
  filter(!(PubMed_ID %in% (ws032$data$s01.tsv %>% filter(Status != "unexamined") %>% session$dataWrangler$extractColumn("PubMed_ID")))) %>% 
  session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

ws032$data$s01.tsv %>% 
  filter(TF_Entrez_ID %in% ws032$m007$neuroTfs) %>% 
  filter(!(PubMed_ID %in% (ws032$data$s01.tsv %>% filter(Status != "unexamined") %>% session$dataWrangler$extractColumn("PubMed_ID")))) %>% 
  session$dataWrangler$extractColumn("PubMed_ID") %>% unique() %>% length()

# how many dtris per paper? 

ws032$data$s03.tsv %>% dplyr::select(PubMed_ID, DTRI_ID) %>% unique() %>% 
  group_by(PubMed_ID) %>% 
  summarize(n = n()) -> x
  
x$n %>% mean()













