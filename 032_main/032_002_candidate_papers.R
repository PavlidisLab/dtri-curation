# run 032_SETUP.R first

ws032$m002 <- list()
ws032$m002$figurePath <- paste0(ws032$figuresDir, "f02_")

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - venn diagram of number of papers provided by external databases vs by pubmed query

pData <- list()

pData$outputFile <- paste0(ws032$m002$figurePath, "01.png")

pData$candidates <- ws032$data$s01.tsv

pData$main <- pData$candidates %>% 
  dplyr::select(PubMed_ID, Database) %>% unique() %>% 
  mutate(Database = (Database == "PubMedQuery")) %>% 
  mutate(Database = Database %>% sapply(function(currDatabase) {
    if (currDatabase) { "PubMed" }
    else { "External" }
  })) %>% 
  unique()

pData$Databases <- pData$main$Database %>% unique() %>% session$dataWrangler$attachNames()
pData$Databases <- pData$Databases %>% lapply(function(currDatabase) {
  pData$main %>% 
    filter(Database == currDatabase) %>% 
    session$dataWrangler$extractColumn("PubMed_ID")
}) 

print(paste0("printing: ", pData$outputFile))
venn.diagram (
  x = pData$Databases,
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
# FIGURE - venn diagram of number of papers provided by external databases vs by pubmed query, filter = nd_tf

pData <- list()

pData$outputFile <- paste0(ws032$m002$figurePath, "02.png")

pData$candidates <- ws032$data$s01.tsv
pData$tfs <- ws032$data$s02.tsv
pData$neuroTfs <- pData$tfs %>% filter(Neurodev_TF == "TRUE") %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique()

pData$main <- pData$candidates %>%
  filter(TF_Entrez_ID %in% pData$neuroTfs) %>% 
  dplyr::select(PubMed_ID, Database) %>% unique() %>% 
  mutate(Database = (Database == "PubMedQuery")) %>% 
  mutate(Database = Database %>% sapply(function(currDatabase) {
    if (currDatabase) { "PubMed" }
    else { "External" }
  })) %>% 
  unique()

pData$Databases <- pData$main$Database %>% unique() %>% session$dataWrangler$attachNames()
pData$Databases <- pData$Databases %>% lapply(function(currDatabase) {
  pData$main %>% 
    filter(Database == currDatabase) %>% 
    session$dataWrangler$extractColumn("PubMed_ID")
}) 

print(paste0("printing: ", pData$outputFile))
venn.diagram (
  x = pData$Databases,
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
  fill = rev(c("#FDE8D9", "#DCF0E7"))
)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - number of papers per Database (without considering overlaps between Databases), simple bar, fill = nd_tf

pData <- list()

pData$outputFile <- paste0(ws032$m002$figurePath, "03.png")

pData$candidates <- ws032$data$s01.tsv
pData$tfs <- ws032$data$s02.tsv
pData$neuroTfs <- pData$tfs %>% filter(Neurodev_TF == "TRUE") %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique()

pData$pmid <- pData$candidates %>% 
  dplyr::select(PubMed_ID, TF_Entrez_ID) %>% unique() %>% 
  mutate(neurodev = TF_Entrez_ID %>% sapply(function(currTf) { currTf %in% pData$neuroTfs })) %>% 
  dplyr::select(PubMed_ID, neurodev) %>% 
  unique() %>% 
  group_by(PubMed_ID) %>% 
  summarize(neurodev_str = paste0(neurodev, collapse = ",")) %>% 
  mutate(neurodev = grepl("TRUE", neurodev_str)) %>% 
  dplyr::select(-neurodev_str)

pData$bar <- pData$candidates %>% 
  left_join(pData$pmid, by = "PubMed_ID")

pData$bar <- pData$bar %>% 
  dplyr::select(Database, neurodev, PubMed_ID) %>% 
  unique() %>% 
  group_by(Database, neurodev) %>% 
  summarize(n_paper = n()) %>% 
  ungroup() %>% 
  arrange(desc(n_paper))

pData$sum <- pData$bar %>% 
  group_by(Database) %>% 
  summarize(n_paper = sum(n_paper)) %>% 
  arrange(desc(n_paper))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = Database, y = n_paper)) +
  geom_bar(aes(fill = neurodev), stat = "identity") +
  scale_fill_manual(values = c("grey80", "#3182BC")) +
  geom_text(data = pData$sum, aes(label = n_paper), hjust = -0.2, size = 8) +
  ylim(0, 14000) + 
  coord_flip() + 
  scale_x_discrete(limits = rev(c("PubMedQuery", pData$sum$Database[!grepl("PubMedQuery", pData$sum$Database)]))) +
  ylab("Number of Candidate Papers") +
  xlab("") +
  theme(legend.position = "none", axis.text = element_text(size = 20))


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, width = 10, height = 6)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - bubble chart, pairwise between all pairs of dbs, for pmids

pData <- list()

pData$outputFile <- paste0(ws032$m002$figurePath, "04.png")

pData$candidates <- ws032$data$s01.tsv

pData$pmids <- pData$candidates %>%
  dplyr::select(Database, PubMed_ID) %>%
  unique()

pData$Databases <- pData$pmids$Database %>% unique() %>% session$dataWrangler$attachNames()

pData$main <- do.call(rbind, pData$Databases %>% session$collectionUtils$lapply(function(currDb) {
  do.call(rbind, pData$Databases %>% session$collectionUtils$lapply(function(targetDb) {
    if (currDb == targetDb) {
      return(tibble(curr_db = currDb, target_db = targetDb, n_intersect = NA, frac_intersect = NA))
    }
    currSet <- pData$pmids %>% filter(Database == currDb) %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique()
    targetSet <- pData$pmids %>% filter(Database == targetDb) %>% session$dataWrangler$extractColumn("PubMed_ID") %>% unique()
    nIntersect <- currSet %>% intersect(targetSet) %>% length()
    fracTarget <- nIntersect / length(targetSet)
    return(tibble(curr_db = currDb, target_db = targetDb, n_intersect = nIntersect, frac_intersect = round(fracTarget, 2)))
  }, reporter = pData$outputFile))
}, reporter = pData$outputFile)) 

pData$Databases <- pData$pmids %>% 
  group_by(Database) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  session$dataWrangler$extractColumn("Database")

pData$Databases <- c("PubMedQuery", pData$Databases[!grepl("PubMedQuery", pData$Databases)])

pData$main %>% 
  session$graphingUtils$ggplot(aes(x = target_db, y = curr_db)) +
  geom_point(aes(size = frac_intersect)) +
  geom_text(data = pData$main %>% filter(frac_intersect >= 0.05), aes(label = frac_intersect), vjust = -1, hjust = 0.2, size = 5) +
  scale_x_discrete(limits = pData$Databases)  +
  scale_y_discrete(limits = pData$Databases) +
  session$graphingUtils$tiltX(angle = 90) +
  ylab(label = "Base") + 
  xlab(label = "Target") +
  theme(legend.position = "none", axis.text = element_text(size = 18))

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, width = 8, height = 6)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - semi sankey showing pairing of Tfs & papers

pData <- list()

pData$outputFile <- paste0(ws032$m002$figurePath, "05.png")

pData$candidates <- ws032$data$s01.tsv
pData$tfs <- ws032$data$s02.tsv
pData$neuroTfs <- pData$tfs %>% filter(Neurodev_TF == "TRUE") %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique()

pData$mapping <- pData$candidates %>% 
  dplyr::select(TF_Entrez_ID, PubMed_ID) %>% 
  unique()

pData$paperCount <- pData$mapping %>%
  group_by(TF_Entrez_ID) %>% 
  summarize(n_paper = n()) %>% 
  arrange(desc(n_paper)) 

# 2211 "transcription factors" in the inclusive sense
pData$paperCount <- pData$paperCount %>% 
  full_join(pData$tfs %>% dplyr::select(TF_Entrez_ID), by = "TF_Entrez_ID") %>% 
  session$dataWrangler$fillNa("n_paper", 0) %>% 
  arrange(desc(n_paper)) %>% 
  unique()

pData$paperCount <- pData$paperCount %>% 
  mutate(rank_paper = 1:nrow(pData$paperCount)) %>% 
  mutate(rank_paper = rank_paper / max(rank_paper))

pData$tfsCount <- pData$mapping %>%
  group_by(PubMed_ID) %>% 
  summarize(n_tfs = n()) %>% 
  arrange(desc(n_tfs))

pData$tfsCount <- pData$tfsCount  %>% 
  mutate(rank_tf = 1:nrow(pData$tfsCount)) %>% 
  mutate(rank_tf = rank_tf / max(rank_tf))

pData$mapping <- pData$mapping %>% 
  mutate(rela_id = paste0("pr_", formatC(1:nrow(pData$mapping), digits = 5, flag = "0"))) %>% 
  left_join(pData$paperCount %>% dplyr::select(TF_Entrez_ID, rank_paper), by = "TF_Entrez_ID") %>% 
  left_join(pData$tfsCount %>% dplyr::select(PubMed_ID, rank_tf), by = "PubMed_ID")

pData$format <- function(tibble) {
  tibble %>% 
    dplyr::select(rela_id, rank_paper, rank_tf) %>% 
    reshape2::melt(id = "rela_id") %>% 
    as_tibble()
}


pData$highlight <- pData$mapping %>% filter(TF_Entrez_ID %in% c("10716")) %>% pData$format()

pData$neurodevRanks <- pData$mapping %>% filter(TF_Entrez_ID %in% pData$neuroTfs) %>% dplyr::select(TF_Entrez_ID, rank_paper) %>% unique()

pData$allRankPaper <- tibble(rank_paper = seq(0.01, 1, 0.001))

pData$mapping %>% 
  pData$format() %>% 
  session$graphingUtils$ggplot(aes(x = variable, y = value), axes = FALSE) +
  geom_line(aes(group = rela_id), alpha = 0.01) +
  geom_line(data = pData$highlight, aes(group = rela_id), color = "red", alpha = 0.6) +
  geom_rug(data = pData$allRankPaper, aes(x = "rank_paper", y = rank_paper), color = "grey80", size = 0.5) +
  geom_rug(data = pData$neurodevRanks, aes(x = "rank_paper", y = rank_paper), color = "#4682B7", size = 0.2) +
  scale_x_discrete(expand = expand_scale(add = 0.05)) + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) +
  xlab("") + ylab("")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, width = 6, height = 6)


# in-text analysis - overlap of top TFs in either pubmed or external

pData$computeCumNPaper <- function(pmids) {
  
  mapping <- pmids %>% unique()
  
  mapping <- mapping %>% 
    group_by(tf_gene_h) %>% 
    summarize(n_paper = n(), 
              pmids = list(pmid)) %>% 
    arrange(desc(n_paper))
  
  # include all the TFs here
  
  mapping <- mapping %>% 
    right_join(ws032$data$s02.tsv %>% dplyr::select(tf_gene_h = TF_Entrez_ID) %>% unique(), by = "tf_gene_h") %>% unique() %>% 
    session$dataWrangler$fillNa("n_paper", 0) 
  
  cumNPaper <- do.call(rbind, mapping$tf_gene_h %>% 
                         session$collectionUtils$lapply(function(currTf) {
                           currRecord <- mapping %>% filter(tf_gene_h == currTf)
                           currNPaper <- currRecord$n_paper
                           currRecords <- mapping %>% filter(n_paper >= currNPaper)
                           currCumNPaper <- currRecords$pmids %>% unlist() %>% unique() %>% length()
                           tibble(tf_gene_h = currTf, cum_n_paper = currCumNPaper)
                         }, "computing cumu"))
  
  mapping <- mapping %>% 
    left_join(cumNPaper, by = "tf_gene_h") %>% 
    mutate(frac_n_paper = cum_n_paper / max(cum_n_paper)) %>% 
    arrange(cum_n_paper)
  
  mapping <- mapping %>% 
    mutate(rank_n_paper = 1:nrow(mapping))
  
}

# first for pubmed: 

pData$cumuPubmed <- ws032$data$s01.tsv %>% 
  filter(Database == "PubMedQuery") %>% 
  dplyr::select(tf_gene_h = TF_Entrez_ID, pmid = PubMed_ID) %>% 
  unique() %>%
  pData$computeCumNPaper()

pData$topTfsPubmed <- pData$cumuPubmed %>% filter(frac_n_paper <= 0.9)

pData$cumuExternal <- ws032$data$s01.tsv %>% 
  filter(Database != "PubMedQuery") %>% 
  dplyr::select(tf_gene_h = TF_Entrez_ID, pmid = PubMed_ID) %>% 
  unique() %>%
  pData$computeCumNPaper()

pData$topTfsExternal <- pData$cumuExternal %>% filter(frac_n_paper <= 0.9)

(pData$topTfsPubmed$tf_gene_h %>% intersect(pData$topTfsExternal$tf_gene_h) %>% length()) / 
  length((union(pData$topTfsPubmed$tf_gene_h, pData$topTfsExternal$tf_gene_h)))


# 364 / 2235; 16% of TFs cover 90% of the papers

