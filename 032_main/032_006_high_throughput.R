# run 032_SETUP.R first

ws032$m006 <- list()

ws032$m006$figurePath <- paste0(ws032$figuresDir, "f06_")

ws032$m006$s03 <- ws032$data$s03.tsv %>% 
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


ws032$m006$dea <- ws032$data$s10.tsv %>% 
  session$dataWrangler$colsToNumeric(c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"))

# dea ranking

ws032$m006$dea$rankingPvalue <- ws032$m006$dea %>% 
  mutate(rank = rank(P.Value)) %>% 
  group_by(entrez) %>% 
  summarize(rank = min(rank)) %>% 
  arrange(rank) %>% 
  session$dataWrangler$extractColumn("entrez") %>% 
  unique() %>% 
  as.character()


# first, code a system to come up with approriate targets (construct an easily queri-able tibble)
# - use the curation sheet (experiment types, context types, age, etc)
# - use the ontology tree (CNS, brain, nervous system, pancreas, etc)
# - use the network to come up with indirect targets

ws032$m006$directTargets <- ws032$m006$s03 %>% 
  dplyr::select(tri_id = DTRI_ID, 
                tf_gene_h = TF_Entrez_ID_Human, 
                target_gene_h = Target_Entrez_ID_Human, 
                species = Species,  
                experiment_type = Experiment_Type, 
                context_type = Context_Type, 
                tissue_or_cell = Cell_Type, 
                age = Age) %>% 
  filter(tf_gene_h == "5080") %>% 
  filter(target_gene_h != "5080") %>% # get rid of auto-regulation interaction
  mutate(is_embryonic = grepl("E", age)) %>% 
  unique()


ws032$m006$directTargets <- ws032$m006$directTargets %>% 
  left_join(session$ontologyUtils$getLabels(ws032$m006$directTargets$tissue_or_cell, scope = session$ontologyUtils$ontologies$UBERON) %>% 
              full_join(session$ontologyUtils$getLabels(ws032$m006$directTargets$tissue_or_cell, scope = session$ontologyUtils$ontologies$CL), by = "node") %>% 
              session$dataWrangler$mergeColumnsXy("label") %>% 
              full_join(session$ontologyUtils$getLabels(ws032$m006$directTargets$tissue_or_cell, scope = session$ontologyUtils$ontologies$CLO), by = "node") %>% 
              session$dataWrangler$mergeColumnsXy("label") %>% 
              dplyr::select(tissue_or_cell = node, label), 
            by = "tissue_or_cell")

ws032$m006$terms$cns <- session$ontologyUtils$getChildNodeRecursive(node = "UBERON_0001017", scope = session$ontologyUtils$ontologies$UBERON) %>% 
  c(session$ontologyUtils$getChildNodeRecursive(node = "UBERON_0001017", scope = session$ontologyUtils$ontologies$CL)) %>% 
  unique()

ws032$m006$terms$eye <- session$ontologyUtils$getChildNodeRecursive(node = "UBERON_0000970", scope = session$ontologyUtils$ontologies$UBERON) %>% 
  c(session$ontologyUtils$getChildNodeRecursive(node = "UBERON_0000970", scope = session$ontologyUtils$ontologies$CL)) %>% 
  unique()

ws032$m006$directTargets <- ws032$m006$directTargets %>% 
  group_by(target_gene_h) %>% 
  mutate(is_cns = tissue_or_cell %in% ws032$m006$terms$cns, 
         is_eye = tissue_or_cell %in% ws032$m006$terms$eye) %>% 
  summarize(species_str = paste0(species, collapse = ","), 
            experiment_type_str = paste0(experiment_type, collapse = ","), 
            context_type_str = paste0(context_type, collapse = ","), 
            is_cns = any(is_cns), 
            is_eye = any(is_eye), 
            is_embryonic = any(is_embryonic))


ws032$m006$directTargets <- ws032$m006$directTargets %>% 
  filter(target_gene_h %in% ws032$m006$dea$rankingPvalue)


# construct a ws027$targetSets to keep track of different sets of targets under different names

ws032$m006$targetSet <- list()

ws032$m006$targetSet$all <- ws032$m006$directTargets$target_gene_h %>% unique()

# BREAKDOWN by experiment type

ws032$m006$targetSet$expTypesMultiple <- ws032$m006$directTargets %>% 
  mutate(n_experiment_type = experiment_type_str %>% sapply(function(currStr) {
    currStr %>% strsplit(",") %>% unlist() %>% unique() %>% length()
  })) %>% 
  filter(n_experiment_type > 1) %>% 
  session$dataWrangler$extractColumn("target_gene_h") %>% 
  unique()

ws032$m006$targetSet$expTypesSingle <- ws032$m006$directTargets %>% 
  mutate(n_experiment_type = experiment_type_str %>% sapply(function(currStr) {
    currStr %>% strsplit(",") %>% unlist() %>% unique() %>% length()
  })) %>% 
  filter(n_experiment_type <= 1) %>% 
  session$dataWrangler$extractColumn("target_gene_h") %>% 
  unique()

# BREAKDOWN by tissue

ws032$m006$targetSet$tissueCns <- ws032$m006$directTargets %>% 
  filter(is_cns) %>% 
  session$dataWrangler$extractColumn("target_gene_h") %>% 
  unique()

ws032$m006$targetSet$tissueEye <- ws032$m006$directTargets %>% 
  filter(is_eye, !is_cns) %>% 
  session$dataWrangler$extractColumn("target_gene_h") %>% 
  unique()

ws032$m006$targetSet$tissueOther <- ws032$m006$directTargets %>% 
  filter(!is_eye, !is_cns) %>% 
  session$dataWrangler$extractColumn("target_gene_h") %>% unique()


ws032$m006$targetSetGenes <- ws032$m006$targetSet %>% lapply(function(currSet) {
  ws032$genesLabs %>% filter(entrez %in% currSet)
})


# process stuff that takes time here, before building more visuals

# write functions to help with hypergeometric tests, ROC, boodstrapping, etc

# ws032$m006$rocPValue <- do.call(rbind, ws032$m006$targetSet %>% session$collectionUtils$lapplyWithName(function(currName, currTargets) {
#   session$evaluationUtils$roc(ws032$m006$dea$rankingPvalue, currTargets) %>% mutate(target_set = currName)
# }, reporter = "computing rocPValue"))
# 
# ws032$m006$rocPValue %>% saveRDS(paste0(ws032$workspaceDir, "misc/", "cache/", "032_006_roc_curves.rds"))
ws032$m006$rocPValue <- readRDS(paste0(ws032$workspaceDir, "misc/", "cache/", "032_006_roc_curves.rds"))


# ws032$m006$rocCiPValue <-  do.call(rbind, ws032$m006$targetSet %>% session$collectionUtils$lapplyWithName(function(currName, currTargets) {
#   session$evaluationUtils$aurocCI(ws032$m006$dea$rankingPvalue, currTargets, iteration = 1000) %>% mutate(target_set = currName)
# }, reporter = "computing rocCiPValue"))
# 
# ws032$m006$rocCiPValue %>% saveRDS(paste0(ws032$workspaceDir, "misc/", "cache/", "032_006_roc_ci.rds"))
ws032$m006$rocCiPValue <- readRDS(paste0(ws032$workspaceDir, "misc/", "cache/", "032_006_roc_ci.rds"))



# hypergeometric tests here

pData <- list()

pData$hits <- ws032$m006$dea %>% filter(adj.P.Val <= 0.1) %>% session$dataWrangler$extractColumn("entrez")

ws032$m006$targetSet %>% lapply(function(currSet) { 
  session$evaluationUtils$computeOra(hits = pData$hits, background = ws032$m006$dea$rankingPvalue, trueSet = currSet)
})



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# in-text calculations

pData <- list()

# all targets - AUROC & p-value

pData$ranksTibble <- tibble(gene = ws032$m006$dea$rankingPvalue, rank = 1:length(ws032$m006$dea$rankingPvalue)) %>% 
  mutate(curated = gene %in% ws032$m006$targetSet$all) %>% 
  mutate(curated = as.numeric(curated))
pData$result <- wilcox.test(rank ~ curated, data = pData$ranksTibble) 
pData$ns <- pData$ranksTibble %>% group_by(curated) %>% summarize(n = n()) %>% arrange(curated)
# auc
pData$result$statistic / (pData$ns$n[1] * pData$ns$n[2])
# pvalue
pData$result$p.value

# cns targets - AUROC & p-value

pData$ranksTibble <- tibble(gene = ws032$m006$dea$rankingPvalue, rank = 1:length(ws032$m006$dea$rankingPvalue)) %>% 
  mutate(curated = gene %in% ws032$m006$targetSet$tissueCns) %>% 
  mutate(curated = as.numeric(curated))
pData$result <- wilcox.test(rank ~ curated, data = pData$ranksTibble) 
pData$ns <- pData$ranksTibble %>% group_by(curated) %>% summarize(n = n()) %>% arrange(curated)
# auc
pData$result$statistic / (pData$ns$n[1] * pData$ns$n[2])
# pvalue
pData$result$p.value


# eye targets - AUROC & p-value

pData$ranksTibble <- tibble(gene = ws032$m006$dea$rankingPvalue, rank = 1:length(ws032$m006$dea$rankingPvalue)) %>% 
  mutate(curated = gene %in% ws032$m006$targetSet$tissueEye) %>% 
  mutate(curated = as.numeric(curated))
pData$result <- wilcox.test(rank ~ curated, data = pData$ranksTibble) 
pData$ns <- pData$ranksTibble %>% group_by(curated) %>% summarize(n = n()) %>% arrange(curated)
# auc
pData$result$statistic / (pData$ns$n[1] * pData$ns$n[2])
# pvalue
pData$result$p.value


# other targets - AUROC & p-value

pData$ranksTibble <- tibble(gene = ws032$m006$dea$rankingPvalue, rank = 1:length(ws032$m006$dea$rankingPvalue)) %>% 
  mutate(curated = gene %in% ws032$m006$targetSet$tissueOther) %>% 
  mutate(curated = as.numeric(curated))
pData$result <- wilcox.test(rank ~ curated, data = pData$ranksTibble) 
pData$ns <- pData$ranksTibble %>% group_by(curated) %>% summarize(n = n()) %>% arrange(curated)
# auc
pData$result$statistic / (pData$ns$n[1] * pData$ns$n[2])
# pvalue
pData$result$p.value


# multiple exp targets - AUROC & p-value

pData$ranksTibble <- tibble(gene = ws032$m006$dea$rankingPvalue, rank = 1:length(ws032$m006$dea$rankingPvalue)) %>% 
  mutate(curated = gene %in% ws032$m006$targetSet$expTypesMultiple) %>% 
  mutate(curated = as.numeric(curated))
pData$result <- wilcox.test(rank ~ curated, data = pData$ranksTibble) 
pData$ns <- pData$ranksTibble %>% group_by(curated) %>% summarize(n = n()) %>% arrange(curated)
# auc
pData$result$statistic / (pData$ns$n[1] * pData$ns$n[2])
# pvalue
pData$result$p.value


# single exp targets - AUROC & p-value

pData$ranksTibble <- tibble(gene = ws032$m006$dea$rankingPvalue, rank = 1:length(ws032$m006$dea$rankingPvalue)) %>% 
  mutate(curated = gene %in% ws032$m006$targetSet$expTypesSingle) %>% 
  mutate(curated = as.numeric(curated))
pData$result <- wilcox.test(rank ~ curated, data = pData$ranksTibble) 
pData$ns <- pData$ranksTibble %>% group_by(curated) %>% summarize(n = n()) %>% arrange(curated)
# auc
pData$result$statistic / (pData$ns$n[1] * pData$ns$n[2])
# pvalue
pData$result$p.value


# ORA at FDR = 0.1

ws032$m006$dea$hits$fdr %>% length()
ws032$m006$dea$rankingPvalue %>% length()

session$evaluationUtils$computeOra(ws032$m006$dea$hits$fdr, trueSet = ws032$m006$targetSet$all, background = ws032$m006$dea$rankingPvalue)
session$evaluationUtils$computeOra(ws032$m006$dea$hits$fdr, trueSet = ws032$m006$targetSet$tissueEye, background = ws032$m006$dea$rankingPvalue)

session$evaluationUtils$computeOra(ws032$m006$dea$hits$fdr, trueSet = ws032$m006$targetSet$tissueCns, background = ws032$m006$dea$rankingPvalue)
session$evaluationUtils$computeOra(ws032$m006$dea$hits$fdr, trueSet = ws032$m006$targetSet$tissueOther, background = ws032$m006$dea$rankingPvalue)

session$evaluationUtils$computeOra(ws032$m006$dea$hits$fdr, trueSet = ws032$m006$targetSet$expTypesMultiple, background = ws032$m006$dea$rankingPvalue)
session$evaluationUtils$computeOra(ws032$m006$dea$hits$fdr, trueSet = ws032$m006$targetSet$expTypesSingle, background = ws032$m006$dea$rankingPvalue)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - confidence intervals for AUROC P-value

pData <- list()

pData$outputFile <- paste0(ws032$m006$figurePath, "01.png")

pData$main <- ws032$m006$rocCiPValue

pData$main <- pData$main %>% 
  reshape2::melt(id = "target_set") %>% 
  as_tibble() 

pData$main %>% 
  session$graphingUtils$ggplot(aes(x = target_set, y = value)) +
  geom_line(aes(group = target_set)) + 
  geom_point(data = pData$main %>% filter(variable == "mean"), size = 3) +
  session$graphingUtils$tiltX(angle = 90) + 
  geom_hline(yintercept = 0.5, color = "grey70", linetype = "dashed") +
  ylim(0, 1) +
  scale_x_discrete(limits = c("all", "tissueCns", "tissueEye", "tissueOther",  "expTypesMultiple", "expTypesSingle"), 
                   labels = c("All", "CNS", "Eye", "Other", "Multiple", "Single")) +
  xlab("Target Category") + 
  ylab("AUROC") + 
  theme(legend.position = "none", axis.text = element_text(size = 15))

print(paste0("printing: ", pData$outputFile))

ggsave(pData$outputFile, 
       width = 6, height = 5)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - ROC curve with all PAX6 targets

pData <- list()

pData$outputFile <- paste0(ws032$m006$figurePath, "02.png")

pData$main <- ws032$m006$rocPValue

pData$main %>% 
  filter(target_set == "all") %>% 
  session$graphingUtils$ggplot(aes(x = false_positive_rate, y = recall)) + 
  geom_line(color = "grey20", size = 0.8) + 
  geom_abline(slope = 1, color = "grey70", linetype = "dashed") +
  xlab("False Positive Rate") + 
  ylab("Recall")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 10, height = 7.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# f05_S02 - ROC curves with PAX6 targets by tissue

pData <- list()

pData$outputFile <- paste0(ws032$m006$figurePath, "03.png")

pData$main <- ws032$m006$rocPValue

pData$main %>% 
  filter(target_set == "all") %>% 
  session$graphingUtils$ggplot(aes(x = false_positive_rate, y = recall)) + 
  geom_line(color = "grey20", size = 0.8) + 
  geom_line(data = pData$main %>% filter(grepl("tissue", target_set)), aes(color = target_set), size = 0.8) +
  scale_color_manual(values = c("red", "#6AA84F", "blue")) +
  geom_abline(slope = 1, color = "grey70", linetype = "dashed") +
  xlab("False Positive Rate") + 
  ylab("Recall") +
  theme(legend.position = "none")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 10, height = 7.5)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# f05_S03 - ROC curves with PAX6 targets by quality

pData <- list()

pData$outputFile <- paste0(ws032$m006$figurePath, "04.png")

pData$main <- ws032$m006$rocPValue

pData$main %>% 
  filter(target_set == "all") %>% 
  session$graphingUtils$ggplot(aes(x = false_positive_rate, y = recall)) + 
  geom_line(color = "grey20", size = 0.8) + 
  geom_line(data = pData$main %>% filter(grepl("expTypes", target_set)), aes(color = target_set), size = 0.8) +
  scale_color_manual(values = c("red", "blue")) +
  geom_abline(slope = 1, color = "grey70", linetype = "dashed") +
  xlab("False Positive Rate") + 
  ylab("Recall") +
  theme(legend.position = "none")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 10, height = 7.5)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# f05_S04 - volcano plot with CNS targets highlighted

pData <- list()

pData$outputFile <- paste0(ws032$m006$figurePath, "05.png")

pData$main <- do.call(rbind, ws032$m006$targetSet %>% session$collectionUtils$lapplyWithName(function(currName, currSet) {
  tibble(target_set = currName, 
         entrez = currSet)
}, reporter = pData$outputFile))

pData$dea <- ws032$m006$dea %>%
  dplyr::select(entrez, qvalue = adj.P.Val, logfc = logFC) %>% 
  mutate(qvalue_log = -log10(qvalue))

pData$main <- pData$dea %>%
  left_join(pData$main, by = "entrez")

pData$main %>% 
  session$graphingUtils$ggplot(aes(x = logfc, y = qvalue_log)) +
  geom_point(color = "grey80") +
  geom_point(data = pData$main %>% filter(target_set == "tissueOther"), color = "blue", alpha = 0.6, size = 3) +
  geom_point(data = pData$main %>% filter(target_set == "tissueCns"), color = "red", alpha = 0.6, size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  xlab("Fold Change (log2)") + 
  ylab("Q-Value (log10)")


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 10, height = 8)







