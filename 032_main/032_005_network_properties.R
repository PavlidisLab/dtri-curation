# run 032_SETUP.R first

ws032$m005 <- list()

ws032$m005$figurePath <- paste0(ws032$figuresDir, "f05_")

ws032$m005$neuroTfs <- ws032$data$s02.tsv %>% filter(Neurodev_TF == "TRUE") %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique()


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - bubble tally of DTRIs across mode of reuglation and tss distance

pData <- list()

pData$outputFile <- paste0(ws032$m005$figurePath, "01.png")

pData$dtris <- ws032$data$s04.tsv %>% 
  dplyr::select(tri_id = DTRI_ID, mode = Mode, tfbs_position = TFBS) %>% 
  unique() 

pData$bubbles <- pData$dtris %>% 
  group_by(mode, tfbs_position) %>% 
  summarize(n = n()) %>% 
  ungroup()

pData$bubbles$mode <- pData$bubbles$mode %>% factor(levels = c("activation", "repression", "both", "unknown"))
pData$bubbles$tfbs_position <- pData$bubbles$tfbs_position %>% factor(levels = c("proximal", "distal", "both", "unknown"))

pData$bubbles %>% 
  session$graphingUtils$ggplot(aes(x = mode, y = tfbs_position)) + 
  geom_point(aes(size = n)) + 
  geom_text(aes(label = n), vjust = -1.5, size = 7) +
  theme(legend.position = "none", axis.text = element_text(size = 18)) +
  scale_x_discrete(limits = c("activation", "repression", "both", "unknown"), 
                   labels = c("Activation", "Repression", "Both", "Unknown"))  +
  scale_y_discrete(limits = c("proximal", "distal", "both", "unknown"), 
                   labels = c("Proximal", "Distal", "Both", "Unknown")) +
  xlab(label = "Mode of Regulation") + 
  ylab(label = "TFBS Position") +
  session$graphingUtils$tiltX(angle = 90)


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - tally of mode of regulation for DTRIs per TF 

pData <- list()

pData$outputFile <- paste0(ws032$m005$figurePath, "02.png")

pData$dtris <- ws032$data$s04.tsv %>% 
  dplyr::select(tf_gene_h = TF_Entrez_ID_Human, tri_id = DTRI_ID, mode = Mode, tfbs_position = TFBS) %>% 
  unique() 

pData$bar <- pData$dtris %>% 
  group_by(tf_gene_h, mode) %>% 
  summarize(n_dtri = n()) %>% 
  ungroup()

pData$bar$mode <- pData$bar$mode %>% factor(levels = rev(c("activation", "repression", "both", "unknown")))

pData$sums <- pData$bar %>% 
  group_by(tf_gene_h) %>% 
  summarize(n_dtri = sum(n_dtri)) %>%
  arrange(desc(n_dtri)) %>% 
  left_join(ws032$genesLabs %>% dplyr::select(tf_gene_h = entrez, label = gene), by = "tf_gene_h")

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = tf_gene_h, y = n_dtri)) + 
  geom_bar(aes(fill = mode), stat = "identity") +
  scale_fill_manual(values = rev(c("#1F78B4", "#A6CEE3", "#B2DF8A", "grey90"))) +
  scale_x_discrete(limits = pData$sums$tf_gene_h[1:50], labels = pData$sums$label[1:50]) + 
  session$graphingUtils$tiltX(angle = 90) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 10, linetype = "dashed") 

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 22, height = 10)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# f06_S03 - tally of tfbs position for DTRIs per TF 

pData <- list()

pData$outputFile <- paste0(ws032$m005$figurePath, "03.png")

pData$dtris <- ws032$data$s04.tsv %>% 
  dplyr::select(tf_gene_h = TF_Entrez_ID_Human, tri_id = DTRI_ID, mode = Mode, tfbs_position = TFBS) %>% 
  unique() 

pData$bar <- pData$dtris %>% 
  group_by(tf_gene_h, tfbs_position) %>% 
  summarize(n_dtri = n()) %>% 
  ungroup()

pData$bar$tfbs_position <- pData$bar$tfbs_position %>% factor(levels = rev(c("proximal", "distal", "both", "unknown")))

pData$sums <- pData$bar %>% 
  group_by(tf_gene_h) %>% 
  summarize(n_dtri = sum(n_dtri)) %>%
  arrange(desc(n_dtri)) %>% 
  left_join(ws032$genesLabs %>% dplyr::select(tf_gene_h = entrez, label = gene), by = "tf_gene_h")

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = tf_gene_h, y = n_dtri)) + 
  geom_bar(aes(fill = tfbs_position), stat = "identity") +
  scale_fill_manual(values = rev(c("#1F78B4", "#A6CEE3", "#B2DF8A", "grey90"))) +
  scale_x_discrete(limits = pData$sums$tf_gene_h[1:50], labels = pData$sums$label[1:50]) + 
  session$graphingUtils$tiltX(angle = 90) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  geom_hline(yintercept = 10, linetype = "dashed") 

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 22, height = 10)

## in-text here ## ***

pData$bar %>% filter(tfbs_position == "distal") %>% 
  arrange(desc(n_dtri))



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# EXPERIMENT - network connectivity

# construct the curated network
# construct random networks, preserving degree, ~100?
# work out distances, etc

pData <- list()

pData$network$curated <- ws032$data$s04.tsv %>% 
  dplyr::select(source = TF_Entrez_ID_Human, target = Target_Entrez_ID_Human) %>% 
  unique() %>% 
  filter(source != target)

pData$graph$curated <- pData$network$curated %>% 
  mutate(edge = paste0(source, ",", target)) %>% 
  session$dataWrangler$extractColumn("edge") %>% 
  paste0(collapse = ",") %>% 
  strsplit(",") %>% 
  unlist() %>% 
  make_graph()

set.seed(1)
pData$network$random <- do.call(cbind, paste0("rand_", formatC(1:1000, digits = 3, flag = "0")) %>% session$collectionUtils$lapply(function(currIteration) {
  currTargets <- tibble(sample(pData$network$curated$target))
  names(currTargets) <- currIteration
  return(currTargets)
}, reporter = "random networks")) %>% as_tibble()

pData$network$all <- pData$network$curated %>% 
  dplyr::select(curated = target) %>% 
  cbind(pData$network$random) %>% 
  as_tibble()

pData$networkIds <- pData$network$all %>% names()
names(pData$networkIds) <- pData$networkIds

# construct the igraphs

pData$graphs <- pData$networkIds %>% session$collectionUtils$lapply(function(currId) {
  sourceNodes <- pData$network$curated$source
  currRandTargets <- pData$network$all[[currId]]
  currNetwork <- tibble(source = sourceNodes, target = currRandTargets)
  currNetwork %>% 
    mutate(edge = paste0(source, ",", target)) %>% 
    session$dataWrangler$extractColumn("edge") %>% 
    paste0(collapse = ",") %>% 
    strsplit(",") %>% 
    unlist() %>% 
    make_graph()
}, reporter = "graphs")

# compute all shortest distances



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - number of finite paths / "reachable" targets

fData <- list()

fData$outputFile <- paste0(ws032$m005$figurePath, "04.png")

# fData$main <- do.call(rbind, fData$graphs %>% session$collectionUtils$lapplyWithName(function(currName, currGraph) {
#   allDistances <- currGraph %>% distances(mode = "out") %>% as.vector()
#   val <- allDistances[allDistances != Inf] %>% length()
#   return(tibble(network = currName, value = val))
# }, reporter = "computing number of finite paths"))
# fData$main %>% saveRDS(paste0(ws032$workspaceDir, "misc/", "cache/", "032_005_n_finite_paths.rds"))

fData$main <- readRDS(paste0(ws032$workspaceDir, "misc/", "cache/", "032_005_n_finite_paths.rds"))

fData$curatedValue <- fData$main %>% filter(network == "curated") %>% session$dataWrangler$extractColumn("value")

fData$randValues <- fData$main %>% 
  filter(network != "curated") %>% 
  session$dataWrangler$extractColumn("value") 

fData$mean <- fData$randValues %>% mean()

fData$sd <- fData$main %>% 
  filter(network != "curated") %>% 
  session$dataWrangler$extractColumn("value") %>% 
  sd()

fData$pvalue <- max(1, (fData$randValues[which(fData$randValues > fData$curatedValue)] %>% length())) / length(fData$randValues)

fData$main %>%
  filter(network != "curated") %>% 
  session$graphingUtils$ggplot(aes(x = 0, y = value)) +
  geom_violin() + 
  geom_point(alpha = 0.1, position = "jitter") +
  geom_point(data = fData$main %>% filter(network == "curated"), size = 3, color = "red") + 
  theme(legend.position = "none") +
  ggtitle(paste0(paste0("N Finite Distances = ", fData$curatedValue, "\npvalue < ", fData$pvalue), 
                 "\n", paste0("Null: mean = ", round(fData$mean), ", sd = ", round(fData$sd)))) +
  ylab("Number of Finite Distances") + 
  xlab("") + 
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


print(paste0("printing: ", fData$outputFile))
ggsave(fData$outputFile, 
       width = 5, height = 12)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE -  Number of cliques of size 3 or larger

fData <- list()  

fData$outputFile <- paste0(ws032$m005$figurePath, "05.png")

fData$main <- do.call(rbind, pData$graphs %>% session$collectionUtils$lapplyWithName(function(currName, currGraph) {
  val <- currGraph %>% cliques(min = 3) %>% length()
  tibble(network = currName, value = val)
}, reporter = "computing number of cliques")) %>% 
  mutate(curated = grepl("curated", network))

fData$curatedValue <- fData$main %>% filter(curated) %>% session$dataWrangler$extractColumn("value")

fData$randValues <- fData$main %>% 
  filter(!curated) %>% 
  session$dataWrangler$extractColumn("value") 

fData$mean <- fData$randValues %>% mean()

fData$sd <- fData$main %>% 
  filter(!curated) %>% 
  session$dataWrangler$extractColumn("value") %>% 
  sd()

fData$pvalue <- max(1, (fData$randValues[which(fData$randValues > fData$curatedValue)] %>% length())) / length(fData$randValues)

fData$main %>%
  filter(!curated) %>% 
  session$graphingUtils$ggplot(aes(x = 0, y = value)) +
  geom_violin() + 
  geom_point(alpha = 0.1, position = "jitter") +
  geom_point(data = fData$main %>% filter(curated), size = 3, color = "red") + 
  theme(legend.position = "none") +
  ggtitle(paste0(paste0("N Cliques of Size 3+ = ", fData$curatedValue, "\npvalue < ", fData$pvalue), 
                 "\n", paste0("Null: mean = ", round(fData$mean), ", sd = ", round(fData$sd)))) +
  ylab("Number of Cliques of Size 3+") + 
  xlab("") + 
  theme(plot.title = element_text(size = 15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


print(paste0("printing: ", fData$outputFile))
ggsave(fData$outputFile, 
       width = 5, height = 12)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE -  Number of components

fData <- list() 

fData$outputFile <- paste0(ws032$m005$figurePath, "06.png")

fData$main <- do.call(rbind, pData$graphs %>% session$collectionUtils$lapplyWithName(function(currName, currGraph) {
  val <- currGraph %>% components() %>% igraph::groups() %>% length()
  tibble(network = currName, value = val)
}, reporter = "computing the number of components")) %>% 
  mutate(curated = grepl("curated", network))

fData$curatedValue <- fData$main %>% filter(curated) %>% session$dataWrangler$extractColumn("value")

fData$randValues <- fData$main %>% 
  filter(!curated) %>% 
  session$dataWrangler$extractColumn("value") 

fData$mean <- fData$randValues %>% mean()

fData$sd <- fData$main %>% 
  filter(!curated) %>% 
  session$dataWrangler$extractColumn("value") %>% 
  sd()

fData$pvalue <- max(1, (fData$randValues[which(fData$randValues > fData$curatedValue)] %>% length())) / length(fData$randValues)

fData$main %>%
  filter(!curated) %>% 
  session$graphingUtils$ggplot(aes(x = 0, y = value)) +
  geom_violin() + 
  geom_point(alpha = 0.1, position = "jitter") +
  geom_point(data = fData$main %>% filter(curated), size = 3, color = "red") + 
  theme(legend.position = "none") +
  ggtitle(paste0(paste0("Number of Components = ", fData$curatedValue, "\npvalue < ", fData$pvalue), 
                 "\n", paste0("Null: mean = ", round(fData$mean), ", sd = ", round(fData$sd)))) +
  ylab("Number of Components") + 
  xlab("") + 
  theme(plot.title = element_text(size = 15), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


print(paste0("printing: ", fData$outputFile))
ggsave(fData$outputFile, 
       width = 5, height = 12)









