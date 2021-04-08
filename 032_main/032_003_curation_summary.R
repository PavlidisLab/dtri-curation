# run 032_SETUP.R first

ws032$m003 <- list()
ws032$m003$figurePath <- paste0(ws032$figuresDir, "f03_")

ws032$m003$neuroTfs <- ws032$data$s02.tsv %>% filter(Neurodev_TF == "TRUE") %>% session$dataWrangler$extractColumn("TF_Entrez_ID") %>% unique()

ws032$m003$s03 <- ws032$data$s03.tsv %>% 
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
# FIGURE - N experiment, bubble chart, per experiment type & context type

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "01.png")

pData$records <- ws032$data$s03.tsv %>% 
  mutate(species_str = paste0(TF_Species, ",", Target_Species)) %>% 
  mutate(species = species_str %>% sapply(function(currSpecies) {
    mouse <- grepl("mouse", currSpecies)
    human <- grepl("human", currSpecies)
    if (mouse & human) { "mixed" }
    else if (mouse) { "mouse" }
    else if (human) { "human" }
    else { currSpecies }
  })) %>% dplyr::select(-TF_Species, -Target_Species, -species_str) %>% 
  unique() %>% 
  dplyr::select(tri_id = DTRI_ID, record_id = Experiment_ID, experiment_type = Experiment_Type, context_type = Context_Type, species) %>% 
  unique()

pData$experimentTypes <- pData$records$experiment_type %>% unique()
names(pData$experimentTypes) <- pData$experimentTypes

pData$contextTypes <- pData$records$context_type %>% unique()
names(pData$contextTypes) <- pData$contextTypes

pData$nOverlap <- do.call(rbind, pData$experimentTypes %>% session$collectionUtils$lapply(function(currExperimentType) {
  do.call(rbind, pData$contextTypes %>% session$collectionUtils$lapply(function(currContextType) {
    setA <- pData$records %>% filter(experiment_type == currExperimentType) %>% session$dataWrangler$extractColumn("record_id") %>% unique()
    setB <- pData$records %>% filter(context_type == currContextType) %>% session$dataWrangler$extractColumn("record_id") %>% unique()
    setIntersect <- setA %>% intersect(setB)
    speciesInfo <- pData$records %>% 
      dplyr::select(record_id, species) %>% 
      filter(record_id %in% setIntersect) %>% 
      unique() %>% 
      group_by(species) %>% 
      summarize(n_experiment = n()) %>% # ***I think this is the problem, should count the number of unique record_ids here instead
      full_join(tibble(species = c("human", "mouse", "mixed")), by = "species") %>% 
      session$dataWrangler$fillNa("n_experiment", 0)
    nIntersect <- setIntersect %>% length()
    return(tibble(experiment_type = currExperimentType, context_type = currContextType, species = speciesInfo$species, n_experiment = speciesInfo$n_experiment))
  }, reporter = "inner"))
}, reporter = "outer")) 

pData$nOverlap <- pData$nOverlap %>% 
  mutate(case = paste0(experiment_type, ".", context_type))

pData$nOverlap <- pData$nOverlap %>% 
  mutate(species = factor(species, levels = c("human", "mouse", "mixed")))

pData$nOverlap <- pData$nOverlap %>% 
  arrange(case, species)


# calculate the start and end angles for each pie
pData$nOverlap <- pData$nOverlap %>% 
  group_by(case) %>%
  mutate(n_experiment_total = sum(n_experiment)) %>% 
  mutate(nb_frac = 2*pi*cumsum(n_experiment)/n_experiment_total,
         start = lag(nb_frac, default = 0)) %>% 
  ungroup()


# convert into numeric for showing, attach labels later
pData$nOverlap <- pData$nOverlap %>% 
  mutate(coord_experiment_type = as.numeric(factor(experiment_type, levels = c("perturbation", "binding", "reporter"))), 
         coord_context_type = as.numeric(factor(context_type, levels = c("primary_tissue", "primary_cell", "cell_line", "in-vitro"))))


pData$expTypeLabels <- pData$nOverlap %>% dplyr::select(coord_experiment_type, experiment_type) %>% unique() %>% arrange(coord_experiment_type)

pData$contextTypeLabels <- pData$nOverlap %>% dplyr::select(coord_context_type, context_type) %>% unique() %>% arrange(coord_context_type)

scale = .5/sqrt(max(pData$nOverlap$n_experiment_total))

pData$labels <- pData$nOverlap %>% 
  group_by(case) %>% 
  summarize(coord_experiment_type = dplyr::first(coord_experiment_type), 
            coord_context_type = dplyr::first(coord_context_type), 
            n_experiment_total = dplyr::first(n_experiment_total)) %>% 
  ungroup()

pData$nOverlap %>% 
  session$graphingUtils$ggplot() +
  geom_arc_bar(aes(x0 = coord_context_type, y0 = coord_experiment_type, r0 = 0, r = sqrt(n_experiment_total)*scale, 
                   start = start, end = nb_frac, fill = species), 
               color = "grey40") +
  scale_x_discrete(name = "", limits = pData$contextTypeLabels$coord_context_type, labels = c("Primary Tissue", "Primary Cell", "Cell Line", "In-Vitro")) + 
  scale_y_discrete(name = "", limits = pData$expTypeLabels$coord_experiment_type, labels = c("TF Perturbation", "TF-DNA Binding", "TF-Reporter")) +
  geom_text(data = pData$labels,
            aes(label = n_experiment_total, 
                x = coord_context_type + sqrt(n_experiment_total + 10)*scale, 
                y = coord_experiment_type + sqrt(n_experiment_total + 10)*scale),
            size = 7) +
  theme(axis.text = element_text(size = 25)) +
  session$graphingUtils$tiltX(angle = 90)

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 11, height = 7.5)

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - number of targets per TF, for top TFs, fill = has_primary, cns

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "02.png")

pData$cnsTerms <- session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$UBERON) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$CL)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$UBERON)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$CL)) %>% 
  unique()

pData$tfTargetCount <- ws032$data$s03.tsv %>% 
  dplyr::select(tf_gene_h = TF_Entrez_ID_Human, tri_id = DTRI_ID, context_type = Context_Type, tissue_or_cell = Cell_Type) %>% 
  unique()

pData$tfTargetCount <- pData$tfTargetCount %>% 
  mutate(has_primary = grepl("primary", context_type), 
         has_cns = tissue_or_cell %in% pData$cnsTerms) %>% 
  group_by(tf_gene_h, tri_id) %>% 
  summarize(has_primary = any(has_primary), 
            has_cns = any(has_cns)) %>% 
  ungroup() %>% 
  unique()


pData$tfTargetCount <- pData$tfTargetCount %>% 
  mutate(context_str = paste0(has_primary, has_cns)) %>% 
  mutate(context = context_str %>% sapply(function(currContext) {
    if (currContext == "TRUETRUE") { "CNS" }
    else if (currContext == "TRUEFALSE") { "Primary tissue or cells" }
    else { "Cell line" }
  }))

pData$tfTargetCount <- pData$tfTargetCount %>% 
  group_by(tf_gene_h, context) %>% 
  summarize(n_target = n()) %>% 
  ungroup()

pData$tfTargetCount$context <- factor(pData$tfTargetCount$context, levels = rev(c("CNS", "Primary tissue or cells", "Cell line")))

pData$sums <- pData$tfTargetCount %>% 
  group_by(tf_gene_h) %>% 
  summarize(n_target = sum(n_target)) %>% 
  arrange(desc(n_target)) %>% 
  left_join(ws032$genesLabs %>% dplyr::select(tf_gene_h = entrez, gene))

pData$ndTfs <- tibble(tf_gene_h = ws032$m003$neuroTfs, n_target = -1)  

pData$tfTargetCount %>% 
  session$graphingUtils$ggplot(aes(x = tf_gene_h, y = n_target)) + 
  geom_bar(aes(fill = context), stat = "identity") +
  geom_point(data = pData$ndTfs, color = "#5382B2", size = 3) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_x_discrete(limits = pData$sums$tf_gene_h[1:50], labels = pData$sums$gene[1:50]) +
  scale_fill_manual(values = rev(c("#4682B7", "#8EB2D9", "grey85"))) +
  xlab("") +
  ylab("Number of Targets")


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 22, height = 6)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - number of TFs per target, for top 14 TFs, fill = has_primary, cns

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "03.png")

pData$cnsTerms <- session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$UBERON) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0001017", scope = session$ontologyUtils$ontologies$CL)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$UBERON)) %>% 
  c(session$ontologyUtils$getChildNodeRecursive("UBERON_0016879", scope = session$ontologyUtils$ontologies$CL)) %>% 
  unique()

pData$targetTfCount <- ws032$m003$s03 %>% 
  dplyr::select(target_gene_h = Target_Entrez_ID_Human, tri_id = DTRI_ID, context_type = Context_Type, tissue_or_cell = Cell_Type) %>% 
  filter(!session$dataWrangler$isNotApplicable(target_gene_h)) %>% 
  unique()

pData$targetTfCount <- pData$targetTfCount %>% 
  mutate(has_primary = grepl("primary", context_type), 
         has_cns = tissue_or_cell %in% pData$cnsTerms) %>% 
  group_by(target_gene_h, tri_id) %>% 
  summarize(has_primary = any(has_primary), 
            has_cns = any(has_cns)) %>% 
  ungroup() %>% 
  unique()

pData$targetTfCount <- pData$targetTfCount %>% 
  mutate(context_str = paste0(has_primary, has_cns)) %>% 
  mutate(context = context_str %>% sapply(function(currContext) {
    if (currContext == "TRUETRUE") { "CNS" }
    else if (currContext == "TRUEFALSE") { "Primary tissue or cells" }
    else { "Cell line" }
  }))

pData$targetTfCount <- pData$targetTfCount %>% 
  group_by(target_gene_h, context) %>% 
  summarize(n_tf = n()) %>% 
  ungroup()

pData$targetTfCount$context <- factor(pData$targetTfCount$context, levels = rev(c("CNS", "Primary tissue or cells", "Cell line")))

pData$sums <- pData$targetTfCount %>% 
  group_by(target_gene_h) %>% 
  summarize(n_tf = sum(n_tf)) %>% 
  arrange(desc(n_tf)) %>% 
  left_join(ws032$genesLabs %>% dplyr::select(target_gene_h = entrez, gene))

pData$ndTfs <- tibble(target_gene_h = ws032$m003$neuroTfs , n_tf = -1)  

pData$targetTfCount %>% 
  session$graphingUtils$ggplot(aes(x = target_gene_h, y = n_tf)) + 
  geom_bar(aes(fill = context), stat = "identity") +
  geom_point(data = pData$ndTfs, color = "#5382B2", size = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 23)) +
  geom_hline(yintercept = 10, linetype = "dashed") +
  scale_x_discrete(limits = pData$sums$target_gene_h[1:50], labels = pData$sums$gene[1:50]) +
  scale_fill_manual(values = rev(c("#4682B7", "#8EB2D9", "grey85"))) +
  xlab("") +
  ylab("Number of TF Regulators") + 
  theme(legend.position = "none")


print(paste0("printing: ", pData$outputFile))
ggsave(paste0(pData$outputFile), 
       width = 22, height = 6)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# CONSOLIDATING DTRI COMBINATION OF EXPERIMENTAL APPROACHES

pData <- list()

pData$perturbation <- ws032$data$s05.tsv %>% 
  left_join(ws032$m003$s03, by = "Experiment_ID") %>% 
  dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, species = Species, context_type = Context_Type, effect = Effect, induced = Induced) %>% 
  unique()

pData$perturbation <- pData$perturbation %>% 
  group_by(tri_id) %>% 
  summarize(effect_str = paste0(effect, collapse = ","), 
            induced_str = paste0(induced, collapse = ",")) %>% 
  mutate(perturbation = 1) %>% 
  dplyr::select(-effect_str, -induced_str)

pData$binding <- ws032$data$s06.tsv %>% 
  left_join(ws032$m003$s03, by = "Experiment_ID") %>% 
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
  left_join(ws032$m003$s03, by = "Experiment_ID")  %>% 
  dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, species = Species, context_type = Context_Type, 
                mutated = Mutated, binding_verified = Binding_Verified) %>% 
  unique()
  
pData$reporter <- pData$reporter %>% 
  group_by(tri_id) %>% 
  summarize(mutated_str = paste0(mutated, collapse = ",")) %>% 
  mutate(reporter = 1) %>% 
  dplyr::select(-mutated_str)

pData$overall <- ws032$m003$s03 %>% 
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


### in-text stuff ###

pData$master %>% 
  filter(perturbation == 1, binding == 1, reporter == 1, overall.cns == 1) %>% 
  session$dataWrangler$extractColumn("tri_id") %>%
  unique() %>% 
  length()

pData$master %>% 
  filter(perturbation == 1, binding == 1, reporter == 1, overall.primary == 1) %>% 
  session$dataWrangler$extractColumn("tri_id") %>%
  unique() %>% 
  length()



# ------------------------------------------------------------------------------
# FIGURE - Number of DTRIs done for different evidence combinations

pData$outputFile <- paste0(ws032$m003$figurePath, "04")

pData$plot <- pData$master %>% 
  dplyr::select(perturbation, binding, reporter, overall.primary, overall.cns) %>% 
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
pdf(file = paste0(pData$outputFile, ".pdf"), 
    width = 15, height = 10) # or other device
pData$plot
dev.off()


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# f03_S01 - distribution of TFs by n_target, out degree

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "05.png")

pData$a <- ws032$data$s04.tsv %>% 
  dplyr::select(TF_Entrez_ID_Human, DTRI_ID) %>% 
  unique() 

pData$a <- pData$a %>% 
  group_by(TF_Entrez_ID_Human) %>% 
  summarize(n_target = n()) %>% 
  arrange(desc(n_target))  

pData$a %>% 
  group_by(n_target) %>% 
  summarize(n_tf = n()) %>% 
  session$graphingUtils$ggplot(aes(x = n_target, y = n_tf)) +
  geom_bar(stat = "identity") +
  xlab("Number of Targets") +
  ylab("Number of TFs")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 10, height = 8)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# f03_S02 - distribution of targets by n_tf, in-degree

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "06.png")

pData$a <- ws032$data$s04.tsv %>% 
  dplyr::select(Target_Entrez_ID_Human, DTRI_ID) %>% 
  unique() 

pData$a <- pData$a %>% 
  group_by(Target_Entrez_ID_Human) %>% 
  summarize(n_tf = n()) %>% 
  arrange(desc(n_tf)) 

pData$a %>% 
  group_by(n_tf) %>% 
  summarize(n_target = n()) %>% 
  session$graphingUtils$ggplot(aes(x = n_tf, y = n_target)) +
  geom_text(aes(label = n_target), vjust = -1, size = 6) +
  geom_bar(stat = "identity") +
  ylim(0, 550) +
  xlab("Number of TFs") +
  ylab("Number of Targets")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 10, height = 8)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - PERTURBATION - n experiment by species

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "07.png")

pData$experiments <- ws032$data$s05.tsv %>% 
  dplyr::select(record_id = Experiment_ID, effect = Effect, induced = Induced) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$bar <- pData$experiments %>% 
  group_by(context_type, species) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(species) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$species <- factor(pData$bar$species, levels = rev(c("human", "mouse", "mixed")))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line")))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = species, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  scale_fill_manual(values = (c("#D9D9D9", "#9AB2D2", "#5D82AF"))) +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  scale_x_discrete(name = "", limits = c("human", "mouse", "mixed"), labels = c("Human", "Mouse", "Mixed")) + 
  ylim(0, 760) +
  session$graphingUtils$tiltX(angle = 90) +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Species")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - PERTURBATION - n experiment by context types

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "08.png")

pData$experiments <- ws032$data$s05.tsv %>% 
  dplyr::select(record_id = Experiment_ID, effect = Effect, induced = Induced) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$bar <- pData$experiments %>% 
  group_by(context_type, species) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(context_type) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$species <- factor(pData$bar$species, levels = rev(c("human", "mouse", "mixed")))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line")))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = context_type, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  scale_fill_manual(values = (c("#D9D9D9", "#9AB2D2", "#5D82AF"))) +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  scale_x_discrete(name = "", limits = c("primary_tissue", "primary_cell", "cell_line"), labels = c("Primary Tissue", "Primary Cell", "Cell Line")) + 
  ylim(0, 760) +
  session$graphingUtils$tiltX(angle = 90) +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Context Type")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - PERTURBATION - n experiment by perturbation mode

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "09.png")

pData$experiments <- ws032$data$s05.tsv %>% 
  dplyr::select(record_id = Experiment_ID, effect = Effect, induced = Induced) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$bar <- pData$experiments %>% 
  mutate(effect = effect %>% sapply(function(currEffect) {
    if (grepl("knockout", currEffect)) { "knockout" }
    else { currEffect }
  })) %>% 
  group_by(effect, context_type) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(effect) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line")))

pData$bar$effect <- factor(pData$bar$effect, levels = c("knockout", "knockdown", "overexpression"))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = effect, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  scale_fill_manual(values = (c("#D9D9D9", "#9AB2D2", "#5D82AF"))) +
  scale_x_discrete(name = "", limits = c("knockout", "knockdown", "overexpression"), labels = c("Knock Out", "Knock Down", "Overexpression")) + 
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  ylim(0, 530) +
  session$graphingUtils$tiltX(angle = 90) +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Mode")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - PERTURBATION - n experiment by perturbation effect

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "10.png")

pData$experiments <- ws032$data$s05.tsv %>% 
  dplyr::select(record_id = Experiment_ID, effect = Effect, induced = Induced) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$experiments <- pData$experiments %>% filter(!session$dataWrangler$isNotApplicable(induced))

pData$bar <- pData$experiments %>% 
  group_by(induced, context_type) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(induced) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line")))

pData$bar$induced <- factor(pData$bar$induced, levels = c("TRUE", "FALSE"))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = induced, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  scale_fill_manual(values = (c("#9AB2D2", "#5D82AF"))) +
  ylim(0, 420) +
  xlab("") +
  ylab("") +
  session$graphingUtils$tiltX(angle = 90) +
  scale_x_discrete(name = "", c("TRUE", "FALSE"), labels = c("Induced", "Germline")) + 
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Effect (Primary tissue or cells only)")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)





# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - BINDING - n experiment by species

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "11.png")

pData$experiments <- ws032$data$s06.tsv %>% 
  dplyr::select(record_id = Experiment_ID, method = Method, tf_source_type = TF_Source_Type) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$bar <- pData$experiments %>% 
  group_by(context_type, species) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(species) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$species <- factor(pData$bar$species, levels = rev(c("human", "mouse", "mixed")))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line", "in-vitro")))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = species, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  ylim(0, 850) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = (c("grey85", "grey70", "#9AB2D2", "#5D82AF"))) +
  session$graphingUtils$tiltX(angle = 90) +
  scale_x_discrete(name = "", limits = c("human", "mouse", "mixed"), labels = c("Human", "Mouse", "Mixed")) + 
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Species")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - BINDING - n experiment by context types

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "12.png")

pData$experiments <- ws032$data$s06.tsv %>% 
  dplyr::select(record_id = Experiment_ID, method = Method, tf_source_type = TF_Source_Type) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$bar <- pData$experiments %>% 
  group_by(context_type, species) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(context_type) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$species <- factor(pData$bar$species, levels = rev(c("human", "mouse", "mixed")))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line", "in-vitro")))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = context_type, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  ylim(0, 750) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = (c("grey85", "grey70", "#9AB2D2", "#5D82AF"))) +
  session$graphingUtils$tiltX(angle = 90) +
  scale_x_discrete(name = "", limits = c("primary_tissue", "primary_cell", "cell_line", "in-vitro"), labels = c("Primary Tissue", "Primary Cell", "Cell Line", "In-Vitro")) + 
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Context Type")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - BINDING - n experiment by method

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "13.png")

pData$experiments <- ws032$data$s06.tsv %>% 
  dplyr::select(record_id = Experiment_ID, method = Method, tf_source_type = TF_Source_Type) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$bar <- pData$experiments %>% 
  group_by(method, context_type) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(method) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line", "in-vitro")))

pData$bar$method <- factor(pData$bar$method, levels = c("chip", "emsa"))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = method, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  ylim(0, 1100) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = (c("grey85", "grey70", "#9AB2D2", "#5D82AF"))) +
  session$graphingUtils$tiltX(angle = 90) +
  scale_x_discrete(name = "", limits = c("chip", "emsa"), labels = c("ChIP-Assay", "EMSA")) + 
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Method")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - BINDING - n experiment tf source TF for EMSA

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "14.png")

pData$experiments <- ws032$data$s06.tsv %>% 
  dplyr::select(record_id = Experiment_ID, method = Method, tf_source_type = TF_Source_Type) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$experiments <- pData$experiments %>% filter(method == "emsa")

pData$bar <- pData$experiments %>% 
  group_by(tf_source_type, context_type) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(tf_source_type) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line", "in-vitro")))

pData$bar$tf_source_type <- factor(pData$bar$tf_source_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line", "unknown")))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = tf_source_type, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  ylim(0, 350) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = (c("grey85", "grey70", "#9AB2D2", "#5D82AF"))) +
  session$graphingUtils$tiltX(angle = 90) +
  scale_x_discrete(name = "", limits = c("primary_tissue", "primary_cell", "cell_line", "unknown"), labels = c("Primary Tissue", "Primary Cell", "Cell Line", "Unknown")) + 
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("TF Source Type (EMSA only)")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - REPORTER - n experiment by species

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "15.png")

pData$experiments <- ws032$data$s07.tsv %>% 
  dplyr::select(record_id = Experiment_ID, mutated = Mutated, verify = Binding_Verified) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$bar <- pData$experiments %>% 
  group_by(context_type, species) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(species) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$species <- factor(pData$bar$species, levels = rev(c("human", "mouse", "mixed")))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line")))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = species, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  scale_fill_manual(values = (c("#D9D9D9", "#9AB2D2", "#5D82AF"))) +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  scale_x_discrete(name = "", limits = c("human", "mouse", "mixed"), labels = c("Human", "Mouse", "Mixed")) + 
  ylim(0, 550) +
  session$graphingUtils$tiltX(angle = 90) +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Species")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - REPORTER - n experiment by context types

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "16.png")

pData$experiments <- ws032$data$s07.tsv %>% 
  dplyr::select(record_id = Experiment_ID, mutated = Mutated, verify = Binding_Verified) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$bar <- pData$experiments %>% 
  group_by(context_type, species) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(context_type) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$species <- factor(pData$bar$species, levels = rev(c("human", "mouse", "mixed")))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line")))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = context_type, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  scale_fill_manual(values = (c("#D9D9D9", "#9AB2D2", "#5D82AF"))) +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  scale_x_discrete(name = "", limits = c("primary_tissue", "primary_cell", "cell_line"), labels = c("Primary Tissue", "Primary Cell", "Cell Line")) + 
  ylim(0, 1000) +
  session$graphingUtils$tiltX(angle = 90) +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Context Type")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - REPORTER - n experiment by mutated or no

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "17.png")

pData$experiments <- ws032$data$s07.tsv %>% 
  dplyr::select(record_id = Experiment_ID, mutated = Mutated, verify = Binding_Verified) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$bar <- pData$experiments %>% 
  group_by(mutated, context_type) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(mutated) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line")))

pData$bar$mutated <- factor(pData$bar$mutated, levels = c("TRUE", "FALSE"))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = mutated, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  scale_fill_manual(values = (c("#D9D9D9", "#9AB2D2", "#5D82AF"))) +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  scale_x_discrete(name = "", limits = c("TRUE", "FALSE"), labels = c("True", "False")) + 
  ylim(0, 600) +
  session$graphingUtils$tiltX(angle = 90) +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Mutated")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - REPORTER - n experiment by verification 

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "18.png")

pData$experiments <- ws032$data$s07.tsv %>% 
  dplyr::select(record_id = Experiment_ID, mutated = Mutated, verify = Binding_Verified) %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, species = Species, context_type = Context_Type) %>% unique(), by = "record_id") %>% 
  unique()

pData$experiments <- pData$experiments %>% filter(!session$dataWrangler$isNotApplicable(verify))

pData$bar <- pData$experiments %>% 
  group_by(verify, context_type) %>% 
  summarize(n_experiment = n()) %>% 
  ungroup()

pData$sums <- pData$bar %>% 
  group_by(verify) %>% 
  summarize(n_experiment = sum(n_experiment)) %>% 
  arrange(desc(n_experiment))

pData$bar$context_type <- factor(pData$bar$context_type, levels = rev(c("primary_tissue", "primary_cell", "cell_line")))

pData$bar$verify <- factor(pData$bar$verify, levels = c("emsa", "putative"))

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = verify, y = n_experiment)) + 
  geom_bar(aes(fill = context_type), stat = "identity") +
  scale_fill_manual(values = (c("#D9D9D9", "#9AB2D2", "#5D82AF"))) +
  geom_text(data = pData$sums, aes(label = n_experiment), vjust = -1, size = 9) +
  scale_x_discrete(name = "", limits = c("emsa", "putative"), labels = c("EMSA", "Putative")) + 
  ylim(0, 320) +
  session$graphingUtils$tiltX(angle = 90) +
  xlab("") +
  ylab("") +
  theme(axis.text = element_text(size = 24), legend.position = "none") +
  ggtitle("Mutation Binding Verified\n(Mutated Only)")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6.5, height = 6.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# TISSUE DISTRIBUTION 


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Starting tissue ontology plots here - 


# SETUP FUNCTIONS / PROCESSING OF TISSUE DATA. ==== RUN ONCE AND THEN JUST GET IT FROM FILES ON DISK

ws032$m003$utils$getTissueExperiments <- function(experiments, scope, subset) {
  
  terms <- session$ontologyUtils$definitions[[scope]]$node %>% unique() %>% na.omit() %>% session$dataWrangler$attachNames()

  resultTibble <- do.call(rbind, 
                          terms %>% session$collectionUtils$lapplyWithName(function(currName, currTerm) {
                            allChildTerms <- session$ontologyUtils$getChildNodeRecursive(currTerm, scope)
                            records <- experiments %>% filter(Cell_Type %in% allChildTerms)
                            recordIds <- records %>% 
                              session$dataWrangler$extractColumn("Experiment_ID") %>% 
                              unique() 
                            tibble(node = currName, n_experiment = length(recordIds), records = recordIds)
                          }, subset))
  
  resultTibble <- resultTibble %>% 
    left_join(session$ontologyUtils$definitions[[scope]], by = "node") %>% 
    dplyr::select(node, label, n_experiment, records) %>% 
    filter(n_experiment > 0) %>% 
    arrange(desc(n_experiment))
  
  return(resultTibble)
}


# ws032$m003$tissueExperiments$primaryTissue <- ws032$data$s03.tsv %>%
#   filter(Context_Type == "primary_tissue") %>%
#   ws032$m003$utils$getTissueExperiments(scope = session$ontologyUtils$ontologies$UBERON, subset = "primary_tissue")

# ws032$m003$tissueExperiments$primaryTissue %>%
#   dplyr::select(record_id = records, everything()) %>%
#   write_tsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "primary_tissue.tsv"))
ws032$m003$tissueExperiments$primaryTissue <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "primary_tissue.tsv"))

# ws032$m003$tissueExperiments$primaryCell <-  ws032$data$s03.tsv %>%
#   filter(Context_Type == "primary_cell") %>%
#   ws032$m003$utils$getTissueExperiments(scope = session$ontologyUtils$ontologies$CL, subset = "primary_cell")

# ws032$m003$tissueExperiments$primaryCell %>%
#   dplyr::select(record_id = records, everything()) %>%
#   write_tsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "primary_cell.tsv"))
ws032$m003$tissueExperiments$primaryCell <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "primary_cell.tsv"))

ws032$m003$tissueExperiments$primaryCombined <- ws032$m003$tissueExperiments$primaryTissue %>%
  rbind(ws032$m003$tissueExperiments$primaryCell) %>%
  dplyr::select(node, label, record_id) %>%
  unique() %>%
  group_by(node) %>%
  mutate(n_experiment = n()) %>%
  ungroup() %>%
  dplyr::select(node, label, n_experiment, record_id)

# ws032$m003$tissueExperiments$primaryCombined %>%
#   dplyr::select(record_id, everything()) %>%
#   write_tsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "primary_combined.tsv"))
ws032$m003$tissueExperiments$primaryCombined <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "primary_combined.tsv"))

# map "future central nervous system" -UBERON_0016879 to "central nervous system" -UBERON_0001017
ws032$m003$tissueExperiments$primaryCombined <- ws032$m003$tissueExperiments$primaryCombined %>% 
  mutate(node = node %>% sapply(function(currNode) { if (currNode == "UBERON_0016879") { "UBERON_0001017" } else { currNode } }), 
         label = label %>% sapply(function(currLabel) { if (currLabel == "future central nervous system") { "central nervous system" } else { currLabel } }))




# ws032$m003$tissueExperiments$cellLine <- ws032$data$s03.tsv %>%
#   filter(Context_Type == "cell_line") %>%
#   ws032$m003$utils$getTissueExperiments(scope = session$ontologyUtils$ontologies$CLO, subset = "cell_line")

# ws032$m003$tissueExperiments$cellLine %>%
#   dplyr::select(record_id = records, everything()) %>%
#   write_tsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "cell_line.tsv"))
ws032$m003$tissueExperiments$cellLine <- session$dataWrangler$readTsv(paste0(ws032$workspaceDir, "misc/", "tissue_annotations/", "cell_line.tsv"))


ws032$m003$terms$anatomicalSystems <- tibble(node = c("UBERON_0001017", # CNS
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


ws032$m003$terms$anatomicalSystems <- ws032$m003$terms$anatomicalSystems %>% 
  left_join(session$ontologyUtils$getLabels(ws032$m003$terms$anatomicalSystems$node, scope = session$ontologyUtils$ontologies$CL), by = "node")
ws032$m003$terms$anatomicalSystems <- ws032$m003$terms$anatomicalSystems %>% 
  rbind(tibble(node = "other", label = "other", scope = "CL"))


ws032$m003$terms$cellLineCategories <- tibble(
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

ws032$m003$terms$cellLineCategories <- ws032$m003$terms$cellLineCategories %>% 
  left_join(session$ontologyUtils$getLabels(ws032$m003$terms$cellLineCategories$node, scope = session$ontologyUtils$ontologies$CLO), by = "node")
ws032$m003$terms$cellLineCategories <- ws032$m003$terms$cellLineCategories %>% 
  rbind(tibble(node = "other", label = "other", scope = "CLO"))

ws032$m003$terms$cellLineCategories <- ws032$m003$terms$cellLineCategories %>% 
  dplyr::select(node, scope, label.y = label) %>% 
  mutate(label.x = label.y %>% str_extract("(?<=immortal).*(?=cell line cell)") %>% trimws()) %>% 
  session$dataWrangler$mergeColumnsXy(colName = "label")

ws032$m003$terms$cellLineCategories <- ws032$m003$terms$cellLineCategories %>% 
  dplyr::select(node, scope, label.y = label) %>% 
  mutate(label.x = label.y %>% str_extract(".*(?=-derived)") %>% trimws()) %>% 
  session$dataWrangler$mergeColumnsXy(colName = "label")


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - primary combined, n_dtri across anatomical systems, fill = nd_tf, simple bar

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "19.png")

# first address coverage

pData$targetExperiments <- ws032$m003$tissueExperiments$primaryCombined %>% filter(node %in% ws032$m003$terms$anatomicalSystems$node) %>% 
  dplyr::select(record_id, node, label) %>% 
  unique()

pData$otherExperiments <- ws032$m003$tissueExperiments$primaryCombined %>% filter(!(record_id %in% pData$targetExperiments$record_id)) %>% unique()
pData$otherExperiments %>% group_by(node) %>% mutate(n_experiment = n()) %>% arrange(n_experiment) %>% ungroup()

pData$otherExperiments <- pData$otherExperiments %>% 
  dplyr::select(record_id) %>% 
  unique() %>% 
  mutate(node = "other", label = "other")

# attach "others"
pData$targetExperiments <- pData$targetExperiments %>% 
  rbind(pData$otherExperiments) %>% 
  unique()

pData$targetInteractions <- pData$targetExperiments %>% 
  left_join(ws032$data$s03.tsv %>% dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, tf_gene_h = TF_Entrez_ID_Human), by = "record_id") %>% 
  dplyr::select(tri_id, tf_gene_h, label, node) %>% 
  mutate(nd_tf = (tf_gene_h %in% ws032$m003$neuroTfs)) %>% 
  unique()


pData$a <- pData$targetInteractions %>% 
  group_by(node, nd_tf) %>%
  summarize(label = dplyr::first(label), n_dtri = n()) %>% 
  arrange(desc(n_dtri)) %>% 
  ungroup()

pData$b <- pData$a %>% 
  group_by(node) %>% 
  summarize(n_dtri = sum(n_dtri))

pData$a %>% 
  session$graphingUtils$ggplot(aes(x = node, y = n_dtri)) +
  geom_bar(stat = "identity", aes(fill = nd_tf)) +
  geom_text(data = pData$b, aes(label = n_dtri), hjust = -0.5, size = 6) +
  coord_flip() +
  scale_x_discrete(breaks = ws032$m003$terms$anatomicalSystems$node, 
                   labels = ws032$m003$terms$anatomicalSystems$label,
                   limits = rev(ws032$m003$terms$anatomicalSystems$node)) +
  scale_fill_manual(values = c("#CCCCCC", "#5382B2")) +
  ylim(0, 350) + 
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_blank(), legend.position = "none")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 5, height = 5.5)


# in-text calculation === 

# level of over representation of ND TFs in the CNS set
pData$background <- ws032$m003$s03$TF_Entrez_ID_Human %>% unique()
pData$positives <- pData$background %>% intersect(ws032$m003$neuroTfs)
pData$hits <- pData$targetInteractions %>% filter(node == "UBERON_0001017") %>% session$dataWrangler$extractColumn("tf_gene_h") %>% unique()
session$evaluationUtils$computeOra(pData$hits, trueSet = pData$positives, background = pData$background)


pData$interactions <- pData$targetInteractions %>% filter(tf_gene_h == "PAX6") %>% session$dataWrangler$extractColumn("tri_id") %>% unique()
pData$tissues <- ws032$m003$cpl$current$master %>% filter(tri_id %in% pData$interactions) %>% 
  dplyr::select(tf_gene_h, target_gene_h, tissue_or_cell, pmid) %>% 
  left_join(session$ontologyUtils$nodeLabels %>% dplyr::select(tissue_or_cell = node, label)) 


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - primary combined, n_dtri across anatomical systems, fill = nd_tf, simple bar, filter = embryonic

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "20.png")

# first address coverage

pData$targetExperiments <- ws032$m003$tissueExperiments$primaryCombined %>% filter(node %in% ws032$m003$terms$anatomicalSystems$node) %>% 
  dplyr::select(record_id, node, label) %>% 
  unique()

pData$otherExperiments <- ws032$m003$tissueExperiments$primaryCombined %>% filter(!(record_id %in% pData$targetExperiments$record_id)) %>% unique()
pData$otherExperiments %>% group_by(node) %>% mutate(n_experiment = n()) %>% arrange(n_experiment) %>% ungroup()

pData$otherExperiments <- pData$otherExperiments %>% 
  dplyr::select(record_id) %>% 
  unique() %>% 
  mutate(node = "other", label = "other")

# attach "others"
pData$targetExperiments <- pData$targetExperiments %>% 
  rbind(pData$otherExperiments) %>% 
  unique()

pData$targetInteractions <- pData$targetExperiments %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, tf_gene_h = TF_Entrez_ID_Human, age = Age), by = "record_id") %>% 
  dplyr::select(tri_id, tf_gene_h, label, node, age) %>% 
  mutate(nd_tf = (tf_gene_h %in% ws032$m003$neuroTfs)) %>% 
  unique()

pData$targetInteractions <- pData$targetInteractions %>% filter(grepl("^E.*", age))

pData$a <- pData$targetInteractions %>% 
  group_by(node, nd_tf) %>%
  summarize(label = dplyr::first(label), n_dtri = n()) %>% 
  arrange(desc(n_dtri)) %>% 
  ungroup()

pData$b <- pData$a %>% 
  group_by(node) %>% 
  summarize(n_dtri = sum(n_dtri))

pData$a %>% 
  session$graphingUtils$ggplot(aes(x = node, y = n_dtri)) +
  geom_bar(stat = "identity", aes(fill = nd_tf)) +
  geom_text(data = pData$b, aes(label = n_dtri), hjust = -0.5, size = 6) +
  coord_flip() +
  scale_x_discrete(breaks = ws032$m003$terms$anatomicalSystems$node, 
                   labels = ws032$m003$terms$anatomicalSystems$label,
                   limits = rev(ws032$m003$terms$anatomicalSystems$node)) +
  scale_fill_manual(values = c("#CCCCCC", "#5382B2")) +
  ylim(0, 350) + 
  xlab("") +
  ylab("") +
  theme(axis.text.x = element_blank(), legend.position = "none")


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 5, height = 5.5)


# in-text calculation === 

# level of over representation of ND TFs in the CNS set
pData$embryonicInteractions <- pData$targetInteractions %>% filter(node == "UBERON_0001017") %>% session$dataWrangler$extractColumn("tri_id") %>% unique()
pData$tissues <- ws032$m003$cpl$current$master %>% filter(tri_id %in% pData$embryonicInteractions) %>% 
  dplyr::select(tf_gene_h, target_gene_h, tissue_or_cell, pmid) %>% 
  left_join(session$ontologyUtils$nodeLabels %>% dplyr::select(tissue_or_cell = node, label)) 


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - primary combined, n experiment across anatomical systems, fill = nd_tf, simple bar, filter = embryonic

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "21.png")

# first address coverage

pData$targetExperiments <- ws032$m003$tissueExperiments$primaryCombined %>% filter(node %in% ws032$m003$terms$anatomicalSystems$node) %>% 
  dplyr::select(record_id, node, label) %>% 
  unique()

pData$otherExperiments <- ws032$m003$tissueExperiments$primaryCombined %>% filter(!(record_id %in% pData$targetExperiments$record_id)) %>% unique()
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
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, tf_gene_h = TF_Entrez_ID_Human, age = Age), by = "record_id") %>% 
  dplyr::select(record_id, tf_gene_h, label, node, age) %>% 
  mutate(nd_tf = (tf_gene_h %in% ws032$m003$neuroTfs)) %>% 
  unique()

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
  scale_x_discrete(breaks = ws032$m003$terms$anatomicalSystems$node, 
                   labels = ws032$m003$terms$anatomicalSystems$label,
                   limits = rev(ws032$m003$terms$anatomicalSystems$node)) +
  scale_fill_manual(values = c("#CCCCCC", "#5382B2")) +
  ylim(0, 350) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), legend.position = "none")


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 5, height = 5.5)


# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# FIGURE - primary combined, n experiment across anatomical systems, fill = nd_tf, simple bar, filter = embryonic

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "22.png")

# first address coverage

pData$targetExperiments <- ws032$m003$tissueExperiments$primaryCombined %>% filter(node %in% ws032$m003$terms$anatomicalSystems$node) %>% 
  dplyr::select(record_id, node, label) %>% 
  unique()

pData$otherExperiments <- ws032$m003$tissueExperiments$primaryCombined %>% filter(!(record_id %in% pData$targetExperiments$record_id)) %>% unique()
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
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, tf_gene_h = TF_Entrez_ID_Human, age = Age), by = "record_id") %>% 
  dplyr::select(record_id, tf_gene_h, label, node, age) %>% 
  mutate(nd_tf = (tf_gene_h %in% ws032$m003$neuroTfs)) %>% 
  unique()

pData$targetExperiments <- pData$targetExperiments %>% filter(grepl("^E.*", age))

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
  scale_x_discrete(breaks = ws032$m003$terms$anatomicalSystems$node, 
                   labels = ws032$m003$terms$anatomicalSystems$label,
                   limits = rev(ws032$m003$terms$anatomicalSystems$node)) +
  scale_fill_manual(values = c("#CCCCCC", "#5382B2")) +
  ylim(0, 350) + 
  xlab("") +
  ylab("") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), legend.position = "none")


print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 5, height = 5.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# f03_S08-A - n_dtri across different cell lines, fill = nd_tf, simple bar

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "23.png")

# first address coverage

pData$targetExperiments <- ws032$m003$tissueExperiments$cellLine %>% filter(node %in% ws032$m003$terms$cellLineCategories$node) %>% unique() %>% 
  dplyr::select(record_id, node, label)

pData$otherExperiments <- ws032$m003$tissueExperiments$cellLine %>% filter(!(record_id %in% pData$targetExperiments$record_id)) %>% unique()
# pData$otherExperiments %>% group_by(node) %>% summarize(label = first(label), n_experiment = n()) %>% ungroup() %>% arrange(n_experiment) 

pData$otherExperiments <- pData$otherExperiments %>% 
  dplyr::select(record_id) %>% 
  unique() %>% 
  mutate(node = "other", label = "other")

# attach "others"
pData$targetExperiments <- pData$targetExperiments %>% 
  rbind(pData$otherExperiments)

pData$targetInteractions <- pData$targetExperiments %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, tf_gene_h = TF_Entrez_ID_Human), by = "record_id") %>% 
  dplyr::select(tri_id, tf_gene_h, label, node) %>% 
  mutate(nd_tf = (tf_gene_h %in% ws032$m003$neuroTfs)) %>% 
  unique()


pData$a <- pData$targetInteractions %>% 
  group_by(node, nd_tf) %>%
  summarize(label = dplyr::first(label), n_dtri = n()) %>% 
  arrange(desc(n_dtri)) %>% 
  ungroup()

pData$b <- pData$a %>% 
  group_by(node) %>% 
  summarize(n_dtri = sum(n_dtri))

pData$a %>% 
  session$graphingUtils$ggplot(aes(x = node, y = n_dtri)) +
  geom_bar(stat = "identity", aes(fill = nd_tf)) +
  geom_text(data = pData$b, aes(label = n_dtri), hjust = -0.5, size = 6) +
  coord_flip() +
  scale_x_discrete(breaks = ws032$m003$terms$cellLineCategories$node, 
                   labels = ws032$m003$terms$cellLineCategories$label,
                   limits = rev(ws032$m003$terms$cellLineCategories$node)) +
  scale_fill_manual(values = c("#CCCCCC", "#5382B2")) +
  ylim(0, 550) + 
  xlab("") +
  ylab("")  +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), legend.position = "none")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6, height = 5.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# FIGURE - n_experiment across different cell lines, fill = nd_tf, simple bar

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "24.png")

# first address coverage

pData$targetExperiments <- ws032$m003$tissueExperiments$cellLine %>% filter(node %in% ws032$m003$terms$cellLineCategories$node) %>% unique() %>% 
  dplyr::select(record_id, node, label)

pData$otherExperiments <- ws032$m003$tissueExperiments$cellLine %>% filter(!(record_id %in% pData$targetExperiments$record_id)) %>% unique()
# pData$otherExperiments %>% group_by(node) %>% summarize(label = first(label), n_experiment = n()) %>% ungroup() %>% arrange(n_experiment) 

pData$otherExperiments <- pData$otherExperiments %>% 
  dplyr::select(record_id) %>% 
  unique() %>% 
  mutate(node = "other", label = "other")

# attach "others"
pData$targetExperiments <- pData$targetExperiments %>% 
  rbind(pData$otherExperiments)

pData$targetExperiments <- pData$targetExperiments %>% 
  left_join(ws032$m003$s03 %>% dplyr::select(record_id = Experiment_ID, tri_id = DTRI_ID, tf_gene_h = TF_Entrez_ID_Human), by = "record_id") %>% 
  dplyr::select(record_id, tf_gene_h, label, node) %>% 
  mutate(nd_tf = (tf_gene_h %in% ws032$m003$neuroTfs)) %>% 
  unique()

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
  scale_x_discrete(breaks = ws032$m003$terms$cellLineCategories$node, 
                   labels = ws032$m003$terms$cellLineCategories$label,
                   limits = rev(ws032$m003$terms$cellLineCategories$node)) +
  scale_fill_manual(values = c("#CCCCCC", "#5382B2")) +
  ylim(0, 550) + 
  xlab("") +
  ylab("")  +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), legend.position = "none")

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 6, height = 5.5)



# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# FIGURE - distribution of n_dtri per pmid

pData <- list()

pData$outputFile <- paste0(ws032$m003$figurePath, "25.png")

pData$dtriCount <- ws032$m003$s03 %>% 
  dplyr::select(pmid = PubMed_ID, tri_id = DTRI_ID) %>% 
  unique()

pData$dtriCount <- pData$dtriCount %>% 
  group_by(pmid) %>% 
  summarize(n_dtri = n()) %>% 
  arrange(desc(n_dtri))

pData$bar <- pData$dtriCount %>% 
  group_by(n_dtri) %>% 
  summarize(n_pmid = n())

pData$mean <- pData$dtriCount$n_dtri %>% mean()

pData$bar %>% 
  session$graphingUtils$ggplot(aes(x = n_dtri, y = n_pmid)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(label = n_pmid), vjust = -1, size = 6) +
  ylim(0, 510) +
  xlab("Number of DTRIs") +
  ylab("Number of Papers") + 
  geom_vline(xintercept = pData$mean, linetype = "dashed") + 
  geom_text(x = 4, y = 400, label = paste0("Mean = ", round(pData$mean, 2)), size = 6)

print(paste0("printing: ", pData$outputFile))
ggsave(pData$outputFile, 
       width = 10, height = 8)







