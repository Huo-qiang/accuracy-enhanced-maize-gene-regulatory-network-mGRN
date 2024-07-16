# Load GAMER.
zmaGO <- readr::read_tsv("https://github.com/timedreamer/public_dataset/raw/master/maize.B73.AGPv4.aggregate.gaf.gz", skip=1) %>% select(term_accession, db_object_id)

# Load GO_ID to GO_annotation mapping.
GOname <- readr::read_tsv("https://github.com/timedreamer/public_dataset/raw/master/agriGOv2_GOConsortium_term_v201608.txt.gz", col_names = c("GO","type","name","number")) %>% select(GO, name)
#
dap_genesets <- readRDS("G:/DAP/dap_genesets.rds")
# Run GO enrichment.

results_list_root <- list()
results_list_stem <- list()
results_list_leaf<- list()
results_list_ear <- list()
results_list_tassel <- list()
results_list_endosperm <- list()
# Perform enricher operations on each gene set
for (i in seq_along(geneSets_root)) {
  ego <- enricher(gene = geneSets_root[[i]],
                  TERM2GENE = zmaGO,
                  TERM2NAME = GOname,
                  pAdjustMethod = "BH")
  
  # Get the name of the current gene set
  gene_set_name <- names(geneSets_root)[i]
  
  # Store results in a list, including gene set names
  results_list_root[[gene_set_name]] <- ego@result
}

# Merge the data in the results list into a single dataframe
combined_results_root <- do.call(rbind, results_list_root)

# save
write.csv(combined_results, file = "enrichment_results_a_grn.csv")
##
