library(tidyverse)
library(curatedMetagenomicData)
library(ape)


se_abs <- curatedMetagenomicData("HMP_2019_ibdmdb.relative_abundance", dryrun = FALSE, counts = TRUE)[[1]]


meta_table <- colData(se_abs) %>%
  as.data.frame()%>%
  filter(visit_number == 1) %>%
  .[, c("age", "disease", "antibiotics_current_use")]

meta_table$disease <- as.factor(meta_table$disease)
meta_table$disease <- relevel(meta_table$disease, "healthy")

otu_table <- as.data.frame(assay(se_abs))[, rownames(meta_table)]
taxonomy_table <- data.frame(taxa_to_genus = character(nrow(otu_table)))
taxonomy_table$taxa_to_genus <- sub(".s__.*", "", rownames(otu_table))

tree_data <- rowTree(se_abs)


pruned_tree <- drop.tip(tree_data, setdiff(tree_data$tip.label, rownames(otu_table)))





# unique(sampleMetadata$study_name)