library(tidyverse)
library(ape)

setwd("Data")
otu_table <- read.csv("d1OTUtable.csv", row.names = 1)
taxonomy_table <- read.csv("d1Taxonomy.csv", row.names = 1)
tree_data <- read.tree("d1Tree.tree")
meta_table <- read.csv("d1Meta.csv", row.names = 1)
setwd("..")



## Data Preparation

taxonomy_table$taxa_to_genus <- apply(taxonomy_table, 1, function(x) paste(x[1:6], collapse = "."))

tree_data$root.edge <- 0
if ((sum(rownames(otu_table) %in% tree_data$tip.label) != length(tree_data$tip.label)) |
    (sum(tree_data$tip.label %in% rownames(otu_table)) != nrow(otu_table))) {
  
  otu_table <- otu_table[rownames(otu_table) %in% tree_data$tip.label, ]
  print(paste0(nrow(otu_table), " OTUs and ", length(tree_data$tip.label), " tip nodes."))
  
  if (nrow(otu_table) < length(tree_data$tip.label)) {
    padding <- matrix(0, ncol = ncol(otu_table), nrow = length(tree_data$tip.label) - nrow(otu_table))
    rownames(padding) <- tree_data$tip.label[!(tree_data$tip.label %in% rownames(otu_table))]
    colnames(padding) <- colnames(otu_table)
    otu_table <- rbind(otu_table, padding)
  }
}

if ((sum(rownames(taxonomy_table) %in% tree_data$tip.label) != length(tree_data$tip.label)) |
    (sum(tree_data$tip.label %in% rownames(taxonomy_table)) != nrow(taxonomy_table))) {
  
  taxonomy_table <- taxonomy_table[rownames(taxonomy_table) %in% tree_data$tip.label, ]
  print(paste0(nrow(taxonomy_table), " rows in taxonomy table, and ", length(tree_data$tip.label), " tip nodes."))
  
  if (nrow(taxonomy_table) < length(tree_data$tip.label)) {
    padding <- matrix(0, ncol = ncol(taxonomy_table), nrow = length(tree_data$tip.label) - sum(tree_data$tip.label %in% rownames(taxonomy_table)))
    rownames(padding) <- tree_data$tip.label[!(tree_data$tip.label %in% rownames(taxonomy_table))]
    colnames(padding) <- colnames(taxonomy_table)
    taxonomy_table <- rbind(taxonomy_table, padding)
  }
}

otu_table <- otu_table[,rownames(meta_table)]
taxonomy_table <- taxonomy_table[rownames(otu_table),]