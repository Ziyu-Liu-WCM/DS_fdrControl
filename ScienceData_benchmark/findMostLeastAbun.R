library(tidyverse)
library(ape)
library(GUniFrac)
library(MiRKAT)
library(vegan)
library(glmnet)
library(parallel)
library(doParallel)


## Read files
source("helper_function.R")
source("preprocess.R")
str(taxonomy_table)

otu_table$taxa_to_genus <- taxonomy_table$taxa_to_genus

genus_count <- otu_table %>%
  group_by(taxa_to_genus) %>%
  summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
  column_to_rownames("taxa_to_genus")

otu_table <- dplyr::select(otu_table, -taxa_to_genus)

genus_abun<-otutabletoabundance(genus_count)
apply(genus_abun, 2, sum)
mostAbun<-names(sort(apply(genus_abun, 1, mean), decreasing = TRUE)[1:20])
leastAbun<-names(sort(apply(genus_abun, 1, mean))[1:20])

save(mostAbun, leastAbun, file = "mostLeastAbun.RData")
