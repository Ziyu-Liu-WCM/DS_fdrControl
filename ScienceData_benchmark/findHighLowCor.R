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

# find the top 10 features that on avearage has the highest correlation with other features and the lowest correlation with other features

cor_matrix <- cor(t(genus_abun), use = "complete.obs")
averageCor<-apply(cor_matrix, 1, function(x) mean(abs(x)))
highCor<-names(sort(averageCor, decreasing = TRUE)[1:20])
lowCor<-names(sort(averageCor, decreasing = FALSE)[1:20])


save(highCor, lowCor, file = "highLowCor.RData")
