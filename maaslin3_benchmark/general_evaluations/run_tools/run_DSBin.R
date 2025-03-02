#!/usr/bin/env Rscript

remove(list = ls())
package_vec = c("devtools", "pkgmaker", "optparse", "ALDEx2",
                "parallel", "stringi", "doParallel", "plyr", "tidyr", "ape", "vegan", "GUniFrac")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))




RandomEffect <- FALSE  # TRUE if using a longitudinal design
metadataType <- "binary"
# nSubjects <- 200
nPerSubject <- 1
nMicrobes <- 100
spikeMicrobes <- 0.1
nMetadata <- 1
# effectSize <- 5
effectPos <- 0.5
qvalMethod <- 'BY'
readDepth <- 50000
noParallel <- FALSE
nIterations <- 100
rSeed <- 1234
nCores <- 10
workingDirectory <- "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/maaslin3_benchmark"
generator <- "SD2"
taxonomy_table = NULL


setwd(workingDirectory)


source("library/Split_fun.R")
source("library/DS_Mann_fun.R")
source("library/helper_function.R")


for (nSubjects in c(100, 200, 500, 1000)) {
  
  
  for (effectSize in c(5)) {
    
    for (nReps in c(50, 100, 200)) {
      
      if (RandomEffect==TRUE){
        inputSubString<-'RandomEffect'
      } else {
        inputSubString<-'noRandomEffect'
      }
      
      options("scipen"=10)
      inputString<-paste(inputSubString, metadataType, nSubjects, nPerSubject, nMicrobes, spikeMicrobes, nMetadata, effectSize, effectPos, readDepth, sep='_')
      options("scipen"=5)
      
      # Create Input Directory
      inputDirectory <- file.path(workingDirectory, 'Input', 'general_evaluations', generator)
      this_params_folder <- file.path(inputDirectory, inputString)
      
      outputs_already_exist <- TRUE
      if (!dir.exists(this_params_folder)) {
        outputs_already_exist <- FALSE
      } else {
        for (i in 1:length(nIterations)) {
          if (!file.exists(paste0(this_params_folder, "/metadata_", i, ".tsv")) |
              !file.exists(paste0(this_params_folder, "/abundance_", i, ".tsv")) |
              !file.exists(paste0(this_params_folder, "/truth_", i, ".tsv"))) {
            outputs_already_exist <- FALSE
          }
        }
      }
      
      # Load Dataset
      if (!outputs_already_exist) stop('The input file does not exist. Generate the dataset first.') 
      
      # How Many Datasets to Run In A List 
      reps<-1:nIterations
      
      # Create Output Directory
      outputDirectory <- file.path(workingDirectory, 'Output', 'general_evaluations', generator)
      outputString <- paste(inputString, qvalMethod, nReps, "DSBin", sep="_")
      this_output_folder <- file.path(outputDirectory, outputString)
      
      outputs_already_exist <- TRUE
      if (!dir.exists(this_output_folder)) {
        outputs_already_exist <- FALSE
        dir.create(this_output_folder, recursive = T, showWarnings = F)
      } else {
        for (i in 1:length(nIterations)) {
          if (!file.exists(paste0(this_output_folder, "/associations_", i, ".tsv"))) {
            outputs_already_exist <- FALSE
          }
        }
      }
      
      # Run the Model Only if the Output Does Not Exist
      if (!outputs_already_exist){
        no_cores <- nCores
        cl <- makeCluster(no_cores)
        registerDoParallel(cl)
        
        f <- foreach(i = reps, .packages = c("MASS", "stringi", 'dplyr', 'parallel', 'doParallel', "ape", "vegan", "GUniFrac", "tidyverse", "glmnet"),
                     .export = c("computeR2","permuteRows", "compDist", "analys", "find_tau","singleSplit")) %dopar% {
                       metadata <- read.csv(paste0(this_params_folder, "/metadata_", i, ".tsv"), sep = "\t")
                       abundance <- read.csv(paste0(this_params_folder, "/abundance_", i, ".tsv"), sep = "\t")
                       truth <- read.csv(paste0(this_params_folder, "/truth_", i, ".tsv"), sep = "\t")
                       
                       # Remove features or samples never present
                       abundance <- abundance[apply(abundance, 1, var) != 0,]
                       abundance <- abundance[,apply(abundance, 2, var) != 0]
                       
                       metadata <- dplyr::select(metadata, -ID)
                       
                       
                       otu_table = abundance
                       meta_table = metadata
                       
                       colnames(meta_table)[1] <- "BinOutcomes"
                       
                       taxonomy_table = NULL
                       
                       selected_taxon <- tryCatch({DSBin(X = otu_table, y = meta_table$BinOutcomes, taxonomy_table, num_split = nReps, qval_bound = 0.05)
                       }, error = function(e) {
                         # If locom threw an error:
                         message("Error in DSBin: ", e$message)
                         return(character(0))   # Return NULL so we can detect it below
                       })
                       
                       pval <- numeric(nrow(abundance))
                       names(pval) <- rownames(abundance)
                       pval <- ifelse(names(pval) %in% selected_taxon, 0.04, 0.8)
                       qval <- pval
                       
                       outputs <- data.frame(taxon = rownames(abundance),
                                             metadata = colnames(metadata)[1],
                                             effect_size = 1,
                                             pval = pval,
                                             qval = qval,
                                             associations = "abundance")
                       return(outputs)
                     }
        
        for (i in 1:nIterations) {
          write.table(f[[i]],
                      file = paste0(this_output_folder, "/associations_", i, ".tsv"),
                      sep = '\t',
                      row.names = F)
        }
        
        # Stop the Cluster 
        stopCluster(cl)
      }
    }
  }
}