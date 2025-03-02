#!/usr/bin/env Rscript

remove(list = ls())
package_vec = c("devtools", "SparseDOSSA2", "pkgmaker", "optparse", 
                "parallel", "stringi", "doParallel", "dplyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))


# # Set parameters manually
RandomEffect <- FALSE  

noParallel <- FALSE
nIterations <- 100
rSeed <- 1234
nCores <- 4
workingDirectory <- "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/maaslin3_benchmark"

# 嵌套多层 for 循环
for (metadataType in c("binary")) {
  
  for (nMicrobes in c(100)) {
    
    for (nSubjects in c(100, 20, 50, 200, 500, 1000)) {
      
      for (spikeMicrobes in c(0.1)) {
        
        for (effectSize in c(5, 2)) {
          
          for (effectPos in c(0.5)) {
            
            for (nMetadata in c(1)) {
              
              for (nPerSubject in c(1)) {
                
                for (readDepth in c(50000)) {
                  


# Load Utility Functions
pkgmaker::source_files(paste(workingDirectory, '/library/SD2_helper_functions.R', sep=''))

if (RandomEffect==TRUE){
  inputSubString<-'RandomEffect'
} else {
  inputSubString<-'noRandomEffect'
}

options("scipen"=10)
inputString<-paste(inputSubString, metadataType, nSubjects, nPerSubject, nMicrobes, spikeMicrobes, nMetadata, effectSize, effectPos, readDepth, sep='_')
options("scipen"=5)

# Create Input Directory
inputDirectory <- file.path(workingDirectory, 'Input', 'general_evaluations', 'SD2')
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

if (!outputs_already_exist){
  dir.create(this_params_folder, recursive = T)
  
  simlist<-trigger_sparseDOSSA2_Simulator(RandomEffect=RandomEffect,
                                          metadataType=metadataType,
                                          nSubjects=nSubjects,
                                          nPerSubject=nPerSubject,
                                          nMicrobes=nMicrobes,
                                          spikeMicrobes = spikeMicrobes, 
                                          nMetadata = nMetadata, 
                                          effectSize = effectSize,
                                          effectPos = effectPos,
                                          readDepth = readDepth,
                                          noParallel = noParallel,
                                          nIterations=nIterations,
                                          rSeed=rSeed,
                                          nCores=nCores)
  
  for (i in 1:length(simlist)) {
    current_sim <- simlist[[i]]
    write.table(current_sim$metadata, 
                file = paste0(this_params_folder, "/metadata_", i, ".tsv"), 
                sep = '\t')
    write.table(current_sim$features, 
                file = paste0(this_params_folder, "/abundance_", i, ".tsv"), 
                sep = '\t')
    write.table(current_sim$truth, 
                file = paste0(this_params_folder, "/truth_", i, ".tsv"), 
                sep = '\t',
                row.names = F)
  }
} else {
  print("Parameter combination already exists")
}

#######################################
# Delete Temporary sparseDOSSA Files  #
#######################################

if (file.exists("SyntheticMicrobiome-Counts.pcl")) file.remove("SyntheticMicrobiome-Counts.pcl")
if (file.exists("SyntheticMicrobiome.pcl")) file.remove("SyntheticMicrobiome.pcl")
if (file.exists("SyntheticMicrobiomeParameterFile.txt")) file.remove("SyntheticMicrobiomeParameterFile.txt")


                }
              }
            }
          }
        }
      }
    }
  }
}
