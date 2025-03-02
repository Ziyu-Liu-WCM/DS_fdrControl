
#######################################################
#Please Disregard the qvalMethod <- 'BY' Parameter!
#######################################################


#!/usr/bin/env Rscript

remove(list = ls())
package_vec = c("devtools", "SparseDOSSA2", "pkgmaker", "optparse", 
                "parallel", "stringi", "doParallel", "dplyr")
invisible(suppressPackageStartupMessages(lapply(package_vec, require, character.only = TRUE)))

# Set parameters manually
RandomEffect <- FALSE  # TRUE if using a longitudinal design
metadataType <- "binary"
# nSubjects <- 100
# nPerSubject <- 1
nMicrobes <- 100
# spikeMicrobes <- 0.1
# nMetadata <- 5
# effectSize <- 5
# effectPos <- 0.5
qvalMethod <- 'BY'
# readDepth <- 50000
noParallel <- FALSE
nIterations <- 100
rSeed <- 1234
nCores <- 4
workingDirectory <- "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/maaslin3_benchmark"
generator <- "SD2"
setwd(workingDirectory)

overall_result_df <- data.frame(matrix(nrow = 0, ncol = 0))

# Load Utility Functions
pkgmaker::source_files(paste(workingDirectory, '/library/run_evaluation_helpers.R', sep=''))


for (nSubjects in c(100, 200)) {

  for (spikeMicrobes in c(0.1)) {
    
    for (effectSize in c(5)) {
      
      for (effectPos in c(0.5)) {
        
        for (nMetadata in c(1)) {
          
          for (nPerSubject in c(1)) {
            
            for (readDepth in c(50000)) {

if (RandomEffect==TRUE){
  inputSubString<-'RandomEffect'
} else {
  inputSubString<-'noRandomEffect'
}
inputString<-paste(inputSubString, metadataType, nSubjects, nPerSubject, nMicrobes, spikeMicrobes, nMetadata, effectSize, effectPos, readDepth, sep='_')


# parameters <- "general_evaluations/data_generation/SD2_me.txt"
# filename_pieces <- unlist(strsplit(parameters, '/'))

# nIterations = 100
# 
# parameter_text <- readLines(parameters)

param_list = list(nMicrobes, nSubjects, spikeMicrobes, effectSize, effectPos, nMetadata, nPerSubject, readDepth)
names(param_list) <- c("nMicrobes", "nSubjects", "spikeMicrobes", "effectSize", "effectPos", "nMetadata", "nPerSubject", "readDepth")
# param_list = list()
# for (parameter_line in parameter_text) {
#   parts = trimws(unlist(strsplit(parameter_line, ':')))
#   vals = unlist(strsplit(parts[2], ' '))
#   param_list[[parts[1]]] = vals
# }

metadata_types = metadataType
# metadata_types = param_list[['metadataType']]
# param_list[['metadataType']] <- NULL


inputDirectory <- file.path('Input', 'general_evaluations', generator)
this_params_folder <- file.path(inputDirectory, inputString)
outputDirectory <- file.path('Output', 'general_evaluations', generator)
tmp_files <- list.files(outputDirectory)
tmp_files <- tmp_files[grepl(paste0(inputString, "_"), tmp_files)]
tools <- gsub(paste0(inputString, "_"), '', tmp_files, perl = T)

this_output_folder <- file.path(outputDirectory, inputString)
tools <- c("BY_50_DSBin", "BY_100_DSBin", "BY_200_DSBin")
general_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
maaslin3_results_df <- data.frame(matrix(nrow = 0, ncol = 0))
metrics <- c("Precision", "Recall", "Precision\n(common taxa)", "Recall\n(common taxa)", 
              "AUC", "AUC (Common taxa)")

for (i in 1:nIterations) {
  possible_error <- tryCatch({
    metadata <- read.csv(paste0(this_params_folder, "/metadata_", i, ".tsv"), sep = "\t")
    metadata <- metadata[,colnames(metadata) != 'ID']
    abundance <- read.csv(paste0(this_params_folder, "/abundance_", i, ".tsv"), sep = "\t")
    truth <- read.csv(paste0(this_params_folder, "/truth_", i, ".tsv"), sep = "\t")
  }, error = function(err) {
    err
  })
  if(inherits(possible_error, "error")) next
  
  truth <- prepare_truth(truth, generator)
  for (tool in tools) {
    possible_error <- tryCatch({
      associations <- read.csv(paste0(this_output_folder, '_', tool, "/associations_", i, ".tsv"), sep = "\t")
    }, error = function(err) {
      err
    })
    if(inherits(possible_error, "error")) next
  
    
    new_row <- c(unweighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                 weighted_precision_recall(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                 pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata),
                 weighted_pval_auc(truth, prepare_associations_general(associations, tool, generator), abundance, metadata)
                 )
  new_row <- lapply(new_row, function(x) {x})
  names(new_row) <- metrics
  
  new_row[['tool']] <- tool
  new_row[['iter']] <- i
  new_row <- c(new_row, param_list)
  general_results_df <- plyr::rbind.fill(general_results_df, data.frame(new_row, check.names = F))
  
  if (grepl('Maaslin3', tool) & generator == 'SD2') {
    tmp_truth <- truth[truth$associations == 'abundance',]
    tmp_associations <- prepare_associations_maaslin3(associations, tool)
    tmp_associations <- tmp_associations[tmp_associations$association == 'abundance',]
    new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                 weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                 effect_size_error_maaslin3(tmp_truth, tmp_associations),
                 effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                 pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                 weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                 issue_prop(associations[associations$associations == 'abundance',], tool))
    new_row <- lapply(new_row, function(x) {x})
    names(new_row) <- metrics
    
    new_row[['tool']] <- tool
    new_row[['iter']] <- i
    new_row[['association_type']] <- 'abundance'
    new_row <- c(new_row, param_list)
    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
    
    tmp_truth <- truth[truth$associations == 'prevalence',]
    tmp_associations <- prepare_associations_maaslin3(associations, tool)
    tmp_associations <- tmp_associations[tmp_associations$association == 'prevalence',]
    new_row <- c(unweighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                 weighted_precision_recall_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                 effect_size_error_maaslin3(tmp_truth, tmp_associations),
                 effect_size_correlation_maaslin3(tmp_truth, tmp_associations),
                 pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                 weighted_pval_auc_maaslin3(tmp_truth, tmp_associations, abundance, metadata),
                 issue_prop(associations[associations$associations == 'prevalence',], tool))
    new_row <- lapply(new_row, function(x) {x})
    names(new_row) <- metrics
    
    new_row[['tool']] <- tool
    new_row[['iter']] <- i
    new_row[['association_type']] <- 'prevalence'
    new_row <- c(new_row, param_list)
    maaslin3_results_df <- plyr::rbind.fill(maaslin3_results_df, data.frame(new_row, check.names = F))
  }
  }
}
overall_result_df <- rbind(overall_result_df, general_results_df)
            }
          }
        }
      }
    }
  }
}

overall_result_df$tool <- gsub("fil|BY_", "", overall_result_df$tool)



save(overall_result_df, file = "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/maaslin3_benchmark/Result/DS_binary.RData")


overall_result_df$tool <- factor(overall_result_df$tool)
levels(overall_result_df$tool) <- c("BY_50_DSBin", "BY_100_DSBin", "BY_200_DSBin")

precisionPlot <- ggplot(overall_result_df, aes(x = tool, y = Precision)) +
  geom_boxplot(
    aes(fill = tool),
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ nSubjects) +  # Separate facet for each main tool group
  labs(x = "Number of Subjects", y = "Precision", title="Precision(1 - FDR) plot of different numbers of splitting at different sample sizes") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
  ) + 
  guides(fill=guide_legend(title="Number of Splitting"))

ggsave(
  filename = "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/maaslin3_benchmark/Plots/precisionPlot_DSBin.pdf",  # The file name, including extension
  plot = precisionPlot,            # The plot to save (optional if it's the last plot created)
  width = 14,                 # Width in inches
  height = 8,                # Height in inches
  dpi = 300                  # Resolution in dots per inch
)

recallPlot <- ggplot(overall_result_df, aes(x = tool, y = Recall)) +
  geom_boxplot(
    aes(fill = tool),
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ nSubjects) +  # Separate facet for each main tool group
  labs(x = "Number of Subjects", y = "Recall", title="Recall(Power) plot of different numbers of splitting at different sample sizes") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
  ) + 
  guides(fill=guide_legend(title="Number of Splitting"))


ggsave(
  filename = "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/maaslin3_benchmark/Plots/recallPlot_DSBin.pdf",  # The file name, including extension
  plot = recallPlot,            # The plot to save (optional if it's the last plot created)
  width = 14,                 # Width in inches
  height = 8,                # Height in inches
  dpi = 300                  # Resolution in dots per inch
)
