

load_essential_packages<-function(){
  if(! require("pacman")) install.packages("pacman", repos='http://cran.us.r-project.org')
  suppressPackageStartupMessages(library("pacman"))
  pacman::p_load("devtools")
  
  ###################
  # Core R Packages #
  ###################
  
  pacman::p_load('tidyverse', 'fdrtool', 'ashr')
  
  ########################
  # Core GitHub Packages #
  ########################
  
  if(! require("GMPR")) {
    devtools::install_github("lichen-lab/GMPR")
  }
  if(! require("swfdr")) {
    devtools::install_github("leekgroup/swfdr")
  }
  library(GMPR)
  library(swfdr)
  
  ##############################
  # Core Bioconductor Packages #
  ##############################
  
  if(! require("genefilter")) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("genefilter")
  }
  library(genefilter)
  
  if(! require("IHW")) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("IHW")
  }
  library(IHW)
}


append_qvalues<-function(paras){
  paras<- paras %>%
    dplyr::mutate(qval_BH = tryCatch(p.adjust(pval, method = 'BH'),
                                     error = function(err){NA}),
                  qval_BY = tryCatch(p.adjust(pval, method = 'BY'),
                                     error = function(err){NA}))
  
  return(paras)
}



#### Apply Filtering To A Dataset ####

SimpleFilter = function(physeq, Threshold_Abundance, Threshold_Prevalence) {
  
  features<-physeq$features
  
  ## filter columns
  filtered_features<-features[,colSums(features > Threshold_Abundance) > nrow(features)*Threshold_Prevalence] 
  
  ## filter rows
  rows_to_keep <- rowSums(filtered_features) > 0
  filtered_features <- filtered_features[rows_to_keep,]
  metadata <- as.matrix(physeq$metadata[rows_to_keep, , drop = FALSE])
  colnames(metadata) <- colnames(physeq$metadata)
  
  return(list(features = filtered_features, metadata = metadata, ID = physeq$ID, libSize=rowSums(filtered_features)))
}

#################################################################################
# Apply SimpleFilter (Prevalence and Abundance Filtering) To A List of Datasets #
#################################################################################

list.SimpleFilter<- function(simlist, Threshold_Abundance, Threshold_Prevalence){
  return(lapply(simlist, SimpleFilter, Threshold_Abundance, Threshold_Prevalence))
}

#### Apply SCRAN Normalization To A Dataset ####

SCRANnorm = function(features) {
 
  #### Extract Features ####
  
  features <- as.matrix(features)
  
  #### SCRAN Normalizing the Data ####
  
  sce <- t(features)
  sizeFactors <- scran::calculateSumFactors(sce)
  
  #### Reset library size ####
  
  libSize <- sizeFactors
  return(libSize)
}

#### CCT P-value Combination ####

CCT <- function(pvals, weights=NULL){
  #### check if there is NA
  pvals <- ifelse(is.na(pvals), 1, pvals)
  # if(sum(is.na(pvals)) > 0){
  #   stop("Cannot have NAs in the p-values!")
  # }
  
  #### check if all p-values are between 0 and 1
  if((sum(pvals<0) + sum(pvals>1)) > 0){
    stop("All p-values must be between 0 and 1!")
  }
  
  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals==0)>=1)
  is.one <- (sum(pvals==1)>=1)
  if(is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if(is.zero){
    return(0)
  }
  if(is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }
  
  #### check the validity of weights (default: equal weights) and standardize them.
  if(is.null(weights)){
    weights <- rep(1/length(pvals),length(pvals))
  }else if(length(weights)!=length(pvals)){
    stop("The length of weights should be the same as that of the p-values!")
  }else if(sum(weights < 0) > 0){
    stop("All the weights must be positive!")
  }else{
    weights <- weights/sum(weights)
  }
  
  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0){
    cct.stat <- sum(weights*tan((0.5-pvals)*pi))
  }else{
    cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
    cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
  }
  
  #### check if the test statistic is very large.
  if(cct.stat > 1e+15){
    pval <- (1/cct.stat)/pi
  }else{
    pval <- 1-pcauchy(cct.stat)
  }
  return(pval)
}

combinePvalCCT<-function(pvalue1, pvalue2){
  
  if (length(pvalue1)!=length(pvalue2)){stop('the lengths of two pvalue vectors should match up!')}
  n<-length(pvalue1)
  combined_p_value<-c()
  for(i in 1:n){
    combined_p_value[i]<-CCT(c(pvalue1[i], pvalue2[i]))
  }
  return(combined_p_value)
}

#### Calculating Performance Metrics ####

PerfMetrics <- function(x, alpha = 0.05){
  
  #### Gathering Findings ####
  
  resi = x$res
  features = x$features
  resi[is.na(resi[, "adjPval"]), "adjPval"] = 1
  wh.pos = resi[resi$adjPval < alpha, ]$feature
  wh.neg = resi[!resi$feature %in% wh.pos, ]$feature
  wh.TP = names(features[, grep("_TP", names(features))])
  FPs = sum(!wh.pos %in% wh.TP)
  TPs = sum(wh.pos %in% wh.TP)
  TNs = sum(!wh.neg %in% wh.TP)
  FNs = sum(wh.neg %in% wh.TP)
  
  #### Getting Basic Measures ####
  
  Sensitivity = TPs/(TPs + FNs)
  Specificity = TNs/(TNs + FPs)
  FDR = ifelse((TPs + FPs) == 0, 0, FPs/(TPs + FPs))
  FScore = 2*TPs/(2*TPs + FNs + FPs) 
  
  #### Matthew's Correlation Coefficient ####
  
  numerator <- (TPs * TNs - FPs * FNs)
  denominator <- sqrt(as.numeric((TPs + FPs))*as.numeric((TPs + FNs))*as.numeric((TNs + FPs))*as.numeric((TNs + FNs)))
  if(denominator == 0) denominator <- 1
  if(is.na(denominator)) denominator <- 1
  MCC <- numerator/denominator
  
  #### Getting AUC ####
  
  all.TPs = grep("_TP", rownames(features))
  wh.truth = (1:nrow(resi) %in% all.TPs)
  wh.pred = (resi[, "adjPval"] < alpha)
  if (all(wh.truth == 'TRUE')) {
    AUC = 1
    pAUC = 1
    fAUC = 1
  } else if (all(wh.truth == 'FALSE')) {
    AUC = 0
    pAUC = 0
  } else {
    pred = prediction(as.numeric(wh.pred), factor(wh.truth, levels = c("TRUE", "FALSE")))
    AUC = round(performance(pred, "auc")@y.values[[1]],3)
    pAUC = round(performance(pred, "auc", fpr.stop = 0.20)@y.values[[1]][[1]], 3)
  }
  Measures <- data.frame(Sensitivity = Sensitivity,
                         FPR = 1 - Specificity,
                         FDR = FDR,
                         FScore = FScore,
                         MCC = MCC,
                         AUC = AUC,
                         pAUC = pAUC,
                         alpha = alpha)
  return(Measures)
}

#### Getting Performance Metrics on Lists of Outputs ####

list.perfMetrics <- function(x){
  
  #### Performance Metrics ####
  
  performance <- PerfMetrics(x)
  ResTmp <- data.frame(Method = method,
                       nGenes = nGenes,
                       nCells = nCells,
                       effectSize = effectSize,
                       dropout = dropout,
                       nDE = (pDE * nGenes),
                       performance,
                       Time.min = x$Time.min,
                       Iteration = nreps)
  return(ResTmp)
}

#### Getting Summary on Outputs ####

df.summary = function(x){
  sum.res = data.frame(Method = unique(x$Method),
                       nGenes = unique(x$nGenes),
                       nCells = unique(x$nCells),
                       effectSize = unique(x$effectSize),
                       dropout = unique(x$dropout),
                       nDE = unique(x$nDE),
                       Sensitivity = median(x$Sensitivity, na.rm = TRUE),
                       FPR = median(x$FPR, na.rm = TRUE),
                       FDR = median(x$FDR, na.rm = TRUE),
                       FScore = median(x$FScore, na.rm = TRUE),
                       MCC = median(x$MCC, na.rm = TRUE),
                       AUC = median(x$AUC, na.rm = TRUE),
                       pAUC = median(x$pAUC, na.rm = TRUE),
                       alpha = unique(x$alpha),
                       Time.min = median(x$Time.min, na.rm = TRUE),
                       Iteration = unique(x$Iteration))
  return(sum.res)
}





eval_res_list = function(resi, alpha = 0.05) {
  
  # resi <- output$ANCOMBC2.TPVC_ZeroInflate_noRandomEffect_UVB_200_1_1000_0.1_1_1_5_50000_1
  qvalmethod <- c('BH', 'BY')
  
  qval <- paste('qval', qvalmethod, sep = '_')
  # Define a q-value column
  
  # Replace Missing Q-values to the Highest Possible Value 1.0
  for (i in qval){
    resi[is.na(resi[, i]), i] <- 1
  }
  
  # Evaluate Detection Performance
  time = mean(resi[,"time"], na.rm=TRUE)
  wh.pred = (resi[, qval] <= alpha)
  
  wh.pos = list()
  wh.neg = list()
  for (col in colnames(wh.pred)){
    wh.pos[[col]] = which(wh.pred[,col])
    wh.neg[[col]] = which(!wh.pred[,col])
  } 
  
  
  wh.TP = grep("[[:print:]]+\\_TP$", resi[, "feature"])
  FPs = lapply(wh.pos, function(x) sum(!x %in% wh.TP))
  TPs = lapply(wh.pos, function(x) sum(x %in% wh.TP))
  TNs = lapply(wh.neg, function(x) sum(!x %in% wh.TP))
  FNs = lapply(wh.neg, function(x) sum(x %in% wh.TP))
  
  
  
  
  # Sensitivity: True Positives Divided by All Positives (Sum of True
  # Positives and False Negatives)
  Sensitivity = mapply(function(x,y) x/(x+y), TPs, FNs, SIMPLIFY = FALSE)
  
  
  # Specificity: True Negatives Divided by All Negatives (Sum of True
  # Negatives and False Positives)
  Specificity = mapply(function(x,y) x/(x+y), TNs, FPs, SIMPLIFY = FALSE)
  
  # False Discovery Rate: False Positives Divided by All Detected Positives
  FDR = mapply(function(x,y) {if ((x + y) == 0) 
    0 else y/(x + y)}, TPs, FPs, SIMPLIFY = FALSE)
  
  # If no true positives, return NA's for irrelevant measures
  
  # FScore
  
  
  FScore <- mapply(function(TP, FN, FP) 2 * TP / (2 * TP + FN + FP), TPs, FNs, FPs)
  
  
  # Calculate Matthew's Correlation Coefficient (MCC)
  MCC <- mapply(function(TP, TN, FP, FN) {
    numerator <- (TP * TN - FP * FN)
    denominator <- sqrt(as.numeric((TP + FP)) * as.numeric((TP + FN)) * as.numeric((TN + FP)) * as.numeric((TN + FN)))
    if (denominator == 0 || is.na(denominator)) denominator <- 1
    MCC <- numerator / denominator
    return(MCC)
  }, TPs, TNs, FPs, FNs)
  
  
  
  # AUC, pAUC (FPR < 0.20)
  wh.truth = (1:nrow(resi) %in% wh.TP)
  AUC <- list()
  pAUC <- list()
  if (all(wh.truth=='TRUE')) {
    for (col in colnames(wh.pred)){
      AUC[[col]]=1
      pAUC[[col]]=1
    }
  } else if (all(wh.truth=='FALSE')) {
    for (col in colnames(wh.pred)){
      AUC[[col]]=0
      pAUC[[col]]=0
    }
  } else {
    for (col in colnames(wh.pred)){
      pred <- prediction(as.numeric(wh.pred[,col]), factor(wh.truth, levels=c("TRUE", "FALSE")))
      AUC[[col]] = performance(pred, "auc")@y.values[[1]]
      pAUC[[col]] = performance(pred, "auc", fpr.stop=0.20)@y.values[[1]][[1]]
    } 
  }
  
  # Departure from Uniformity Under the Null
  AucAocVals <- AucAocFun(resi[!wh.truth, "pval"], plotIt = FALSE, pch = 20, type = "l")
  totalArea = AucAocVals["conservArea"] + AucAocVals["liberalArea"]
  names(totalArea) = "totalArea"  
  
  
  # Return
  combine_results <- function(Sensitivity, Specificity, FDR, FScore, MCC, AUC, pAUC, conservArea, liberalArea, totalArea, time, alpha, qval) {
    return(list(
      Sensitivity = Sensitivity,
      Specificity = Specificity,
      FDR = FDR,
      FScore = FScore,
      MCC = MCC,
      AUC = AUC,
      pAUC = pAUC,
      conservArea = conservArea,
      liberalArea = liberalArea,
      totalArea = totalArea,
      time = time,
      alpha = alpha,
      qval = qval
    ))
  }
  
  # Use mapply to iterate over the elements and combine the results into a list of lists
  return(mapply(function(sens, spec, fdr, fscore, mcc, auc, pauc, consArea, libArea, tArea, t, a, q) {
    combine_results(sens, spec, fdr, fscore, mcc, auc, pauc, consArea, libArea, tArea, t, a, sub(".*_", "", q))
  }, Sensitivity, Specificity, FDR, FScore, MCC, AUC, pAUC, AucAocVals["conservArea"], AucAocVals["liberalArea"], totalArea, time, alpha, qval, SIMPLIFY = FALSE)
  )
}


make_power_df = function(reslist, simparamslabels) {
  # Get the names of the sublists (assuming all elements of reslist have the same structure)
  sublist_names = names(reslist[[1]])
  
  # Initialize an empty list to store the combined results for each sublist
  combined_results = list()
  
  # Iterate over the sublist names and combine the results for each
  for (sublist_name in sublist_names) {
    # Extract the corresponding sublists from each element of reslist
    sublist = lapply(reslist, function(x) data.frame(x[[sublist_name]]))
    
    # Combine the sublist results into a data frame
    powerdf = ldply(sublist, .parallel = TRUE)
    colnames(powerdf)[1] <- "Combinations"
    paramdf = ldply(strsplit(powerdf[, "Combinations"], "_"), .parallel = TRUE)
    colnames(paramdf) <- simparamslabels
    powerdf = cbind(powerdf, paramdf)
    
    # Store the combined data frame in the results list
    combined_results[[sublist_name]] = powerdf
  }
  
  return(combined_results)
}



AucAocFun <- function(pVals, maxVal = 0.25, plotIt = FALSE, ...) {
  ## sum of height differences between the curve and the y=x line, half maximum
  ## = 0.5 * max(pVals) * length(pVals) this is the case of 100% rejection with
  ## any alpha value _onlyTotArea_ controls if only the total area is returned,
  ## and not also the Area Over and Under the y=x line.  halfMaxArea <- 0.5 *
  ## max(maxVal) * length(pVals)
  halfMaxArea <- 0.5 * length(pVals)
  
  pVals <- pVals[!is.na(pVals)]
  estimated <- sort(pVals)
  theoretic <- seq_along(pVals)/length(pVals)
  
  if (plotIt) {
    plot(theoretic, estimated, ...)
    abline(0, 1)
  } else {
  }
  
  diffPVals <- theoretic - estimated
  indConserv <- theoretic <= estimated
  conservArea <- sum(-diffPVals[indConserv])/halfMaxArea
  liberalArea <- sum(diffPVals[!indConserv])/halfMaxArea
  
  c(conservArea = conservArea, liberalArea = liberalArea, totalArea = liberalArea + 
      conservArea)
  
}  # END - function: aucAocFun, AUC/AOC calculation for p-values distribution