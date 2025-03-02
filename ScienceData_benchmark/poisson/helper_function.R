otutabletoabundance <- function(otutable){
  total_counts <- colSums(otutable)
  relative_abundance <- sweep(otutable, 2, total_counts, "/")
  relative_abundance
}

permuteRows <- function(mat) {
  # If matrix has only one row, transpose it to ensure apply works properly
  if(is.vector(mat)){
    sample(mat)
  }else{
    # Apply sample function to each row
    t(apply(mat, 1, sample))
  }
}

## compDist function
compDist<-function(otutable,metric,tree=NULL){
  if(metric %in% c("Weighted UniFrac","Unweighted UniFrac","robust")){
    temp<- GUniFrac(t(otutable),tree, verbose = FALSE)
    if(metric=="Weighted UniFrac"){
      distMat<-as.dist(temp$unifracs[,,"d_1"])         
    }else if(metric=="Unweighted UniFrac"){
      distMat<-as.dist(temp$unifracs[,,"d_UW"])   
    }else{
      distMat1<-as.matrix(as.dist(temp$unifracs[,,"d_UW"]))
      distMat1<-distMat1/max(distMat1)
      distMat2<-as.matrix(vegdist(t(otutable),"bray"))
      distMat2<-distMat2/max(distMat2) 
      distMat<-0.5*distMat1+0.5*distMat2
    }
  }else{
    distMat<-vegdist(t(otutable),metric)    
  }
  return(distMat)
}


# computeR2 <- function(testList, otutable, taxonomy, metaData, tree, metric){
#   distMat <- suppressMessages(compDist(otutable, metric, tree))
#   distResult <- adonis2(distMat ~ metaData[,"BinOutcomes"], permutations = 1)
#   origR2 <- distResult$R2[1]
#     testResult <- NULL
#     for(testIter in 1:length(testList)) {
#       otutemp <- otutable
#       otuUnder <- rownames(otutable)[taxonomy_table$taxa_to_genus %in% testList[testIter]]
#       otutemp[otuUnder,]<-permuteRows(otutable[otuUnder,])
#       distMat_removed <- compDist(otutemp, metric, tree)
#       distResult_removed <- adonis2(distMat_removed ~ metaData[,"BinOutcomes"], permutations = 1)
#       taxaR2 <- distResult_removed$R2[1]
#       deltaR2 <- origR2-taxaR2
#       newrow <- c(origR2, taxaR2, deltaR2)
#       testResult <- rbind(testResult, newrow)
#     }
# 
#   rownames(testResult) <- testList
#   colnames(testResult) <- c("original R2", "After remove taxon R2", "deltaR2")
#   return(testResult)
# }


computeR2 <- function(testList, otutable, taxonomy_table, metaData, tree, metric, nPerm = 10) {
  
  # 1. Original R squared
  distMat <- suppressMessages(compDist(otutable, metric, tree))
  distResult <- adonis2(distMat ~ metaData[,"BinOutcomes"], permutations = 1)
  origR2 <- distResult$R2[1]
  
  testResult <- NULL
  
  # Multiple permutations
  for (testIter in seq_along(testList)) {
    # Find OTUs that belong to the tested taxon
    otuUnder <- rownames(otutable)[taxonomy_table$taxa_to_genus %in% testList[testIter]]
    
    # Store permuted R squared
    permR2_vals <- numeric(nPerm)
    
    for (i in seq_len(nPerm)) {
      otutemp <- otutable
      
      # Permute target taxon
      otutemp[otuUnder, ] <- permuteRows(otutable[otuUnder, ])
      
      # Compute Distance Matrix
      distMat_removed <- compDist(otutemp, metric, tree)
      distResult_removed <- adonis2(distMat_removed ~ metaData[,"BinOutcomes"], permutations = 1)
      permR2_vals[i] <- distResult_removed$R2[1]
    }
    
    # Take mean of multiple permutations
    taxaR2_mean <- mean(permR2_vals)
    deltaR2_mean <- origR2 - taxaR2_mean
    
    # Store Result
    newrow <- c(origR2, taxaR2_mean, deltaR2_mean)
    testResult <- rbind(testResult, newrow)
  }
  
  # Output
  rownames(testResult) <- testList
  colnames(testResult) <- c("original R2", "avg R2 after removal", "avg deltaR2")
  return(testResult)
}


# mm <- M
# ww <- abs(M)


analys <- function(mm, ww, qval_bound = 0.05){
  ### mm: mirror statistics
  ### ww: absolute value of mirror statistics
  ### qval_bound:  FDR control level
  cutoff_set <- max(ww)
  for(t in ww){
    ps <- length(mm[mm >= t])
    ng <- length(na.omit(mm[mm <= -t]))
    rto <- (ng + 1)/max(ps, 1)
    if(rto <= qval_bound){
      cutoff_set <- c(cutoff_set, t)
    }
  }
  cutoff <- min(cutoff_set)
  selected_index <- which(mm >= cutoff)
  
  return(selected_index)
}



### calculate fdp and power
fdp_power <- function(selected_index, signal_index){
  num_selected <- length(selected_index)
  tp <- length(intersect(selected_index, signal_index))
  fp <- num_selected - tp
  fdp <- fp / max(num_selected, 1)
  power <- tp / length(signal_index)
  return(list(fdp = fdp, power = power))
}


find_tau <- function(M, target_fdr = 0.05) {
  # Sort the vector M in descending order
  M_abs_sorted <- sort(abs(M), decreasing = TRUE)
  
  # Initialize variables to store the best tau and current FDP
  best_tau <- NA
  fdp <- c(1)  # Start with an FDP of 1 to find the first one below target_fdr
  
  # Loop over the sorted values of M to consider them as tau candidates
  for (i in 1:length(M_abs_sorted)) {
    tau <- M_abs_sorted[i]
    
    # Calculate the number of false discoveries (numerator)
    numerator <- sum(M < -tau)
    
    # Calculate the number of positives (denominator)
    denominator <- sum(M > tau)
    
    # Calculate the False Discovery Proportion (FDP)
    fdp[i] <- numerator / max(denominator, 1)  # Avoid division by zero
    
  }
  
  # Find the largest index where FDP is below the target FDR
  valid_indices <- which(fdp <= target_fdr)
  
  
  # # Return the best tau found that makes FDP < target_fdr
  if (length(valid_indices) == 0) {
    return(1)  # If no valid tau found, return 1
  } else {
    best_index <- max(valid_indices)  # Largest index where FDP < target_fdr
    best_tau <- unname(M_abs_sorted[best_index])  # Corresponding tau
    return(best_tau)
  }
}