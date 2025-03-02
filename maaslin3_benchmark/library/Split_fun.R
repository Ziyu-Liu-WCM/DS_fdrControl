singleSplit<-function(rep, inclusion_table, otu_table, taxonomy_table, meta_table, tree_data, 
                      metric,  qval_bound = 0.05, nPerm = 1){
  
  data1_ind <- sample(colnames(otu_table), ncol(otu_table)/2, replace = FALSE)
  data2_ind <- setdiff(colnames(otu_table), data1_ind)
  
  data1_otu_table <- otu_table[,data1_ind]
  data2_otu_table <- otu_table[,data2_ind]
  
  data1_meta_table <- meta_table[data1_ind,,drop=FALSE]
  data2_meta_table <- meta_table[data2_ind,,drop=FALSE]
  
  testList<-unique(taxonomy_table$taxa_to_genus)
  
  temp1<- computeR2(testList = testList, otutable = data1_otu_table, taxonomy = taxonomy_table, metaData = data1_meta_table,
                     tree = tree_data, metric = metric, nPerm = 1)
  R21<-temp1[1,1]
  beta1<-temp1[,3]
  
  temp2 <- computeR2(testList = testList, otutable = data2_otu_table, taxonomy = taxonomy_table, metaData = data2_meta_table,
                     tree = tree_data, metric = metric, nPerm = 1)
  R22<-temp2[1,1]
  beta2<-temp2[,3]
  
  M <- sign(beta1 * beta2) * (abs(beta1) * abs(beta2))
  selected_index <- analys(M, abs(M), qval_bound)
  M_selected <- M[selected_index]
  inclusion_rep <- ifelse(inclusion_table$feature %in% names(M_selected), 1, 0)
    
  return(list(
    inclusion_rep = inclusion_rep,
    R21 = R21,
    R22 = R22,
    beta1 = beta1,
    beta2 = beta2
  ))
}


CATSplit_Parallel <- function(otu_table, taxonomy_table, meta_table, tree_data, 
                              metric, nCore, nReps, qval_bound = 0.05, inputParam, iter = 1, nPerm = 1){
  
  if(is.null(tree_data) & metric != "euclidean") stop("Only euclidean distance is available when phylogenetic tree is not provided!")
  if(is.null(taxonomy_table)){
    cat("taxonomy_table not provided. Running per-feature test!\n")
    taxonomy_table <- data.frame(taxa_to_genus = rownames(otu_table))
  }
  
  inclusion_table <- data.frame(feature = unique(taxonomy_table$taxa_to_genus))
  
  # Check if saved data for the given inputParam exists
  save_dir <- paste0("savedData/CATSplit_Internal/", inputParam, "_", metric, "_p", nPerm, "_nReps", nReps)
  save_file <- paste0(save_dir, "/", inputParam, "_", metric, "_p", nPerm, "_nReps", nReps, "_iter", iter, ".RData")
  
  if (file.exists(save_file)) {
    # If the file exists, load the saved data
    cat("Saved R21, R22, beta1, and beta2 detected. Loading data from", save_file, "\n")
    load(save_file)  # Loads R21s, R22s, beta1s, beta2s
    
    results <- matrix(0, nrow = nrow(inclusion_table), ncol = nReps)
    
    for (rep in 1:nReps) {
      beta1 <- beta1s[[rep]]
      beta2 <- beta2s[[rep]]
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      selected_index <- analys(M, abs(M), qval_bound)
      M_selected <- M[selected_index]
      inclusion_rep <- ifelse(inclusion_table$feature %in% names(M_selected), 1, 0)
      results[, rep] <- inclusion_rep
    }
    
  } else {
    # If the file does not exist, run the parallel computation
    cat("No saved data detected. Running singleSplit in parallel...\n")
    # Register parallel backend
    cl <- makeCluster(nCore)
    registerDoParallel(cl)
    
    # Data Splitting
    start_time <- Sys.time()
    
    # Parallel
    comp_results <- foreach(rep = 1:nReps,
                            .packages = c("ape", "vegan", "GUniFrac", "doParallel"),
                            .export = c("computeR2","permuteRows", "compDist", "analys", "find_tau","singleSplit")) %dopar% {
                              
                              singleSplit(rep,inclusion_table,otu_table, taxonomy_table, meta_table, tree_data, metric, qval_bound, nPerm = 1)
                            }
    
    
    end_time <- Sys.time()
    cat("Run time:", end_time - start_time, "\n")
    # Stop the parallel backend
    stopCluster(cl)
    
    # Deregister the parallel backend
    registerDoSEQ()
    
    # Remove objects
    rm(cl)
    
    # Call garbage collection
    gc()
    
    results <- sapply(comp_results, function(x) x$inclusion_rep)
    R21s <- sapply(comp_results, function(x) x$R21)
    R22s <- sapply(comp_results, function(x) x$R22)
    beta1s <- lapply(comp_results, function(x) x$beta1)
    beta2s <- lapply(comp_results, function(x) x$beta2)
    
    if(!is.null(inputParam)){
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      save(R21s, R22s, beta1s, beta2s, file = save_file)
      cat("Saved data to", save_file, "\n")
    }
  }
  
  ## Calculated the inclusion rate
  column_sums <- ifelse(colSums(results)==0, 1, colSums(results))
  inclusion_table$inclusion_rate <- rowMeans(sweep(results, 2, column_sums, "/"))
  
  inclusion_rate <- inclusion_table$inclusion_rate
  names(inclusion_rate) <- inclusion_table$feature
  inclusion_rate <- inclusion_rate[order(inclusion_rate, decreasing = FALSE)]
  
  
  ## find the largest l
  cum_sum <- cumsum(inclusion_rate)
  l <- max(which(cum_sum <=qval_bound))
  
  ## select features and order them in descending order
  selected_features <- rev(names(inclusion_rate[-(1:l)]))

  return(selected_features)
}


CATSplit_noParallel <- function(otu_table, taxonomy_table, meta_table, tree_data, 
                                metric, nReps, qval_bound = 0.05, inputParam, iter = 1, nPerm = 1){
  
  if(is.null(tree_data) & metric != "euclidean") stop("Only euclidean distance is available when phylogenetic tree is not provided!")
  if(is.null(taxonomy_table)){
    cat("taxonomy_table not provided. Running per-feature test!\n")
    taxonomy_table <- data.frame(taxa_to_genus = rownames(otu_table))
  }
  
  inclusion_table <- data.frame(feature = unique(taxonomy_table$taxa_to_genus))
  
  # Check if saved data for the given inputParam exists
  save_dir <- paste0("savedData/CATSplit_Internal/", inputParam, "_", metric, "_nReps", nReps)
  save_file <- paste0(save_dir, "/", inputParam, "_", metric, "_p", nPerm, "_nReps", nReps, "_iter", iter, ".RData")
  
  if (file.exists(save_file)) {
    # If the file exists, load the saved data
    cat("Saved R21, R22, beta1, and beta2 detected. Loading data from", save_file, "\n")
    load(save_file)  # Loads R21s, R22s, beta1s, beta2s
    
    results_matrix <- matrix(0, nrow = nrow(inclusion_table), ncol = nReps)
    
    for (rep in 1:nReps) {
      beta1 <- beta1s[[rep]]
      beta2 <- beta2s[[rep]]
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      selected_index <- analys(M, abs(M), qval_bound)
      M_selected <- M[selected_index]
      inclusion_rep <- ifelse(inclusion_table$feature %in% names(M_selected), 1, 0)
      results_matrix[, rep] <- inclusion_rep
    }
    
  } else {
    # If the file does not exist, run the non-parallel computation
    cat("No saved data detected. Running singleSplit in non-parallel...\n")
    
    results <- list()
    R21s <- numeric(nReps)
    R22s <- numeric(nReps)
    beta1s <- list()
    beta2s <- list()
    
    start_time <- Sys.time()
    
    # Non-parallel
    for(r in 1:nReps){
      comp_results<-singleSplit(iter,inclusion_table,otu_table, taxonomy_table, meta_table, tree_data, metric, qval_bound, nPerm = 1)
      
      results[[r]] <- comp_results$inclusion_rep
      R21s[r] <- comp_results$R21
      R22s[r] <- comp_results$R22
      beta1s[[r]] <- comp_results$beta1
      beta2s[[r]] <- comp_results$beta2
    }
    results_matrix <- do.call(cbind, results)
    
    end_time <- Sys.time()
    cat("Run time:", end_time - start_time, "\n")
    
    if(!is.null(inputParam)){
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      save(R21s, R22s, beta1s, beta2s, file = save_file)
      cat("Saved data to", save_file, "\n")
    }
  }

  ## Calculated the inclusion rate
  column_sums <- ifelse(colSums(results_matrix)==0, 1, colSums(results_matrix))
  inclusion_table$inclusion_rate <- rowMeans(sweep(results_matrix, 2, column_sums, "/"))
  
  # write.csv(inclusion_table, file = "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/inclusion_table.csv", row.names = FALSE)
  # inclusion_table <- read.csv("C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/inclusion_table.csv")
  
  inclusion_rate <- inclusion_table$inclusion_rate
  names(inclusion_rate) <- inclusion_table$feature
  inclusion_rate <- inclusion_rate[order(inclusion_rate, decreasing = FALSE)]
  
  
  ## find the largest l
  cum_sum <- cumsum(inclusion_rate)
  l <- max(which(cum_sum <= qval_bound))
  
  ## select features and order them in descending order
  selected_features <- rev(names(inclusion_rate[-(1:l)]))
  
  return(selected_features)
}


