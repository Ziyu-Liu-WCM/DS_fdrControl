
## Mann Whitney U test
Mann_WhitU <- function(otu_table, taxonomy_table, meta_table, qval_bound = 0.05){
  
  
  if(is.null(taxonomy_table)){
    cat("taxonomy_table not provided. Running per-feature test!\n")
    taxonomy_table <- data.frame(taxa_to_genus = rownames(otu_table))
  }
  
  ## Aggregate at genus level
  otu_table$taxa_to_genus <- taxonomy_table$taxa_to_genus
  
  genus_count <- otu_table %>%
    dplyr::group_by(taxa_to_genus) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
    column_to_rownames("taxa_to_genus")

  otu_table <- dplyr::select(otu_table, -taxa_to_genus)

  ## Per-feature test
  test_pval <- NA
  for(i in 1:nrow(genus_count)){
    target_count <- as.numeric(genus_count[i, ])
    pval <- wilcox.test(target_count ~ meta_table$BinOutcomes, exact = FALSE)$p.value
    test_pval[i] <- ifelse(is.na(pval), 1, pval)
  }
  
  test_qval_BH <- p.adjust(test_pval, method = "BH")
  

  test_result <- data.frame(taxa_to_genus = rownames(genus_count), wilcox_qval_BH = test_qval_BH)
  selected_genus <- test_result[test_result$wilcox_qval_BH <= qval_bound,]$taxa_to_genus
  return(selected_genus)
}

# #for testing
# X <- sim_data[[1]]
# y <- meta_table$BinOutcomes
# num_split = 50
# q = 0.05
# iter = 1

DS <- function(X, y, taxonomy_table, num_split, qval_bound = 0.05){
  
  if(is.null(taxonomy_table)){
    cat("taxonomy_table not provided. Running per-feature test!\n")
    taxonomy_table <- data.frame(taxa_to_genus = rownames(otu_table))
  }
  
  X <- as.data.frame(X)
  X$taxa_to_genus <- taxonomy_table$taxa_to_genus
  
  genus_count <- X %>%
    group_by(taxa_to_genus) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
    column_to_rownames("taxa_to_genus")
  
  X <- t(genus_count)
  
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  
  for(iter in 1:num_split){
    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    
    ### get the penalty lambda for Lasso
    cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 5)
    lambda <- cvfit$lambda.min
    ?cv.glmnet
    ### run Lasso on the first half of the data
    beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta)
    # a <- glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta
    # a@Dimnames[[1]][nonzero_index]
    nonzero_index <- which(beta1 != 0)
    
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
      
      ### calculate the mirror statistics
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      # M <- abs(beta1 + beta2) - abs(beta1 - beta2)
      selected_index <- analys(M, abs(M), qval_bound)
      
      ### number of selected variables
      if(length(selected_index)!=0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
      }
    }
  }
  colnames(inclusion_rate) <- colnames(X)
  # ### single data-splitting (DS) result
  # DS_fdp <- fdp[1]
  # DS_power <- power[1]
  
  ### multiple data-splitting (MDS) result
  inclusion_rate <- apply(inclusion_rate, 2, mean)
  
  ### rank the features by the empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  
  if(length(feature_rank)!=0){
    null_feature <- numeric()
    
    ### backtracking 
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > qval_bound){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    selected_index <- colnames(X)[setdiff(feature_rank, null_feature)]
  }
  else{
    selected_index <- character(0)
  }
  
  return(selected_index)
}



DSBin <- function(X, y, taxonomy_table, num_split, qval_bound = 0.05){
  
  if(is.null(taxonomy_table)){
    cat("taxonomy_table not provided. Running per-feature test!\n")
    taxonomy_table <- data.frame(taxa_to_genus = rownames(otu_table))
  }
  
  X <- as.data.frame(X)
  X$taxa_to_genus <- taxonomy_table$taxa_to_genus
  
  genus_count <- X %>%
    group_by(taxa_to_genus) %>%
    summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE))) %>%
    column_to_rownames("taxa_to_genus")
  
  X <- t(genus_count)
  
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  fdp <- rep(0, num_split)
  power <- rep(0, num_split)
  num_select <- rep(0, num_split)
  
  for(iter in 1:num_split){
    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    
    ### get the penalty lambda for Lasso
    cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], family="binomial", type.measure = "class", nfolds = 5)
    lambda <- cvfit$lambda.min
    if(lambda==min(cvfit$lambda)){
      lambda.sequence <- exp(seq(log(min(cvfit$lambda)*0.01), log(min(cvfit$lambda)), length.out = 100))
      cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], family="binomial", nfolds = 5,lambda=lambda.sequence)
      lambda <- cvfit$lambda.min
    }
    
    if(lambda==max(cvfit$lambda)){
      lambda.sequence <- exp(seq(log(max(cvfit$lambda)), log(max(cvfit$lambda)*100), length.out = 100))
      cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], family="binomial", nfolds = 5,lambda=lambda.sequence)
      lambda <- cvfit$lambda.min
    }
    
    ### run Lasso on the first half of the data
    beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "binomial", alpha = 1, lambda = lambda)$beta)
    # a <- glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta
    # a@Dimnames[[1]][nonzero_index]
    nonzero_index <- which(beta1 != 0)
    
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
      
      ### calculate the mirror statistics
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      # M <- abs(beta1 + beta2) - abs(beta1 - beta2)
      selected_index <- analys(M, abs(M), qval_bound)
      
      ### number of selected variables
      if(length(selected_index)!=0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
      }
    }
  }
  colnames(inclusion_rate) <- colnames(X)
  # ### single data-splitting (DS) result
  # DS_fdp <- fdp[1]
  # DS_power <- power[1]
  
  ### multiple data-splitting (MDS) result
  inclusion_rate <- apply(inclusion_rate, 2, mean)
  
  ### rank the features by the empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  
  if(length(feature_rank)!=0){
    null_feature <- numeric()
    
    ### backtracking 
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > qval_bound){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    selected_index <- colnames(X)[setdiff(feature_rank, null_feature)]
  }
  else{
    selected_index <- character(0)
  }
  
  return(selected_index)
}
