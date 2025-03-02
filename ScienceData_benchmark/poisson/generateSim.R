generate_sim <- function(iter, testTaxon, addedCount, taxonomy_table, otu_table, tree, averageCount, marginal, sampleSize, diriSum){
  
  set.seed(iter)
  OTUsim <- matrix(0, ncol = sampleSize, nrow = nrow(otu_table))
  rownames(OTUsim) <- rownames(otu_table)

  for (i in 1:ncol(OTUsim)) {
    OTUsim[, i] <- rmultinom(1, averageCount, rdirichlet(1, marginal * diriSum))
    if (addedCount > 0) {
      if (i <= (sampleSize / 2)) {
        for(j in 1:length(testTaxon)){
          OTUsim[taxonomy_table$taxa_to_genus %in% testTaxon[j], i] <- OTUsim[taxonomy_table$taxa_to_genus %in% testTaxon[j], i] + as.vector(rmultinom(1, size = addedCount, prob = rep(1, sum(taxonomy_table$taxa_to_genus %in% testTaxon[j]))))
        }
      }
    }
  }
  
  
  colnames(OTUsim) <- c(paste("Responder", 1:(sampleSize / 2)), paste("non-Responder", 1:(sampleSize / 2)))
  
  return(OTUsim)
}

