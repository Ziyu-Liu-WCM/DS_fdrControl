if(! require("pacman")) install.packages("pacman") 
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('MCMCpack', 'glmnet', 'tidyverse', 'ape', 'GUniFrac', 'MiRKAT', 'vegan', 'parallel',
               'doParallel', 'stringr', 'dirmult', "optparse", "pkgmaker")

# Command Line Usage
option_list = list(
  make_option(
    c("-d", "--dataSet"), # Which set of simulation data template
    type = "character"),
  make_option(
    c("-f", "--featureSet"), # Which set of simulation data template
    type = "character"),
  make_option(
    c("-s", "--sampleSize"), # Sample size of simulation data
    type = "integer"),
  make_option(
    c("-a", "--addedCount"), default=100,# Spike-in amount
    type = "numeric"),
  make_option(
    c("-l", "--diriSum"), default=62,# Spike-in amount
    type = "numeric"),
  make_option(
    c("-m", "--metric"), default="euclidean",# Spike-in amount
    type = "character"),
  make_option(
    c("-g", "--parallel"), default=TRUE, # Whether parallel computing
    action = "store_false"),
  make_option(
    c("-t", "--nIterations"), default=200, # t stands for how many times an experiment is repeated
    type = "integer"),
  make_option(
    c("-r", "--rSeed"), default=1234, # r stands for reproducibility index (random seed)
    type = "integer"),
  make_option(
    c("-c", "--nCore"), default=4, # c stands for how many cores to be used in the analysis
    type = "integer"),
  make_option(
    c("-w", "--workingDirectory"), # w stands for working directory
    type = "character"))
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
print(opt)

dataSet<- opt$options$dataSet # High-level parameter
featureSet <- opt$options$featureSet # High-level parameter
sampleSize<-opt$options$sampleSize # High-level parameter
addedCount <- opt$options$addedCount # High-level parameter
diriSum <- opt$options$diriSum # Low-level parameter
metric <- opt$options$metric # Low-level parameter
parallel <- opt$options$parallel # Low-level parameter
nIterations<-opt$options$nIterations # Low-level parameter
rSeed <- opt$options$rSeed # Low-level parameter
nCore <- opt$options$nCore # Low-level parameter
workingDirectory <- opt$options$workingDirectory # Default parameter

# For Testing
dataSet<- "SCI" # High-level parameter
featureSet <- "mostAbun"
sampleSize<- 200 # High-level parameter
addedCount <- 200 # High-level parameter
diriSum <- 62
metric <- "euclidean"
parallel <- TRUE # Low-level parameter
nIterations<- 200 # Low-level parameter
rSeed <- 1234 # Low-level parameter
nCore <- 10 # Low-level parameter
workingDirectory <- "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/DS_fdr/ScienceData_benchmark/poisson" # Default parameter

## Check Input
if(!dataSet %in% c("SCI", "IBD")) stop("dataSet must be one of SCI, IBD")
if(!featureSet %in% c("lowCor", "highCor", "leastAbun", "mostAbun", "leastVar", "mostVar")) stop("featureSet must be one of lowCor, highCor, leastAbun, mostAbun, leastVar, mostVar")
if(!metric %in% c("euclidean", "Weighted UniFrac", "Unweighted UniFrac", "robust")) stop("metric must be one of euclidean, Weighted UniFrac, Unweighted UniFrac, robust")

readData <- paste("preprocess_", dataSet, ".R", sep = "")

if(featureSet %in% c("lowCor", "highCor")){
  featureData <- "../highLowCor.RData"
} else if(featureSet %in% c("leastAbun", "mostAbun")){
  featureData <- "../mostLeastAbun.RData"
} else featureData <- "../mostVarLeastVar.RData"

## Import Needed Function
setwd(workingDirectory)
source("generateSim.R")
source("Split_fun.R")
source("helper_function.R")
source(readData)
load(featureData)


## Parameter Setting
spiked_chain_genera <- get(featureSet)
marginal<-apply(otutabletoabundance(otu_table),1,mean)
averageCount<-mean(apply(otu_table,2,sum))
testTaxon <- spiked_chain_genera


## Initialize Metadata
meta_table <- data.frame(BinOutcomes = c(rep(1000, sampleSize / 2), rep(0, sampleSize / 2)))
rownames(meta_table) <- c(paste("Responder", 1:(sampleSize / 2)), paste("non-Responder", 1:(sampleSize / 2)))

inputParam <- paste(dataSet, "_", featureSet, "_s", sampleSize, "_a", addedCount, "_l", diriSum, sep = "")
sim_data_file <- paste("Input/sim_data_", inputParam, ".RData", sep = "")

## Check Input Directory
inputDirectory <- file.path(workingDirectory, 'Input')
if (!dir.exists(inputDirectory)){
  print("Input directory created!")
  dir.create(inputDirectory)
} else {
  print("Input directory already exists!")
}


cl <- makeCluster(nCore)
registerDoParallel(cl)
if(!file.exists(sim_data_file)){
  sim_data <- foreach(seedNum=1:nIterations,
                      .packages=c("ape","GUniFrac","vegan","dirmult")) %dopar% {
                        generate_sim(seedNum, testTaxon, addedCount, taxonomy_table, otu_table, tree_data, averageCount, marginal, sampleSize, diriSum)
                      }
  save(sim_data, file = sim_data_file)
} else {
  print("Input file already exists. No new data generated!")
  load(sim_data_file)
}
stopCluster(cl)

cat("Simulation datasets at", sim_data_file)

