if(! require("pacman")) install.packages("pacman") 
suppressPackageStartupMessages(library("pacman"))
pacman::p_load('MCMCpack', 'glmnet', 'tidyverse', 'ape', 'GUniFrac', 'MiRKAT', 'vegan', 'parallel',
               'doParallel', 'stringr', 'dirmult', "optparse", "pkgmaker")

# Command Line Usage
option_list = list(
  make_option(
    c("-n", "--methodName"), default="CATSplit", # Which set of simulation data template
    type = "character"),
  make_option(
    c("-z", "--nReps"), default=50, # Which set of simulation data template
    type = "integer"),
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
    c("-a", "--addedCount"), default=100, # Spike-in amount
    type = "numeric"),
  make_option(
    c("-l", "--diriSum"), default=62,
    type = "numeric"),
  make_option(
    c("-m", "--metric"), default="euclidean", # Spike-in amount
    type = "character"),
  make_option(
    c("-p", "--nPerm"), default=1, # p stands for how many permutations for each split data
    type = "integer"),
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

methodName <- opt$options$methodName # Low-level parameter
nReps <- opt$options$nReps # Low-level parameter
dataSet<- opt$options$dataSet # High-level parameter
featureSet <- opt$options$featureSet # High-level parameter
sampleSize<-opt$options$sampleSize # High-level parameter
addedCount <- opt$options$addedCount # High-level parameter
diriSum <- opt$options$diriSum # Low-level parameter
metric <- opt$options$metric # Low-level parameter
nPerm <- opt$options$nPerm # Low-level parameter
parallel <- opt$options$parallel # Low-level parameter
nIterations<-opt$options$nIterations # Low-level parameter
rSeed <- opt$options$rSeed # Low-level parameter
nCore <- opt$options$nCore # Low-level parameter
workingDirectory <- opt$options$workingDirectory # Default parameter


# For Testing
methodName <- "DS"
nReps <- 50
dataSet<- "SCI" # High-level parameter
featureSet <- "mostAbun"
sampleSize<- 200 # High-level parameter
addedCount <- 200 # High-level parameter
diriSum <- 62
metric <- "euclidean"
parallel <- TRUE # Low-level parameter
nIterations<- 200 # Low-level parameter
rSeed <- 1234 # Low-level parameter
nCore <- 12 # Low-level parameter
nPerm <- 1 # Low-level parameter
workingDirectory <- "C:/Users/lclce/Desktop/2nd Sem/Yushu RA/CATSplit/DS_fdr/ScienceData_benchmark/poisson" # Default parameter

## Check Input
if(!methodName %in% c("CATSplit", "MannWhitney", "DS")) stop("methodName must be one of CATSplit, MannWhitney, DS")
if(!dataSet %in% c("SCI", "IBD")) stop("dataSet must be one of SCI, IBD")
if(!featureSet %in% c("lowCor", "highCor", "leastAbun", "mostAbun", "leastVar", "mostVar")) stop("featureSet must be one of lowCor, highCor, leastAbun, mostAbun, leastVar, mostVar")
if(!metric %in% c("euclidean", "Weighted UniFrac", "Unweighted UniFrac", "robust")) stop("metric must be one of euclidean, Weighted UniFrac, Unweighted UniFrac, robust")
if(methodName %in% c("MannWhitney", "DS")) metric <- "euclidean"
  
readData <- paste("preprocess_", dataSet, ".R", sep = "")

if(featureSet %in% c("lowCor", "highCor")){
  featureData <- "../highLowCor.RData"
} else if(featureSet %in% c("leastAbun", "mostAbun")){
  featureData <- "../mostLeastAbun.RData"
} else featureData <- "../mostVarLeastVar.RData"

## Import Needed Function
setwd(workingDirectory)
source("Split_fun.R")
source("DS_Mann_fun.R")
source("helper_function.R")
load(featureData)
spiked_chain_genera <- get(featureSet)

source(readData)
rm(meta_table)
rm(otu_table)

## Initialize Metadata
meta_table <- data.frame(BinOutcomes = c(rep(1000, sampleSize / 2), rep(0, sampleSize / 2)))
rownames(meta_table) <- c(paste("Responder", 1:(sampleSize / 2)), paste("non-Responder", 1:(sampleSize / 2)))


## Create Input Dir If Doesn't Exist
inputDirectory <- file.path(workingDirectory, 'Input')
if (!dir.exists(inputDirectory)){
  print("Input directory created!")
  dir.create(inputDirectory)
} else {
  print("Input directory already exists!")
}


## Load Input Data
inputParam <- paste(dataSet, "_", featureSet, "_s", sampleSize, "_a", addedCount, "_l", diriSum, sep = "")
sim_data_file <- paste("Input/sim_data_", inputParam, ".RData", sep = "")
if(!file.exists(sim_data_file)) stop('The input file does not exist. Generate the dataset first.')
load(sim_data_file)




# Create Output Directory
outputDirectory <- file.path(workingDirectory, 'savedData')
if (!dir.exists(outputDirectory)){
  print("Output directory created!")
  dir.create(outputDirectory)
  dir.create(paste(outputDirectory, "CATSplit_Internal", sep = "/"))
  dir.create(paste(outputDirectory, "Result", sep = "/"))
} else {
  print("Output directory already exists!")
}


res_data_file <- paste("savedData/Result/res_data_", methodName, "_", inputParam, "_", metric, "_p", nPerm, "_nReps", nReps, ".RData", sep = "")

cat(c("Job started at:",date()), "\n")
start.time <- Sys.time()

if(!file.exists(res_data_file)){
  if(methodName == "CATSplit"){
    
    cl <- makeCluster(nCore)
    registerDoParallel(cl)
    
    sim_res <- foreach(iter=1:nIterations, .combine = "rbind",
                       .packages=c("ape","GUniFrac","vegan","tidyverse", "MiRKAT"),
                       .export = c("compDist", "computeR2", "find_tau")) %dopar% {
                         
                         otu_table <- as.data.frame(sim_data[[iter]])
                         if(metric == "euclidean") otu_table <- otutabletoabundance(otu_table)
                         
                         selected_features <- CATSplit_noParallel(otu_table, taxonomy_table, meta_table, tree_data, 
                                                                  metric, nReps, qval_bound = 0.05, inputParam, iter, nPerm = nPerm)
                         
                         fd_pw <- fdp_power(selected_features, spiked_chain_genera)
                         num_selected <- length(selected_features)
                         
                         return(c(fd_pw$fdp, fd_pw$power, num_selected))
                       }
    
    stopCluster(cl)

    save(sim_res, file = res_data_file)
  }
  
  if(methodName == "MannWhitney"){
    
    cl <- makeCluster(nCore)
    registerDoParallel(cl)
    
    sim_res <- foreach(iter=1:nIterations, .combine = "rbind",
                       .packages=c("ape","GUniFrac","vegan","tidyverse", "MiRKAT")) %dopar% {
                         
                         otu_table <- otutabletoabundance(as.data.frame(sim_data[[iter]]))
                         selected_features <- Mann_WhitU(as.data.frame(sim_data[[iter]]), taxonomy_table, meta_table, qval_bound = 0.05)
                         
                         fd_pw <- fdp_power(selected_features, spiked_chain_genera)
                         num_selected <- length(selected_features)
                         
                         return(c(fd_pw$fdp, fd_pw$power, num_selected))
                       }
    
    stopCluster(cl)
    
    save(sim_res, file = res_data_file)
  }
  
  if(methodName == "DS"){
    
    cl <- makeCluster(nCore)
    registerDoParallel(cl)
    
    sim_res <-foreach(iter=1:nIterations, .combine = "rbind",
                      .packages=c("ape","GUniFrac","vegan","tidyverse", "MiRKAT", "glmnet")) %dopar% {
                        
                        otu_table <- otutabletoabundance(as.data.frame(sim_data[[iter]]))
                        selected_features <- DSBin(X = otu_table, y = meta_table$BinOutcomes, taxonomy_table, num_split = 50, qval_bound = 0.05)
                        
                        fd_pw <- fdp_power(selected_features, spiked_chain_genera)
                        num_selected <- length(selected_features)
                        
                        return(c(fd_pw$fdp, fd_pw$power, num_selected))
                      }
    
    stopCluster(cl)
    
    colnames(sim_res) <- c("FDR", "Power", "Number of Selected Features")
    
    save(sim_res, file = res_data_file)
  }
} else{
  print("Output file already exists. No new results generated!")
  load(res_data_file)
}


# Track End Time
stop.time <- Sys.time()
time<-round(difftime(stop.time, start.time, units="min"),3)
cat(c("Job finished at:",date()), "\n");
cat("Computational time:",time,"minutes \n")
cat("The output is in:", paste(outputDirectory, "Result", sep = "/"), fill=TRUE)

cat("Mean FDR:", mean(sim_res[,1]))
cat("Mean Power:", mean(sim_res[,2]))
