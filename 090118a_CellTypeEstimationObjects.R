# Cell Type Deconvolution of Non-CF Healthy BAL: RefFreeEWAS Objects
# Author: David Chen
# Reference: Way, Titus, Johnson, Christensen
# Date: 09/01/2018
# Notes:
# -- Upload DNAm RData, targets CSV & this script to Linux cluster (e.g. Discovery or Andes-Polaris)
# -- Use a PBS script or `nohup Rscript <this-script.R> &;`

rm(list=ls());
library(matrixStats);
library(RefFreeEWAS);
library(doParallel); registerDoParallel(detectCores() - 1);

nCpGs <- 1e4;
K_LIST <- 2:10;
ITERS <- 10;
BOOTSTRAP_ITER <- 5; #more may be necessary 
R <- 500;

setDataPath <- function(verbose=FALSE) {
  #'@description Sets path for previously saved RData file based on operating system
  myMachine <- Sys.info()["sysname"]; 
  if(myMachine == "Darwin") {
    if(verbose) print("Running on local machine...");
    dataPath <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/";
  } else if (myMachine == "Linux") {
    if(verbose) print("Running on Linux cluster...");
    dataPath <- "./data/"; 
  } else {
    stop("Check your operating system!");
  }
  
  if(verbose) {
    print(paste("Current DATA_PATH is:", dataPath ));
    print(paste("Current working directory is:", getwd() )); 
  }
  return(dataPath); 
}

selectMostVariableCpGs <- function(data, k) {
  #'@description Subset a matrix of CpGs based on sample variance
  sele <- order(matrixStats::rowVars(data), decreasing=TRUE)[1:k];
  mat <- data[sele, ];
  return(mat); 
}


Main <- function() {
  print("*********************Process Begins*********************");
  
  DATA_PATH <- setDataPath(verbose=TRUE);
  load(paste0(DATA_PATH,"081518_NonCF_betas.RData")); #loads allHealthyBetas
  targets <- read.csv(paste0(DATA_PATH,"NonCF_updated_sample_sheet_082018.csv"), stringsAsFactors=FALSE); #overrides targets
  Y_full <- allHealthyBetas; #copy
  Y_shortened <- selectMostVariableCpGs(Y_full, k=nCpGs); #subset
  
  ## Steps 1a: Alternate fixing matrices Mu & Omega by iterating from 2-10 cell types:
  RefFree_Array <- RefFreeCellMixArray(
    Y_shortened,
    Klist = K_LIST,
    iters = ITERS
  ); 
  print( sapply(RefFree_Array, deviance, Y=Y_shortened) ); 
  
  ## Step 1b. Use the full beta matrix plus the k most variable probes to infer the cell types:
  RefFree_Array2 <- RefFreeCellMixArray(
    Y_full,
    Klist = K_LIST,
    iters = ITERS,
    Yfinal = Y_shortened
  ); 
  print(sapply(RefFree_Array2, deviance, Y=Y_shortened)); 
  
  ## Step 3: Determine the optimal number of putative cell types K by bootstrapping:
  RefFree_Boots <- RefFreeCellMixArrayDevianceBoots(
    RefFree_Array2,
    Y_shortened,
    R = R, #more bootstraps may be necessary
    bootstrapIterations = BOOTSTRAP_ITER
  ); 
  
  ## Export:
  print("Saving objects...");
  save(
    list = c("RefFree_Array","RefFree_Array2","RefFree_Boots"),
    file = "090118_RefFreeE2_computed_objects.RData", 
    compress = TRUE
  ); 
  
  print(sessionInfo());
  print("*********************Process Complete!*********************");
}

Main(); 
