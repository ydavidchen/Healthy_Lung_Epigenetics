# Cell Type Deconvolution of Non-CF Healthy BAL: Visualize RefFreeEWAS Objects
# Author: David Chen
# Reference: Way, Titus, Johnson, Christensen
# Date: 09/01/2018
# Copyright (c) 2018-19 ydavidchen & Christensen-Lab
# Notes:
# -- Method/Code reference: Way, Titus, Johnson, Christensen
# -- Continuation of RefFree cell-type estimation done in Linux cluster

rm(list=ls());
library(ggplot2);
library(matrixStats);
library(pheatmap);
library(RefFreeEWAS);
library(doParallel); registerDoParallel(detectCores() - 1);

setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");
DATA_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/";

plot_matrices_mu_and_Omega <- function(refFreeArray) {
  #'@description Plot Matrix Mu and Omega for every K measured
  #'@param refFreeArray Object returned by supplying the "shortened" DNAm matrix to RefFreeCellMixArray
  
  parInit <- par(); 
  par(mfrow=c(3,3), oma=c(2,0,2,0));
  
  ## Plot distribution of matrices Mu:
  for(x in 1:length(refFreeArray)){
    hist(refFreeArray[[x]]$Mu, main=paste0("K=",names(refFreeArray)[x]), xlab="", col="royalblue", border=FALSE);
  }
  mtext(paste("Distribution of matrix Mu values for", length(refFreeArray), "presumptive cell types"), outer=TRUE);

  ## Check constrain the range to unit interval and check distribution of matrices Omega:
  for(x in 1:length(refFreeArray)){
    refFreeArray[[x]]$Omega[refFreeArray[[x]]$Omega < 0] <- 0.0;
    refFreeArray[[x]]$Omega[refFreeArray[[x]]$Omega > 1] <- 1.0;
    hist(refFreeArray[[x]]$Omega, main=paste0("K=",names(refFreeArray)[x]), xlab="", col="orange", border=FALSE);
  }
  mtext(paste("Distribution of matrix Omega values for", length(refFreeArray), "presumptive cell types"), outer=TRUE);
  
  par(parInit); #reset
}

getOptimalOmega <- function(refFreeArray, optimalK) {
  #'@description Retrieve the matrix of proportions of K cell types
  #'@param refFreeArray Object returned by supplying the "shortened" DNAm matrix to RefFreeCellMixArray
  #'@param optimalK Desired number of cell types to retrieve
  
  myOmega <- as.data.frame(refFreeArray[[optimalK]]$Omega);
  colnames(myOmega) <- paste0("CellType_",colnames(myOmega));
  myOmega$Sample_Name <- rownames(myOmega);
  rownames(myOmega) <- NULL;
  return(myOmega); 
}


Main <- function() {
  load(paste0(DATA_PATH, "/RefFreeEWAS_NonCF/090118_RefFreeE2_computed_objects.RData"));
  
  ## Step 1. Determine the optimal number of cell types:
  mean_Ks <- apply(RefFree_Boots[-1, ], 2, mean, trim=0.25);
  
  optimalK <- names(mean_Ks)[which.min(mean_Ks)];
  print(optimalK);
  plot_matrices_mu_and_Omega(RefFree_Array);
  
  mean_Ks <- as.data.frame(mean_Ks);
  colnames(mean_Ks) <- "Mean_Variance"; 
  mean_Ks$K <- as.numeric(rownames(mean_Ks));
  plt <- ggplot(mean_Ks, aes(K, Mean_Variance)) +
    geom_point(stat="identity", size=2) +
    myScatterTheme;
  print(plt);
  
  ## Step 2. Obtain matrix of cell types:
  matOmega <- getOptimalOmega(refFreeArray=RefFree_Array, optimalK=optimalK);
  
  ## Step 3. Update existing sample annotation & export:
  targets <- read.csv(paste0(DATA_PATH, "NonCF_updated_sample_sheet_082018.csv"), stringsAsFactors=FALSE);
  targets <- merge(targets, matOmega, by="Sample_Name");
  write.csv(targets, file=paste0(DATA_PATH, "NonCF_updated_sample_sheet_090118.csv"), row.names=FALSE, quote=FALSE);
}

Main();

