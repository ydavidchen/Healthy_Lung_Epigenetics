# RPMM on Non-CF Samples and Update Sample Annotations
# Script author: David Chen
# Date: 08/20/2018
# Copyright (c) 2018-19 ydavidchen & Christensen-Lab
# Notes:

rm(list=ls());
library(RPMM);
library(doParallel); registerDoParallel(detectCores()-1);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");

DATA_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/"; 
VAR_THRESH <- 0.01; 
MAX_LEV <- 2;

getRPMMClustLabels <- function(rpmmObject, Y_inv=NULL) {
  #'@description Extracts RPMM hard cluster labels
  #'@param rpmmObject RPMM object
  #'@param Y_inv Optional. Input matrix for RPMM computation. If provided, the sample name will be updated
  hardLabels <- blcTreeLeafClasses(rpmmObject);
  hardLabels <- as.data.frame(hardLabels);
  colnames(hardLabels) <- "RPMMClusters";
  if(! is.null(Y_inv)) rownames(hardLabels) <- hardLabels$Sample_Name <- rownames(Y_inv);
  return(hardLabels);
}

getRPMMSampOrder <- function(rpmmClusters, Y_inv) {
  #'@description Retrieves sample orders fo heat map visualization
  #'@describeIn Hinoue et al. online tutorial
  #'@param rpmmClusters data.frame with row.names = sample ID & 1 column named RPMM with cluster assignments
  #'@param Y_inv Input matrix for RPMM computation
  sampOrder <- c();
  for(r in names(table(rpmmClusters$RPMM))) {
    samps <- rownames(rpmmClusters)[rpmmClusters$RPMM == r];
    clu <- t(Y_inv[samps, ]);
    s_i <- seriation::seriate(clu, margin=2);
    so_i <- seriation::get_order(s_i);
    sampOrder <- c(sampOrder, samps[so_i]);
  }
  sampOrder <- data.frame(
    Sample_Name = sampOrder,
    RPMMSampleOrder = 1:length(sampOrder)
  );
  return(sampOrder);
}

calcHorvathAge <- function(betas) {
  #'@description Calculates Horvath methylation age
  #'@param betas Genome-wide methylation beta value matrix, row=CpGs, column=samples
  require(wateRmelon);
  DNAmAge <- wateRmelon::agep(betas);
  DNAmAge <- as.data.frame(DNAmAge);
  colnames(DNAmAge) <- "HorvathAge";
  DNAmAge$Sample_Name <- rownames(DNAmAge);
  return(DNAmAge); 
}

calcAgeAcceleration <- function(targets) {
  #'@description Calculates age acceleration: Horvath vs. age regression residual
  #'@describeIn Johnson KC et al. 2017 Breast Cancer Research
  #'@param targets Sample annotation with "HorvathAge" & "Age" columns
  #'@return None
  print("Spearman correlation analysis:");
  print( cor.test(targets$Horvath, targets$AGE, method="spearman") );

  clockRegress <- lm(HorvathAge ~ AGE, data=targets); #order of regression variable matters!
  accel <- summary(clockRegress)$coefficients[2, 4]; #exact P-value
  print(paste("Age acceleration:", accel));
}


Main <- function() {
  load(paste0(DATA_PATH, "081518_NonCF_betas.RData"));
  pdf(paste0(DATA_PATH,"082818_NonCF_RPMM.pdf"), paper="a4r"); #captures graphic outputs
  
  ## Select most variable CpGs for further exploration:
  k <- decideNumberOfMostVariable(allHealthyBetas, varThresh=VAR_THRESH, plot=FALSE);
  
  ## Task 1. RPMM: 
  Y_inv <- t(selectMostVariableCpGs(allHealthyBetas, k=k));
  stopifnot(identical(rownames(Y_inv), targets$Sample_Name)); #checkpoint
  
  rpmm <- blcTree(Y_inv, maxlevel=MAX_LEV, splitCriterion=blcSplitCriterionBIC, verbose=1);
  rpmm
  
  plotImage.blcTree(rpmm); #added to pdf
  plotTree.blcTree(rpmm); #added to pdf
  dev.off();
  
  rpmmClusters <- getRPMMClustLabels(rpmm, Y_inv);
  targets <- merge(targets, rpmmClusters, by="Sample_Name"); #add into annotation
  
  sampOrders <- getRPMMSampOrder(rpmmClusters, Y_inv);
  targets <- merge(targets, sampOrders, by="Sample_Name"); #add into annotation; for heat map ordering

  ## Task 2. Calculate Horvath's DNA methylation age:
  horvathAge <- calcHorvathAge(allHealthyBetas);
  targets <- merge(targets, horvathAge, by="Sample_Name"); #add into annotation
  
  ## Export updated sample annotations & plots:
  targets <- targets[match(colnames(allHealthyBetas), targets$Sample_Name), ];
  rownames(targets) <- NULL; 
  write.csv(targets, file=paste0(DATA_PATH, "NonCF_updated_sample_sheet_082018.csv"), row.names=FALSE, quote=FALSE);
}

Main();
