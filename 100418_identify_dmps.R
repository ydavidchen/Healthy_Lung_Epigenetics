# Differential Methylation Analysis by Lung Lobes
# Author: David Chen
# Last update: 10/04/2018
# Notes:
# -- Run this script as a program in Linux: Rscript <script-name>

rm(list=ls());
if(interactive()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");
library(limma);
library(matrixStats);
REDUCED_VAR_THRESH <- VAR_THRESH / 2; #for differential analysis
EXPORT_PATH <- paste0(DATA_PATH, "limma_DMP_100418/");

constructCustomDesignMatrix <- function(mVals, covarMat) {
  #'@description Generate design matrix & export factor block for paired study design
  #'@param mVals Matrix of M values with row=CpG, col=samples
  #'@param covarMat Covariate table where each row is a subject or sample
  
  myDesign <- model.matrix( ~ LOBE + AGE + GENDER + CellType_1, data=covarMat); 
  stopifnot(identical(covarMat$Sample_Name, colnames(mVals))); #checkpoint
  return(myDesign);
}

fitCustomLimma <- function(mVals, myDesign, myBlock, dupCorrCoef, cpgAnnot, mainOutcome) {
  #'@description Apply limma based on design matrix, pair block & internal correlation
  #'@param mVals Matrix of M values with row=CpG, col=samples
  #'@param myDesign Design matrix
  #'@param myBlock Vector indicating pairwise subject variable. Must match rows of `myDesign`
  #'@param dupCorrCoef Duplicate correlation coefficient calculated using `limma::duplicateCorrelation`
  #'@param cpgAnnot Illumina annotation data.frame
  #'@param mainOutcome String indicating main outcome of interest in design matrix
  
  fit <- lmFit(
    mVals,
    design = myDesign, 
    correlation = dupCorrCoef,
    block = myBlock
  );
  fit <- eBayes(fit);
  DMPs <- topTable(
    fit,
    number = Inf, 
    coef = mainOutcome, 
    genelist = cpgAnnot, 
    adjust.method = "fdr", 
    sort.by = "p"
  );
  return(DMPs);
}

Main <- function() {
  print("*********************** Process Begins ***********************");
  
  ## Load data & annotation:
  load(paste0(DATA_PATH, "081518_NonCF_betas.RData"));
  targets <- read.csv(paste0(DATA_PATH, "NonCF_updated_sample_sheet_090118.csv"), stringsAsFactors=FALSE); #frequently updated file
  annot.850kb3 <- loadEPICannotationFile();
  
  ## Select most variable CpGs, match samples against annotation, & convert to M-values:
  png(paste0(EXPORT_PATH,"CpG_universe_variance.png"), height=8.27, width=11.69, units="in", res=200);
  k <- decideNumberOfMostVariable(allHealthyBetas, varThresh=REDUCED_VAR_THRESH, plot=TRUE);
  dev.off();
  print(paste("Number of CpGs for differential methylation analysis:", k)); 
  betaVals <- selectMostVariableCpGs(allHealthyBetas, k=k);
  if(! identical(colnames(betaVals), targets$Sample_Name)) {
    betaVals <- betaVals[ , match(targets$Sample_Name, colnames(betaVals))];
  }
  mVals <- minfi::logit2(betaVals);
  stopifnot(identical(colnames(mVals), targets$Sample_Name)); #checkpoint
  
  ## Construct elements for LIMMA:
  myBlock <- targets$Subject; #block for calculating duplicate correlation for limma
  print(myBlock);
  
  myDesign <- constructCustomDesignMatrix(mVals, targets);
  print(myDesign);
  
  ## Calculate duplicate correlation using BETA values:
  dupCor <- limma::duplicateCorrelation(betaVals, myDesign, block=myBlock);
  consensusCorr <- dupCor$consensus.correlation;
  print(paste("Inter-sample correlation using universe beta-values:", consensusCorr));
  
  ## Build LIMMA linear model fit object:
  ann850KSub <- annot.850kb3[match(rownames(mVals),annot.850kb3$Name), ]; #subset annotation
  DMPs <- fitCustomLimma(mVals, myDesign, myBlock, consensusCorr, ann850KSub, mainOutcome="LOBERUL"); 
  DMPs$isSignifAt0.05 <- DMPs$adj.P.Val < FDR_THRESH;
  DMPs$isSignifAt0.10 <- DMPs$adj.P.Val < 2*FDR_THRESH;
  
  print("Proportion meeting FDR threshold of:", FDR_THRESH);
  print(mean(DMPs$isSignifAt0.05)); 
  print(sum(DMPs$isSignifAt0.05));
  
  print("Proportion meeting FDR threshold of:", 2*FDR_THRESH);
  print(mean(DMPs$isSignifAt0.10)); 
  print(sum(DMPs$isSignifAt0.10));
  
  ## Export:
  write.csv(DMPs, file=paste0(EXPORT_PATH, "100418_NonCF_DMPs.csv"), row.names=FALSE, quote=FALSE);
  
  print("*********************** Process Complete! ***********************");
}

if(! interactive()) Main();
