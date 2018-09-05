# Assembling minfi sample sheet for QC
# Script author: David Chen
# Date: 08/15/2018; 08/17/2018
# Notes:

rm(list=ls());
library(gdata);
COVAR_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/NONCF_COVARIATES/";

Main <- function() {
  setwd(COVAR_PATH);
  
  ## Sheet for array:
  arrayInfo <- read.xls("1506(24)+1457(1)_combined.xls", stringsAsFactors=FALSE); #original file
  arrayInfo <- arrayInfo[ , 1:8];
  arrayInfo$Plate <- NULL;
  colnames(arrayInfo) <- gsub(".", "_", colnames(arrayInfo), fixed=TRUE);
  colnames(arrayInfo) <- gsub("Terminus", "Sentrix_position", colnames(arrayInfo));
  
  ## Clinical information:
  patientInfo <- read.csv("covariates_by_email_26NonCF_081718.csv", stringsAsFactors=FALSE); #assembled from email
  patientInfo$GENDER <- gsub("[0-9]", "", patientInfo$GENDER);
  
  ## Join & reformat based on minfi requirement:
  sampSheet <- merge(arrayInfo, patientInfo, by="SAMPLE_ID");
  colnames(sampSheet) <- gsub("SAMPLE_NAME", "Subject", colnames(sampSheet), ignore.case=FALSE);
  colnames(sampSheet) <- gsub("SAMPLE_ID", "Sample_ID", colnames(sampSheet), ignore.case=FALSE);
  colnames(sampSheet) <- gsub("Bead_chip__", "Bead_Chip", colnames(sampSheet));
  sampSheet$Sample_Name <- paste(sampSheet$Subject, sampSheet$LOBE, sep="_");
  sampSheet$Sample_Name <- gsub("-", "", sampSheet$Sample_Name, fixed=TRUE);
  
  ## Export:
  write.csv(sampSheet, "../IDAT_FILES/081718_minfi_sample_sheet.csv", row.names=FALSE, quote=FALSE);
  print("*********************Process Complete!*********************");
}

if(! interactive()) Main();

