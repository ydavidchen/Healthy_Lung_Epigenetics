# Assembling minfi sample sheet for QC
# Script author: David Chen
# Date: 08/15/2018
# Notes:

rm(list=ls());
library(gdata);
COVAR_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/NONCF_COVARIATES/";

Main <- function() {
  setwd(COVAR_PATH);
  
  ## Sheet for array:
  arrayInfo <- read.xls("1506 (Armstrong-Ashare-24).xls", stringsAsFactors=FALSE); #original file
  arrayInfo <- arrayInfo[ , 1:8];
  arrayInfo$Plate <- NULL; 
  colnames(arrayInfo) <- gsub(".","_",colnames(arrayInfo),fixed=TRUE);
  colnames(arrayInfo) <- gsub("Terminus", "Sentrix_position", colnames(arrayInfo));
  
  ## Clinical information:
  patientInfo <- read.csv("covariates_by_email_24NonCF.csv", stringsAsFactors=FALSE); #assembled from email
  patientInfo$Sex <- gsub("[0-9]", "", patientInfo$GENDER);
  patientInfo$subjectMatch <- gsub("[A-Z]", "", patientInfo$GENDER);
  stopifnot( all(table(patientInfo$Sex) == 12) & all(table(patientInfo$subjectMatch) == 4)); #checkpoint
  patientInfo$GENDER <- NULL; 
  
  ## Join & reformat based on minfi requirement:
  sampSheet <- merge(arrayInfo, patientInfo, by="SAMPLE_ID");
  colnames(sampSheet) <- gsub("SAMPLE_NAME", "Subject", colnames(sampSheet), ignore.case=FALSE);
  colnames(sampSheet) <- gsub("SAMPLE_ID", "Sample_ID", colnames(sampSheet), ignore.case=FALSE);
  colnames(sampSheet) <- gsub("Bead_chip__", "Bead_Chip", colnames(sampSheet));
  sampSheet$Sample_Name <- paste(sampSheet$Subject, sampSheet$LOBE, sep="_");
  sampSheet$Sample_Name <- gsub("-", "", sampSheet$Sample_Name, fixed=TRUE);
  
  ## Export:
  write.csv(sampSheet, "../IDAT_FILES/081518_minfi_sample_sheet.csv", row.names=FALSE, quote=FALSE);
  print("*********************Process Complete!*********************");
}

if(! interactive()) Main();

