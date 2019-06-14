# Export Data for Gene Expression Omnibus (GEO)
# Script author: David Chen
# Last update: 06/09/2019
# Copyright (c) 2018-19 ydavidchen & Christensen-Lab
# Notes:
# -- Prepare for GEO submission upon manuscript acceptance

rm(list=ls());
if(interactive()) setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");

## CONSTANTS:
OUTPUT_DIR <- "../geo_submissions/"; #gitignore'd
NUM_NEWLY_MEASURED <- (13-1) * 2;

main <- function() {
  load(paste0(DATA_PATH, "081518_NonCF_betas.RData"));
  targets <- read.csv(paste0(DATA_PATH, "NonCF_updated_sample_sheet_090118.csv"), stringsAsFactors=FALSE); #frequently updated file
  
  print("Exporting processed beta-value matrix...");
  write.csv(allHealthyBetas, file=paste0(OUTPUT_DIR,"processed_beta_matrix.csv"), row.names=TRUE, quote=FALSE);
  
  print("Preparing sample sheet/record...");
  write.csv(targets, file=paste0(OUTPUT_DIR, "sample_sheet_for_geo.csv"), row.names=FALSE, quote=FALSE);

  ## Check IDAT files:
  ## Error is thrown if any checkpoint fails
  print("Checking IDAT files...");
  setwd("../geo_submissions/IDAT_files/"); 
  idats <- list.files(); 
  uIdats <- unique(substr(idats, 1, 19));
  stopifnot(length(uIdats) == NUM_NEWLY_MEASURED) #12 newly measured subjects, 2 samples/subject
  stopifnot(sum(targets$Complete_Barcode %in% uIdats) == NUM_NEWLY_MEASURED);
  stopifnot(all(uIdats %in% targets$Complete_Barcode));
}

if(! interactive()) main();
