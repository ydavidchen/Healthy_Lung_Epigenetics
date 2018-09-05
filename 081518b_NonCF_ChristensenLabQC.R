# Functional Implementation of Christensen Lab Minfi & ENmix QC
# Script author: David Chen
# Date: 08/15/2018
# Notes:

rm(list=ls());
library(doParallel); registerDoParallel(detectCores() - 1);
library(ENmix);
library(matrixStats);
library(minfi);

## User variables:
IDAT_DIR <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/IDAT_FILES/";
OUTPUT_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/NonCF_EPIC_QC/082018_Combined_Set_QC";
DET_P <- 0.000001; #ENmix default
SAMP_FRACTION <- 0.20;
DPI <- 200;

if(Sys.getenv("RSTUDIO")=="1") setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");

## Pipelines & helper functions:
ChristensenLabMethArrayQCs <- function(rgSet, detP, sampThresh, dpi, outputPath, what2Return=c("ENmix","minfi"), usedFFPE=TRUE) {
  #'@description Phase 1 of 3: Run a set of existing MethylationEPIC QC pipelines
  #'@param rgSet RGChannelSet object of IDATs loaded by minfi
  #'@param detP Min detection P-value
  #'@param sampThresh Max proportion of samples allowing to exceed detection P
  #'@param dpi png figure resolution passed into `png` function
  #'@param outputPath Path to the output directory
  #'@param what2Return User choice of minfi or ENmix object to return
  #'@param usedFFPE Boolean indicating whether FFPE Restoration Control should be saved as a plot
  
  if(getwd() != outputPath) setwd(outputPath);
  print(paste("Current working directory has been set to:", getwd()));
  
  ## 1) Plot ENmix control plots:
  ENmix::plotCtrl(rgSet);
  
  ## 2) Run ENmix QC: 
  qc <- ENmix::QCinfo(rgSet, detPthre=detP, samplethre=sampThresh);
  
  ## 3) Generate minfi QC report:
  minfi::qcReport(rgSet=rgSet, pdf="minfi_QC_report.pdf");
  
  ## 4) minfi identification of outliers based on meth & unmeth intensities
  Mset <- minfi::preprocessRaw(rgSet);
  Mset <- minfi::minfiQC(Mset, fixOutliers=TRUE, verbose=TRUE);
  png("minfi_OutliersByIntensity.png", height=8.27, width=11.69, units="in", res=dpi);
  minfi::plotQC(Mset$qc);
  title("Poor-performing Outlier Identification");
  dev.off();
  
  ## 6) Plot FFPE restoration probe:
  if(usedFFPE) {
    png("minfi_FFPEcontrol.png", height=8.27, width=11.69, units="in", res=dpi);
    controlStripPlot(rgSet, controls=c("RESTORATION"));
    dev.off();
  }
  
  if(what2Return == "ENmix") { 
    return(qc);
  } else if(what2Return == "minfi") {
    return(Mset);
  }
}

ChristensenLabMinfiAdoption <- function(rgSet, detP, sampThresh) {
  #'@description Phase 2 of 3: Executes minfi preprocessing by the Funnorm approach & excludes
  #'@param rgSet RGChannelSet object of IDATs loaded by minfi
  #'@param detP Detection P-value threshold
  #'@param sampThresh Max proportion of samples allowed to fail detection P. Defaults to 20%
  
  ## Step 1. Executes Funnorm  & methylumi.noob background correction:
  genomRatSet <- preprocessFunnorm(rgSet);
  print(genomRatSet); 
  
  ## Step 2. Remove probes with low-call rate:
  print(paste(c("Detection P-value:","Sample threshold:"), c(detP, sampThresh))); 
  
  pvals <- detectionP(rgSet);
  failedP <- (pvals > detP); 
  
  print("Number of CpGs w/ failed P-values:");
  print(sum(failedP)); 
  
  print("Proportion of CpGs w/ failed P-values:");
  print(mean(failedP)); 
  
  print(paste("Number of CpGs w/ failed P-values in >", sampThresh, "of samples"));
  print(sum(rowMeans(failedP) > sampThresh)); 
  
  print(paste("Proportion of CpGs w/ failed P-values in >", sampThresh, "of samples"));
  print(mean(rowMeans(failedP) > sampThresh)); 
  
  failedProbes <- rownames(failedP)[rowMeans(failedP) > sampThresh];
  genomRatSet <- genomRatSet[! rownames(genomRatSet) %in% failedProbes];
  
  ## Step 3. Remove non-CpGs, control SNP probes, and polymorphic SNP probes:
  genomRatSet <- dropMethylationLoci(genomRatSet); #drop technical probes
  genomRatSet <- dropLociWithSnps(genomRatSet); #drop SNPs w/ default MAF=0
  
  return(genomRatSet);
}

extractMethMatrix <- function(genomRatSet, what=c("beta","M"), dropXY=FALSE, targets=NULL) {
  #'@description Phase 3 of 3: Extract beta-values from processed DNAm minfi object
  #'@param genomRatSet GRChannelSet minfi object
  #'@param what Methylation data type
  #'@param dropXY Should probes on chr X & Y be dropped?
  #'@param targets Optional data.frame for barcode-to-sample name conversion. A SAMPLE_NAME column needed
  
  if(what == "beta") {
    methMat <- minfi::getBeta(genomRatSet); 
  } else if(what == "M") {
    methMat <- minfi::getM(genomRatSet); 
  }
  
  if(dropXY) {
    print("Dropping chrX & Y probes...");
    sexProbes <- getXandorYProbes(chroms=c("chrX","chrY"));
    methMat <- methMat [! rownames(methMat) %in% sexProbes, ];
  }
  
  ## Map array barcode to sample ID:
  if(! is.null(targets) & ! any(is.na(targets$Sample_Name)) ) {
    if(identical(colnames(methMat), targets$Complete_Barcode) & ! all(duplicated(targets$Sample_Name)) ) {
      print("Switching barcode to sample ID based on minfi sample sheet...");
      colnames(methMat) <- targets$Sample_Name; #need to be unique
    } else {
      stop("Sample barcode must be matched before ID conversion!");
    }
  }
  
  return(methMat); 
}

getXandorYProbes <- function(chroms=c("chrX","chrY")) {
  #'@description Helper function to a character vector of CpGs on chr X and/or Y
  #'@param chroms Chromosomes to drop
  annot.850kb3 <- loadEPICannotationFile();
  xyProbes <- annot.850kb3$Name[annot.850kb3$chr %in% chroms]; 
  return(xyProbes);
}

## Execute pipelines as a single program **EDIT AS NEEDED**
## Alternatively, run commands in `Main` line by line
Main <- function() {
  #------------------------------------Step 1. Load IDAT files as a RGChannelSet------------------------------------
  setwd(IDAT_DIR);
  targets <- read.metharray.sheet(getwd());
  
  rgSet <- read.metharray.exp(targets=targets, extended=TRUE);
  print(rgSet@annotation);
  
  #------------------------------------Step 2. Execute custom sets of QCs (minfi & ENmix)------------------------------------
  qc <- ChristensenLabMethArrayQCs(rgSet, detP=DET_P, sampThresh=SAMP_FRACTION, dpi=DPI, outputPath=OUTPUT_PATH, what2Return="ENmix", usedFFPE=FALSE);
  print( paste("Total number of bad samples:", length(qc$badsample)) );
  save(list=c("rgSet","targets", "qc"), file="../../081518_NonCF_BAL_raw_rgSet_and_qc.RData", compress=TRUE); #saved for future queries
  
  ## Plot ENmix overall intensity:
  png("ENmix_overall_intensities.png", height=8.27, width=11.69, units="in", res=DPI);
  par(mar=c(10.5,4,4,2));
  barplot(qc$bisul, las=2, border=NA, ylab="Intensity", cex.names=0.75, main="Total Intensities by Sample, Before Preprocessing");
  dev.off();
  
  #------------------------------------Step 3. Basic raw-data exploration------------------------------------
  if(getwd() != OUTPUT_PATH) setwd(OUTPUT_PATH);
  
  ## Raw methylation density curves: 
  png("minfi_raw_densityCurves.png", height=8.27, width=11.69, units="in", res=DPI);
  par(mar=c(5,4,4,2));
  densityPlot(rgSet, main="Raw beta-value densities", sampGroups=targets$Bead_Chip);
  dev.off();
  
  png("minfi_raw_beanPlots.png", height=15, width=11.69, units="in", res=DPI);
  par(mar=c(5,10.5,4,2));
  densityBeanPlot(rgSet, main="Raw beta-value densities by sample", sampNames=targets$Sample_Name, sampGroups=targets$Bead_Chip);
  dev.off();
  
  #------------------------------------Step 4. Run the customized minfi pipeline------------------------------------
  normalizedEPIC <- ChristensenLabMinfiAdoption(rgSet, detP=DET_P, sampThresh=SAMP_FRACTION);
  allHealthyBetas <- extractMethMatrix(normalizedEPIC, what="beta", dropXY=TRUE, targets=targets);
  
  #------------------------------------Step 5. Plot final distributions & save------------------------------------
  save(list=c("allHealthyBetas","targets"), file="../../081518_NonCF_betas.RData", compress=TRUE);
  
  ## MDS plots:
  png("minfi_mdsPlotPanels_after_normalization.png", height=11.69, width=11.69, units="in", res=DPI);
  par(mfrow=c(2,2), mar=c(5,4,4,2));
  for(nCpGs in c(1e3, 1e4)) {
    mdsPlot(allHealthyBetas, numPositions=nCpGs, sampGroups=targets$Bead_Chip);
    mdsPlot(allHealthyBetas, numPositions=nCpGs, sampGroups=targets$LOBE);
  }
  dev.off();
  
  ## Methylation density plot:
  png("minfi_processed_density_curves_by_ENmix.png", height=8.27, width=11.69, units="in", res=DPI);
  par(mar=c(5,4,4,2));
  densityPlot(allHealthyBetas, main="Processed Beta-value Densities", sampGroups=targets$Bead_Chip);
  dev.off();
  
  png("minfi_processed_beanPlots_by_ENmix.png", height=15, width=11.69, units="in", res=DPI);
  par(mar=c(5,10.5,4,2));
  densityBeanPlot(allHealthyBetas, main="Processed Densities by Sample", sampGroups=targets$Bead_Chip);
  dev.off();
}

if(! interactive()) Main();
