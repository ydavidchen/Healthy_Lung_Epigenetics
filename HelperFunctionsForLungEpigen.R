# Helper Functions for Lung Epigenetics
# Script author: David Chen
# Script maintainer: David Chen
# Date: (ongoing)
# Copyright (c) 2018 ydavidchen & Christensen-Lab
# Notes:
# -- This script consists of 1) universal constants/paths, 2) methods, 3) ggplot2 themes

require(ggplot2);
require(doParallel); registerDoParallel(detectCores()-1);

TWO_GROUP_COLORS <- c("royalblue","olivedrab3");  

#------------------------------------- CONSTANTS & PATHS -------------------------------------
DATA_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/"; 
VAR_THRESH <- 0.01;
MVAl_THRESH <- 0.5; 
FDR_THRESH <- 0.05;

#------------------------------------- Data Loading Methods -------------------------------------
load_epigen_data <- function(path="~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/healthy_epigendx_validation/190125 EpigenDx Healthy Epic Data.csv",
                             totalReads=FALSE, dropExtraCols=TRUE) {
  dat <- data.table::fread(path, stringsAsFactors=FALSE, strip.white=TRUE, data.table=FALSE);
  if(totalReads) {
    dat <- subset(dat, Category == "Total Reads")
  } else {
    dat <- subset(dat, Category == "Percent Methylation");
  }
  dat$Sample_Name <- paste(gsub("-","",dat$Subject, fixed=TRUE), dat$Lobe, sep="_");
  if(dropExtraCols) dat$Barcode <- dat$`Customer ID` <- dat$`EpigenDx ID` <- dat$DNA_conc <- dat$Tube_label <- dat$Lobe <- dat$Subject <- dat$Category <- NULL;
  return(dat);
}

loadEPICannotationFile <- function() {
  #'@description Loas MethylationEPIC 850K annotation as a data.frame
  require(IlluminaHumanMethylationEPICanno.ilm10b3.hg19);
  annot.850kb3 <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b3.hg19));
  annot.850kb3$isEnhancer <- annot.850kb3$X450k_Enhancer=="TRUE" | annot.850kb3$Phantom4_Enhancers != "" | annot.850kb3$Phantom5_Enhancers != ""; 
  annot.850kb3$isPromoter <- grepl("TSS", annot.850kb3$UCSC_RefGene_Group); 
  return(annot.850kb3); 
}

createCpGTrackingBars <- function() {
  #'@description Creates pheatmap CpG annotation data.frame for tracking bars
  annot.850kb3 <- loadEPICannotationFile();
  row_annot <- data.frame(
    row.names = annot.850kb3$Name, 
    Island = ifelse(annot.850kb3$Relation_to_Island=="Island", "Yes", "No"),
    Promoter = ifelse(annot.850kb3$isPromoter, "Yes", "No"),
    Enhancer = ifelse(annot.850kb3$isEnhancer, "Yes", "No")
  );
  return(row_annot); 
}

#------------------------------------- Exploration & Statistical Methods -------------------------------------
decideNumberOfMostVariable <- function(data, varThresh, plot=TRUE) {
  #'@description Subset a matrix of CpGs based on sample variance
  #'@param data Matrix of CpGs with rows=CpGs, columns=samples
  #'@param varThresh Variance threshold for color in red
  #'@param plot Should a ranked variance distribution be shown?
  require(matrixStats);
  
  vars <- matrixStats::rowVars(data);
  vars <- sort(vars, decreasing=TRUE);
  bool <- (vars >= varThresh); 
  k <- sum(bool);
  
  if(plot) {
    plot(
      vars,
      col = ifelse(bool, "red", "black"),
      cex = 0.3,
      bty = "l",
      xlab = "CpGs",
      ylab = "Variance",
      main = "Inter-sample Variance Distribution"
    );
    abline(h=varThresh, lty=2);
    text(1.5e5, 0.10, paste(k, "CpGs"), col="red");
  }
  
  return(k); 
}

selectMostVariableCpGs <- function(data, k) {
  #'@description Subset a matrix of CpGs based on sample variance
  #'@param data Matrix of CpGs with rows=CpGs, columns=samples
  #'@param k Number of most variable CpGs to select
  require(matrixStats); 
  sele <- order(matrixStats::rowVars(data), decreasing=TRUE)[1:k];
  mat <- data[sele, ];
  return(mat); 
}

calculateDeltaBetas <- function(betas, group1, group2, g1Name=NULL, g2Name=NULL) {
  #'@description Calculates mean beta-values by group and then the differences
  #'@param betas Matrix of beta values. Rows=CpGs, Columns=Samples
  #'@param group1,group2 String vector of sample (column) names
  #'@param g1Name,g2Name Optional. Names for `group1`,`group2` columns
  
  mean_betas <- data.frame(
    Name = rownames(betas), 
    group1 = matrixStats::rowMeans2(betas[ , colnames(betas) %in% group1]),
    group2 = matrixStats::rowMeans2(betas[ , colnames(betas) %in% group2])
  );
  
  mean_betas$g1_minus_g2 <- mean_betas$group1 - mean_betas$group2;
  
  if(! is.null(g1Name) & ! is.null(g2Name)) {
    colnames(mean_betas)[colnames(mean_betas) == "group1"] <- g1Name;
    colnames(mean_betas)[colnames(mean_betas) == "group2"] <- g2Name;
    colnames(mean_betas)[colnames(mean_betas) == "g1_minus_g2"] <- paste(g1Name,"minus",g2Name,sep="_");
  }
  
  return(mean_betas);
}

#------------------------------------- Graphic Themes -------------------------------------
myScatterTheme <- theme_classic() + 
  theme(axis.text.x=element_text(size=20,color="black"), axis.title.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"), axis.title.y=element_text(size=20,color="black"),
        strip.text.x=element_text(size=20,colour="black",face="bold"),
        legend.position="top", legend.title=element_blank(), legend.text=element_text(size=15,color="black")); 

myBarplotTheme <- theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=20,color="black"),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20, color="black"),
        legend.position="top",legend.title=element_blank(),legend.text=element_text(size=12,color="black") ); 

myBoxplotTheme <- theme_classic() +
  theme(axis.text.x=element_text(size=21,color="black"), axis.text.y=element_text(size=21,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=22,color="black"),
        strip.text.x=element_text(size=20,colour="black",face="bold"),
        legend.position="top", legend.title=element_text(size=20), legend.text=element_text(size=20,color="black") );

myVolcanoTheme <- theme_classic() +
  theme(axis.text.x=element_text(color="black",size=16),axis.title.x=element_text(size=21, color="black"), 
        axis.text.y=element_text(color="black",size=16),axis.title.y=element_text(size=21, color="black"),
        legend.title=element_blank(), legend.text=element_blank(), legend.position="none"); 

