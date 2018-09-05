# Helper Functions for Lung Epigenetics
# Author: David Chen
# Date: 08/15/2018
# Notes:

#------------------------------------- Data Loading Procedures -------------------------------------
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

#------------------------------------- Graphic Themes -------------------------------------
require(ggplot2);
myScatterTheme <- theme_classic() + 
  theme(axis.text.x=element_text(size=20,color="black"), axis.title.x=element_text(size=20,color="black"),
        axis.text.y=element_text(size=20,color="black"), axis.title.y=element_text(size=20,color="black"),
        strip.text.x=element_text(size=12,colour="black",face="bold"),
        legend.position="top", legend.title=element_blank(), legend.text=element_text(size=15,color="black")); 

myBarplotTheme <- theme_classic() +
  theme(axis.text.x=element_blank(), axis.text.y=element_text(size=20,color="black"),
        axis.ticks.x=element_blank(),
        axis.title=element_text(size=20, color="black"),
        legend.position="top",legend.title=element_blank(),legend.text=element_text(size=12,color="black") ); 

myBoxplotTheme <- theme_classic() +
  theme(axis.text.x=element_text(size=20,color="black"), axis.text.y=element_text(size=21,color="black"),
        axis.title.x=element_blank(), axis.title.y=element_text(size=20,color="black"),
        legend.position="top", legend.title=element_text(size=20), legend.text=element_text(size=15,color="black") );

