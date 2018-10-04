# Visualize DMPs in RUL
# Author: David Chen
# Last update: 10/04/2018
# Notes:
# -- Run this script in RStudio

rm(list=ls());
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");
library(ggrepel);
library(limma);
library(matrixStats);

drawLabeledVolcanoPlot <- function(DMPs, pThresh, xThreshDown=-0.01, xThreshUp=0.5) {
  ## Direction of change:
  DMPs$dir[DMPs$adj.P.Val < pThresh & DMPs$logFC > 0] <- "Positive";
  DMPs$dir[DMPs$adj.P.Val < pThresh & DMPs$logFC < 0] <- "Negative";
  DMPs$dir[is.na(DMPs$dir)] <- "";
  
  nPos <- sum(DMPs$dir == "Positive", na.rm=TRUE);
  nNeg <- sum(DMPs$dir == "Negative", na.rm=TRUE);
  
  ## Gene names as labels:
  DMP.genes <- as.character(DMPs$UCSC_RefGene_Name); 
  DMP.genes <- strsplit(DMP.genes, split=";");
  DMPs$Gene <- rep(NA, length(DMP.genes));
  for(k in 1:length(DMP.genes)) {
    DMP.genes[[k]] <- unique(DMP.genes[[k]]); 
    if(length(DMP.genes[[k]]) > 1) {
      DMPs$Gene[k] <- paste(DMP.genes[[k]], collapse="; "); 
    } else if(length(DMP.genes[[k]]) == 1) {
      DMPs$Gene[k] <- DMP.genes[[k]]; 
    } else if(length(DMP.genes[[k]]) == 0) {
      DMPs$Gene[k] <- NA;
    }
  }
  DMPs$Label <- DMPs$Gene;
  DMPs$Label[(DMPs$logFC <= xThreshUp & DMPs$logFC >= xThreshDown) | DMPs$adj.P.Val >= pThresh] <- NA; #higher threshold for gene labeling
  
  ## Plot:
  plt <- ggplot(DMPs, aes(x=logFC, y=-log10(adj.P.Val), color=dir)) +
    geom_point(aes(size=dir, alpha=dir)) + #override
    scale_x_continuous(expand=c(0,0), limits=c(-1.75,1.75), breaks=seq(-1.5,1.5,by=0.5)) +
    scale_color_manual(values = c("gray","royalblue","olivedrab3")) +
    scale_size_manual(values=c(1, 1.1, 1.1)) +
    scale_alpha_manual(values=c(0.9, 1, 1)) +
    geom_text_repel(aes(label=Label), color="black", size=4) +
    geom_hline(yintercept=-log10(pThresh), linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed") +
    labs(x="log2FC(M-value)", y="-log10(FDR)") +
    annotate("text",1,3, label=paste(nPos,"HYPERmethyl. \n in RUL"),size=5,color="olivedrab3") +
    annotate("text",-1,3,label=paste(nNeg,"HYPOmethyl. \n in RUL"),size=5,color="royalblue") +
    myVolcanoTheme;
  
  print(plt);
}

Main <- function() {
  load(paste0(DATA_PATH, "081518_NonCF_betas.RData"));
  targets <- read.csv(paste0(DATA_PATH, "NonCF_updated_sample_sheet_090118.csv"), stringsAsFactors=FALSE);
  DMPs <- read.csv(paste0(DATA_PATH, "limma_DMP_100418/100418_NonCF_DMPs.csv"), stringsAsFactors=FALSE);
  
  png("~/Downloads/Diff_CpGs_in_RUL.png", res=300, units="in", height=8.27, width=11.69);
  drawLabeledVolcanoPlot(DMPs, FDR_THRESH, xThreshDown=-0.5, xThreshUp=0.45);
  dev.off();
}

Main();
