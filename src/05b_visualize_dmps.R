# Visualize DMPs in RUL
# Author: David Chen
# Last update: 10/04/2018
# Copyright (c) 2018-19 ydavidchen & Christensen-Lab
# Notes:
# -- Run this script in RStudio
# -- Recal that RUL has been, and will continue to be, the reference point

rm(list=ls());
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");
library(ggrepel);
library(limma);
library(matrixStats);

drawLabeledVolcanoPlot <- function(DMPs, pThresh=FDR_THRESH, xThreshDown=-0.5, xThreshUp=0.45) {
  #'@description Draws volcano plot for differentially methylated CpGs
  #'@param pThresh Threshold for adjusted P-values
  #'@param xThreshDown,xThreshUp Thresholds for coloring data points along the x-axis (logFC) dimension
  
  DMPs$negLog10Pval <- -log10(DMPs$adj.P.Val); 
  
  ## Direction of change:
  DMPs$dir[DMPs$adj.P.Val < pThresh & DMPs$logFC > 0] <- "Positive";
  DMPs$dir[DMPs$adj.P.Val < pThresh & DMPs$logFC < 0] <- "Negative";
  DMPs$dir[is.na(DMPs$dir)] <- "";
  
  ## Numbers to label:
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
  
  ## First pass for gene labels: systematically remove genes taht didn't meet threshold:
  DMPs$Label[(DMPs$logFC <= xThreshUp & DMPs$logFC >= xThreshDown) | DMPs$adj.P.Val >= pThresh] <- NA; #higher threshold for gene labeling
  
  ## Second pass for gene labels: manually add back or further remove labels as requested:
  boolAddBack <- DMPs$Gene %in% c("NR4A1","CLIP4","RUNX2","SNX10","SCNN1A") & DMPs$adj.P.Val < pThresh; 
  boolAddBack <- boolAddBack | DMPs$negLog10Pval >= 2; #a few more points
  DMPs$Label[boolAddBack] <- DMPs$Gene[boolAddBack];
  
  boolFurtherRemove <- DMPs$Gene %in% c("NFATC3","GFOD1","DENND1A");
  DMPs$Label[boolFurtherRemove] <- NA; 
  
  ## Third pass for gene labels: simplify/reformat view on plot:
  DMPs$Label <- gsub("; CMSS1; FILIP1L", "", DMPs$Label);
  DMPs$Label <- gsub("; RNF103-CHMP3", "", DMPs$Label);
  
  plt <- ggplot(DMPs, aes(x=logFC, y=negLog10Pval, color=dir)) +
    geom_point(aes(size=dir, alpha=dir)) + #override
    scale_x_continuous(expand=c(0,0), limits=c(-1.75,1.75), breaks=seq(-1.5,1.5,by=0.5)) +
    scale_y_continuous(limits=c(0,2.5)) +
    scale_color_manual(values = c("gray","royalblue","olivedrab3")) +
    scale_size_manual(values=c(1, 1.1, 1.1)) +
    scale_alpha_manual(values=c(0.9, 1, 1)) +
    geom_text_repel(aes(label=Label), color="black", size=4) +
    geom_hline(yintercept=-log10(pThresh), linetype="dashed") +
    geom_vline(xintercept=0, linetype="dashed") +
    labs(x="log2FC(M-value)", y="-log10(FDR)") +
    annotate("text",1.25,2, label=paste(nPos,"HYPERmethyl. \n in RUL"),size=7,color="olivedrab3") +
    annotate("text",-1.25,2,label=paste(nNeg,"HYPOmethyl. \n in RUL"),size=7,color="royalblue") +
    myVolcanoTheme;
  
  return(plt);
}

Main <- function() {
  ## Results & original data loading:
  DMPs <- read.csv(paste0(DATA_PATH, "limma_DMP_100418/100418_NonCF_DMPs.csv"), stringsAsFactors=FALSE);
  
  ## Volcano plot:
  png("~/Downloads/Figure3_Diff_CpGs_in_RUL.png", res=300, units="in", height=8.27, width=11.69);
  drawLabeledVolcanoPlot(DMPs, FDR_THRESH, xThreshDown=-0.5, xThreshUp=0.45);
  dev.off();
  
  ## Calculate mean delta beta-values & update DMP spreadsheet:
  samps_rul <- targets$Sample_Name[targets$LOBE == "RUL"];
  samps_rll <- targets$Sample_Name[targets$LOBE == "RLL"];
  stopifnot(all(colnames(allHealthyBetas) %in% c(samps_rul, samps_rll))); #checkpoint
  meanBetasByGroup <- calculateDeltaBetas(
    allHealthyBetas,
    group1=samps_rul, group2=samps_rll, 
    g1Name="RUL_mean_beta", g2Name="RLL_mean_beta"
  );
  dmps_with_delta_beta <- merge(DMPs, meanBetasByGroup, by="Name");
  write.csv(dmps_with_delta_beta, file="~/Downloads/100418_NonCF_DMPs_with_Mean_Delta_Betas.csv", row.names=FALSE, quote=FALSE);  
}

Main();
