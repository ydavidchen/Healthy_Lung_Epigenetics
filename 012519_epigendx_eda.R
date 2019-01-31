# EpigenDX validation: EDA
# Script author: David Chen
# Script maintainer: David Chen
# Last update: 01/25/2019
# Copyright (c) 2018 ydavidchen & Christensen-Lab
# Notes:

rm(list=ls());
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");
library(matrixStats);
library(reshape2);

## CpG annotation loading:
epigen_cg_annot <- read.csv("~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/healthy_epigendx_validation/CpG_annotations.csv", stringsAsFactors=FALSE);
epigen_cg_annot <- subset(epigen_cg_annot, ! is.na(EPIC_ID));

## Data loading:
epigen_val <- load_epigen_data(totalReads=FALSE);
epigen_qual <- load_epigen_data(totalReads=TRUE);

load(paste0(DATA_PATH, "081518_NonCF_betas.RData"));
targets <- read.csv(paste0(DATA_PATH, "NonCF_updated_sample_sheet_090118.csv"), stringsAsFactors=FALSE);
DMPs <- read.csv(paste0(DATA_PATH, "limma_DMP_100418/100418_NonCF_DMPs.csv"), stringsAsFactors=FALSE);

epigen_val <- merge(epigen_val, targets, by="Sample_Name");
epigen_qual <- merge(epigen_qual, targets, by="Sample_Name");

#---------------------------------- Data visualization ----------------------------------
## Total reads check:
plt_qual <- epigen_qual[ , grepl("LOBE", colnames(epigen_qual)) | grepl("CpG", colnames(epigen_qual))]
plt_qual <- melt(plt_qual, variable.name="CpG", value.name="totalReads");
ggplot(plt_qual, aes(LOBE, totalReads, color=LOBE)) +
  geom_boxplot(outlier.colour=NA, outlier.shape=NA) +
  geom_jitter(width=0.25) +
  scale_color_manual(values=TWO_GROUP_COLORS) +
  ylab("Total Reads") + 
  facet_wrap(~ CpG, scales="free_y", nrow=2) +
  myBoxplotTheme;

## % Methylation:
plt_meth <- epigen_val[ , grepl("LOBE", colnames(epigen_val)) | grepl("CpG", colnames(epigen_val))];
plt_meth <- melt(plt_meth, variable.name="CpG", value.name="Methylation");
ggplot(plt_meth, aes(LOBE, Methylation, color=LOBE)) +
  geom_boxplot(outlier.colour=NA, outlier.shape=NA) +
  geom_jitter(width=0.25) +
  scale_color_manual(values=TWO_GROUP_COLORS) +
  facet_wrap(~CpG, scales="free_y", nrow=2) +
  ylab("% Methylation") +
  myBoxplotTheme;

#-------------------------------------- Compute delta betas by group --------------------------------------
rul_val <- rll_val <- epigen_val[ , grepl("CpG", colnames(epigen_val)) | grepl("Sample_Name", colnames(epigen_val)) | grepl("Subject", colnames(epigen_val))];

rul_val <- subset(rul_val, grepl("RUL", Sample_Name));
rll_val <- subset(rll_val, grepl("RLL", Sample_Name)); 

## Important: Check if CpGs & subjects match between matrices. 
## If not, manually match before proceeding
stopifnot(identical(rul_val$Subject, rll_val$Subject)); 
stopifnot(identical(colnames(rul_val), colnames(rll_val)));

## Perform element-wise subtractions:
rownames(rul_val) <- rul_val$Subject;
rownames(rll_val) <- rll_val$Subject;

rul_val$Subject <- rll_val$Subject <- rul_val$Sample_Name <- rll_val$Sample_Name <- NULL;

rul_val <- as.matrix(rul_val);
rll_val <- as.matrix(rll_val);

epigen_db <- rul_val - rll_val;
epigen_db <- rbind(epigen_db, Median=colMedians(epigen_db), Mean=colMeans2(epigen_db), SD=colSds(epigen_db));

# write.csv(epigen_db, file="~/Downloads/013119_Delta_Betas_RUL_minus_RLL.csv", row.names=TRUE, quote=FALSE);

#--------------------------------------  Comparison with EPIC array --------------------------------------
# SITES <- c("cg26118047","cg17323256","cg06086177","cg16711835","cg13427753","cg00889217");
# sele_betas_epic <- allHealthyBetas[rownames(allHealthyBetas) %in% SITES, ];

sele_betas_epic <- allHealthyBetas[rownames(allHealthyBetas) %in% epigen_cg_annot$EPIC_ID, ];
sele_betas_epic <- as.data.frame(t(sele_betas_epic));
sele_betas_epic$Sample_Name <- rownames(sele_betas_epic);
sele_betas_epic <- merge(sele_betas_epic, epigen_val, by="Sample_Name");
sele_betas_epic <- sele_betas_epic[ , grepl("cg",colnames(sele_betas_epic)) | grepl("CpG",colnames(sele_betas_epic))];

plt_scatter_list <- list();
for(i in 1:nrow(epigen_cg_annot)) {
  plt_scatter_list[[i]] <- ggplot(sele_betas_epic, aes_string(epigen_cg_annot$EPIC_ID[i], epigen_cg_annot$EpigenDx_ID[i])) +
    geom_point() +
    geom_smooth(method="lm") +
    myScatterTheme  
}
gridExtra::grid.arrange(grobs=plt_scatter_list, ncol=2);
