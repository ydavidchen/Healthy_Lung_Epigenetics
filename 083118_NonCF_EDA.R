# Exploratory Data Analysis of Healthy BAL DNAm
# Author: David Chen
# Date: 08/31/2018
# Notes:

rm(list=ls());
library(ggbiplot);
library(ggpubr);
library(gridExtra); 
library(matrixStats);
library(pheatmap); 
library(tableone);
library(doParallel); registerDoParallel(detectCores()-1);

setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");

DATA_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2018/NonCF_Healthy_EPIC/"; 
VAR_THRESH <- 0.01; 
load(paste0(DATA_PATH, "081518_NonCF_betas.RData"));
targets <- read.csv(paste0(DATA_PATH, "NonCF_updated_sample_sheet_082018.csv"), stringsAsFactors=FALSE);

## Select most variable CpGs for further exploration:
k <- decideNumberOfMostVariable(allHealthyBetas, varThresh=VAR_THRESH, plot=TRUE);
k
mat <- selectMostVariableCpGs(allHealthyBetas, k=k);
stopifnot(identical(colnames(mat), targets$Sample_Name)); #checkpoint

#-----------------------------------------PCA-----------------------------------------
topKsites.pca <- prcomp(t(mat));

pcaPlotList <- list();
pcaPlotList[[1]] <- ggbiplot(topKsites.pca, var.scale=1, varname.size=0, obs.scale=1, var.axes=FALSE, ellipse=FALSE, circle=FALSE, 
                             groups=as.factor(targets$Bead_Chip)) +
  ggtitle("(A) Chip ID (batch)") +
  myScatterTheme; 

pcaPlotList[[2]] <- ggbiplot(topKsites.pca, var.scale=1, varname.size=0, obs.scale=1, var.axes=FALSE, ellipse=FALSE, circle=FALSE, 
                             groups=targets$GENDER) +
  ggtitle("(B) Gender") +
  myScatterTheme;

pcaPlotList[[3]] <- ggbiplot(topKsites.pca, var.scale=1, varname.size=0, obs.scale=1, var.axes=FALSE, ellipse=FALSE, circle=FALSE, 
                             groups=targets$LOBE) +
  ggtitle("(C) Lobe") +
  myScatterTheme; 

pcaPlotList[[1]]$labels$x <- pcaPlotList[[3]]$labels$x <- ""; 
pcaPlotList[[2]]$labels$y <- pcaPlotList[[3]]$labels$y <- "";
grid.arrange(grobs=pcaPlotList, nrow=1);

#-----------------------------------------Heat map-----------------------------------------
targets <- targets[order(targets$RPMMSampleOrder), ];

## Sample annotation:
heatAnnot <- data.frame(
  row.names = targets$Sample_Name,
  Lobe = targets$LOBE,
  RPMM = targets$RPMMClusters,
  stringsAsFactors = FALSE
);

## CpG annotation:
rowAnnot <- createCpGTrackingBars(); 

## Annotation colors:
annColors <- list(
  Lobe = c(`RUL`="black", `RLL`="lightgray"),
  RPMM = c(rLL="black", rLR="gray25", rRL="gray55", rRR="lightgray"),
  Promoter = c(Yes="black", No="lightgray"),
  Enhancer = c(Yes="black", No="lightgray"),
  Island = c(Yes="black", No="lightgray")
);

## Select most variable CpGs & draw heat map with RPMM order:
sampOrders <- targets$Sample_Name[order(targets$RPMMSampleOrder)];
pheatmap(
  mat[ , match(sampOrders, colnames(mat))],
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_cols = FALSE,
  gaps_col = c(4,8,16),
  # gaps_col = c(4,8,10,12,16,18,22,26), # max=Inf
  annotation_col = heatAnnot,
  annotation_row = rowAnnot,
  annotation_colors = annColors,
  color = colorRampPalette(c("yellow", "black", "blue"))(512),
  border_color = NA, 
  fontsize = 12,
  fontsize_col = 8
)


#-----------------------------------------Table One for Patients-----------------------------------------
targets_rul <- subset(targets, LOBE=="RUL");
colnames(targets_rul)[colnames(targets_rul)=="HorvathAge"] <- "Horvath.RUL";

targets_rll <- subset(targets, LOBE=="RLL");
colnames(targets_rll)[colnames(targets_rll)=="HorvathAge"] <- "Horvath.RLL";

df4Table1 <- merge(
  targets_rll[ , c("Subject","AGE","GENDER","STATUS","RPMMClusters","Horvath.RLL")],
  targets_rul[ , c("Subject","Horvath.RUL")],
  by = "Subject"
)

myFactorVars <- c("GENDER","STATUS","RPMMClusters");
myVars <- c(myFactorVars,"AGE","Horvath.RLL", "Horvath.RUL");
myTableOne <- CreateTableOne(
  data = df4Table1,
  vars = myVars,
  factorVars = myFactorVars,
  test = TRUE,
  includeNA = TRUE
);

y <- print(myTableOne, showAllLevels=TRUE); 
# write.csv(y, file="~/Downloads/083118_NonCF_table1.csv",quote=F)

#-----------------------------------------Methylation age comparison-----------------------------------------
summary(lm(HorvathAge ~ AGE, data=targets));
cor(targets$HorvathAge, targets$AGE, method="pearson");

label <- expression(atop(
  'Pearson Corr Coef. ='~0.66~',',
  'Adj. R'^2~'='~0.41~', P ='~2.7E-4
));

ggscatter(
  targets, x="HorvathAge", y="AGE",
  xlab="Horvath Methyl. Age", ylab = "Actual Age", 
  size = 3, add = "reg.line", conf.int=TRUE #, cor.coef=TRUE, cor.coeff.args=list(method="pearson", label.x=38, label.sep=", ")  , 
) + myScatterTheme +
  annotate("text", 35, 35, label=label, size=7)

