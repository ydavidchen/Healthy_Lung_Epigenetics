# Exploratory Data Analysis (EDA) of Healthy BAL DNAm
# Author: David Chen
# Last update: 03/17/2019
# Notes:
# -- This file is dedicated for EDA and is frequently updated

rm(list=ls());
library(ggbiplot);
library(ggpubr);
library(gridExtra); 
library(matrixStats);
library(pheatmap);
library(reshape2);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsForLungEpigen.R");

load(paste0(DATA_PATH, "081518_NonCF_betas.RData"));
targets <- load_most_recent_covar(); #override original target

## Select most variable CpGs for further exploration:
k <- decideNumberOfMostVariable(allHealthyBetas, varThresh=VAR_THRESH, plot=TRUE);
k
mat <- selectMostVariableCpGs(allHealthyBetas, k=k);

## Ensure samples are matched:
if(! identical(colnames(mat), targets$Sample_Name)) {
  print("Proceed to matching samples...");
  mat <- mat[ , match(targets$Sample_Name, colnames(mat))]; 
} else {
  print("Already matched!");
}
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
                             groups=targets$Lobe) +
  ggtitle("(C) Lobe") +
  myScatterTheme; 

pcaPlotList[[1]]$labels$x <- pcaPlotList[[3]]$labels$x <- ""; 
pcaPlotList[[2]]$labels$y <- pcaPlotList[[3]]$labels$y <- "";
grid.arrange(grobs=pcaPlotList, nrow=1);

#-----------------------------------------Heat map-----------------------------------------
mat4Heat <- mat;
stopifnot(identical(targets$Sample_Name, colnames(mat4Heat))); #if fails, match first
colnames(mat4Heat) <- targets$Subject_ID;

## Sample annotation:
heatAnnot <- data.frame(
  row.names = targets$Subject_ID,
  Lobe = targets$Lobe,
  RPMM = targets$RPMM,
  stringsAsFactors = FALSE
);

## CpG annotation:
rowAnnot <- createCpGTrackingBars(); 

## Annotation colors:
annColors <- list(
  Lobe = c(Upper="black", Lower="lightgray"),
  RPMM = c(A1="black", A2="gray25", B1="gray55", B2="lightgray"),
  Promoter = c(Yes="black", No="lightgray"),
  Enhancer = c(Yes="black", No="lightgray"),
  Island = c(Yes="black", No="lightgray")
);

## Draw heat map with RPMM order:
rpmmOrderedSubjecIds <- targets$Subject_ID[order(targets$RPMMSampleOrder)];
mat4Heat <- mat4Heat[ , match(rpmmOrderedSubjecIds, colnames(mat4Heat))]; 

pheatmap(
  mat4Heat,
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_cols = FALSE,
  gaps_col = c(4,8,16),
  # gaps_col = c(4,8,10,12,16,18,22,26), #for max=Inf
  annotation_col = heatAnnot,
  annotation_row = rowAnnot,
  annotation_colors = annColors,
  color = colorRampPalette(c("yellow", "black", "blue"))(512),
  border_color = NA, 
  fontsize = 12,
  fontsize_col = 8
)

## Examine association between RPMM cluster A vs. Upper Lobe:
contTab <- table(
  RPMM = ifelse(targets$RPMM %in% c("A1","A2"), "As", "Bs"),
  Lobe = targets$Lobe
);
contTab <- contTab[ , c(2,1)]; 
contTab

#-----------------------------------------Cell Type Comparison-----------------------------------------
plt_cellType <- melt(targets[ , c("Sample_Name","CellType_1","CellType_2","Lobe")]); 
plt_cellType$variable <- gsub("CellType_1", "Cell Type 1", plt_cellType$variable);
plt_cellType$variable <- gsub("CellType_2", "Cell Type 2", plt_cellType$variable);
plt_cellType$Lobe <- factor(plt_cellType$Lobe, levels=c("Upper","Lower"));

ggplot(plt_cellType, aes(x=variable, y=value, color=Lobe)) +
  geom_boxplot(outlier.colour=NA, outlier.fill=NA, outlier.shape=NA) +
  geom_point(position=position_jitterdodge(jitter.width=0.25)) +
  scale_color_manual(values=c("blue","red")) +
  scale_y_continuous(limits=c(0, 1.05)) +
  labs(y="Proportion") +
  myBoxplotTheme;

#-----------------------------------------Table One for Patients-----------------------------------------
generate_table_one <- function() {
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
  write.csv(y, file="~/Downloads/083118_NonCF_table1.csv",quote=F)
}

# generate_table_one();

#-----------------------------------------Methylation age comparison-----------------------------------------
summary(lm(HorvathAge ~ AGE, data=targets));
cor(targets$HorvathAge, targets$AGE, method="pearson");

label <- expression(atop(
  'Pearson Corr Coef. ='~0.66~',',
  'Adj. R'^2~'='~0.41~', P ='~2.69E-04
));

ggscatter(
  targets, x="HorvathAge", y="AGE",
  xlab="Horvath Methyl. Age", ylab = "Actual Age", 
  size = 3, add = "reg.line", conf.int=TRUE #, cor.coef=TRUE, cor.coeff.args=list(method="pearson", label.x=38, label.sep=", ")  , 
) + myScatterTheme +
  annotate("text", 35, 35, label=label, size=7)

