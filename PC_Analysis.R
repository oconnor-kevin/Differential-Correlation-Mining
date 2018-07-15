# TITLE: PC_Analysis.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 7/10/18
# DATE MODIFIED: 7/10/18

# Switches.
is.test <- TRUE

# Libraries and directories.
library(R.utils)
library(dplyr)
library(limma)

sourceDirectory("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/DCM/R")
sourceDirectory("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/Utils",
                recursive=FALSE)
if(is.test){
  setwd("/Users/kevinoconnor/Documents/Research/DCM/Test")  
} else {
  setwd("/Users/kevinoconnor/Documents/Research/DCM")
}
data.dir <- file.path("/Users/kevinoconnor/Documents/Research/DCM/Data")

## Create output directory.
out.dir <- filePath(data.dir, gsub("-", "_", Sys.time()))
dir.create(out.dir)

# Reading data.
if(!(exists("SeqData") & exists("sampleTab"))){
  SeqData   <- read.delim(file.path(data.dir, "BRCA_1201_FINAL_Cluster.txt") , 
                          header=F)
  sampleTab <- read.delim(file.path(data.dir, "BRCA_1218_pam50scores_FINAL.txt"))
}
data   <- SeqData[-c(1:5),-(1:3)] 
SamBar <- SeqData[,-c(1:3)]
samTab <- filter(sampleTab, Use=="YES")

# Match barcodes in both datasets.
SamBar0 <- unlist(lapply(SamBar[1,], as.character))
m       <- match(SamBar0, as.character(samTab$Barcode))
samTab0 <- samTab[m, ]
dataM0  <- apply(data, 2, as.numeric)

# Filter and transform data.
dataM <- filter_data(dataM0,
                     min.var=100,
                     max.var=10000)

# Create matrix for each group.
tar <- samTab0$Call
wb  <- which(tar=="Basal")
wla <- which(tar=="LumA")
wlb <- which(tar=="LumB")
MAM.basal <- dataM[, wb]
MAM.LA    <- dataM[, wla]
MAM.LB    <- dataM[, wlb]


# Principal Components Analysis
# Compute PC's.
pca.la <- prcomp(t(MAM.LA), center=TRUE, scale=TRUE)
pca.lb <- prcomp(t(MAM.LB), center=TRUE, scale=TRUE)
pca.basal <- prcomp(t(MAM.basal), center=TRUE, scale=TRUE)
pca.la.lb <- prcomp(cbind(MAM.LA, MAM.LB), center=TRUE, scale=TRUE)
# Plot sorted absolute loadings for first three PC's.
## Lum A vs. Lum B observations
par(mfrow=c(3,1))
plot(abs(pca.la.lb$rotation[,1]))
abline(v=ncol(MAM.LA))
plot(abs(pca.la.lb$rotation[,2]))
abline(v=ncol(MAM.LA))
plot(abs(pca.la.lb$rotation[,3]))
abline(v=ncol(MAM.LA))
pdf(file.path(out.dir, "LumA_vs_LumB_Loadings.pdf"))
par(mfrow=c(3,1))
plot(pca.la.lb$rotation[,1])
abline(v=ncol(MAM.LA))
plot(pca.la.lb$rotation[,2])
abline(v=ncol(MAM.LA))
plot(pca.la.lb$rotation[,3])
abline(v=ncol(MAM.LA))
dev.off()
## Compare mean loading between observations in Lum A vs Lum B.
mean(abs(pca.la.lb$rotation[,1][1:ncol(MAM.LA)]))
mean(abs(pca.la.lb$rotation[,1][(ncol(MAM.LA)+1):ncol(MAM.LB)]))
mean(pca.la.lb$rotation[,1][1:ncol(MAM.LA)])
mean(pca.la.lb$rotation[,1][(ncol(MAM.LA)+1):ncol(MAM.LB)])

