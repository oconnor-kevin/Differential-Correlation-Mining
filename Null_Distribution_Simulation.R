# TITLE: Null_Distribution_Simulation.R
# AUTHOR: Kevin O'Connor
# DATE MODIFIED: 6/20/18

# Switches.
is.test <- FALSE

# Libraries and directories.
library(R.utils)
library(dplyr)
library(limma)

sourceDirectory("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/DCM/R")
if(is.test){
  setwd("/Users/kevinoconnor/Documents/Research/DCM/Null Distribution Simulation/Test")  
} else {
  setwd("/Users/kevinoconnor/Documents/Research/DCM/Null Distribution Simulation")
}
data.dir <- file.path("/Users/kevinoconnor/Documents/Research/DCM")


# Reading data.
if(!(exists("SeqData") & exists("sampleTab"))){
  SeqData   <- read.delim(file.path(data.dir, "BRCA.1201.FINAL.Cluster.txt") , header=F)
  sampleTab <- read.delim(file.path(data.dir, "BRCA.1218_pam50scores.FINAL.txt"))
}

data   <- SeqData[-c(1:5),-(1:3)] 
SamBar <- SeqData[,-c(1:3)]

samTab <- filter(sampleTab, Use=="YES")

## Matching barcodes in both datasets.
SamBar0 <- unlist(lapply(SamBar[1,], as.character))
m <- match(SamBar0, as.character(samTab$Barcode))
##m[1:5]
##table(is.na(m))

samTab0 <- samTab[m, ]
# matched with data

dataM=apply(data,2, as.numeric)
design <- model.matrix(~0+samTab0$Call)
##dim(dataM)
## [1] 20531  1201
##hist(dataM[,1] ) # How the reading data was generated? RPKM? There are zeros
##hist(log2(dataM[,1] ))

## Replacing zeroes with minimum value.
minD <- min(apply(dataM, 2,function(x) min(x[x!=0])  ))
###[1] 0.0021
dataM0 <- dataM
dataM[dataM==0] <- minD

## Taking log of data.
logD <- log2(dataM)

## Renaming columns.
colnames(design) <- unlist(lapply(strsplit(colnames(design), "Call"), function(x) x[[2]]))
###colnames(design)
#### "Basal"  "Her2"   "LumA"   "LumB"   "Normal"

i <- apply(dataM==0, 1, all)  #all false
table(i)
# min sample size 82
numZeroPerG <- apply(dataM0, 1, function(x) length(which(x==0)) )
hist(numZeroPerG)

## Create matrix for each group.
tar <- samTab0$Call
wb  <- which(tar=="Basal")
wla <- which(tar=="LumA")
wlb <- which(tar=="LumB")
MAM.basal <- logD[, wb]
MAM.LA    <- logD[, wla]
MAM.LB    <- logD [, wlb] 

###dim(MAM.basal) # 20531 191
###dim(MAM.LA) # 572
###dim(MAM.LB) #219


# Running simulation.
## Create output directory.
out.dir <- filePath(getwd(), gsub("-", "_", Sys.time()))
dir.create(out.dir)
setwd(out.dir)

start.time <- Sys.time()
n.iter <- 10
p.vals.l <- list()
for(i in 1:n.iter){
  ## Take random sample of data under null.
  rand.inds <- sample(1:ncol(MAM.LA), round(ncol(MAM.LA)/2))
  MAM.LA.1 <- MAM.LA[,rand.inds]
  MAM.LA.2 <- MAM.LA[,-rand.inds]
  
  ## Run DCM on null data.
  DCM.Null.LA <- DCM_Kevin(MAM.LA.1, 
                           MAM.LA.2, 
                           max.iter = 2,
                           max.groups = 1,
                           max.time = 100, 
                           alpha = .05,  
                           strict='low',
                           echo=TRUE, 
                           QN=TRUE, 
                           resid.full=TRUE)
  load(file.path(out.dir, "Debug_Output/DCM_1.RData"))
  p.vals <- DCM$"pvals"
  p.vals.l[[i]] <- p.vals
}
save(p.vals.l, file=file.path(out.dir, "Null_P_Vals.RData"))
difftime(Sys.time(), start.time)

# Making histograms of p-values.
load(file.path(out.dir, "Null_P_Vals.RData"))
pdf(file=file.path(out.dir, "Null_P_Vals.pdf"))
for(i in 1:n.iter){
  hist(p.vals.l[[i]], main=paste0("P-Values from Random Sample ", i), xlab="P-Values")
}
dev.off()
