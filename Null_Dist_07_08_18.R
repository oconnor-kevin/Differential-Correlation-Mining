# TITLE: Null_Dist_07_08_18.R
# AUTHOR: Kevin O'Connor
# DATE MODIFIED: 7/8/18

# Switches and settings.
is.test <- TRUE
init.size.v <- c(50, 1500)
n.iter <- 10

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


# Running simulation with RNAseq data.
## Create output directory.
out.dir <- filePath(getwd(), gsub("-", "_", Sys.time()))
dir.create(out.dir)
setwd(out.dir)
for (k in init.size.v){
  ## Make directory for this seed size.
  size.dir <- filePath(out.dir, paste0("Seed Size ", k))
  dir.create(size.dir)
  setwd(size.dir)
  
  start.time <- Sys.time()
  p.vals.l <- list()
  for(i in 1:n.iter){
    ## Take random sample of data under null.
    rand.inds <- sample(1:ncol(MAM.LA), round(ncol(MAM.LA)/2))
    MAM.LA.1 <- MAM.LA[,rand.inds]
    MAM.LA.2 <- MAM.LA[,-rand.inds]
    
    ## Run DCM on null data.
    DCM.Null.LA <- DCM_Kevin(MAM.LA.1, 
                             MAM.LA.2, 
                             max.iter   = 1,
                             max.groups = 1,
                             max.time   = 100,
                             est.size   = k,
                             alpha      = .05,  
                             strict     = 'low',
                             echo       = TRUE, 
                             QN         = TRUE,
                             initialize = FALSE,
                             resid.full = TRUE)
    load(file.path(size.dir, "Debug_Output/DCM_1.RData"))
    p.vals <- DCM$"pvals"
    p.vals.l[[i]] <- p.vals
  }
  save(p.vals.l, file=file.path(size.dir, "Null_P_Vals.RData"))
  difftime(Sys.time(), start.time)
  
  # Making histograms of p-values.
  load(file.path(size.dir, "Null_P_Vals.RData"))
  pdf(file=file.path(size.dir, "Null_P_Vals.pdf"))
  for(i in 1:n.iter){
    hist(p.vals.l[[i]], main=paste0("P-Values from Random Sample ", i), xlab="P-Values")
  }
  dev.off()
}


# Running simulation with simulated data.
## Create output directories.
## Looking at LumA data before and after transformation.
la.vec <- as.vector(MAM.LA)
par(mfrow=c(2,1))
hist(la.vec)
hist(2^la.vec)
### Will remove extremely large values and replot.
hist(la.vec)
hist((2^la.vec)[(2^la.vec) < 1000])
### Looking at variances of each gene.
vars <- apply(MAM.LA, 1, var)
par(mfrow=c(1,1))
hist(vars)
hist(vars[vars > 0])
hist(vars[vars<=0.1])
sum(vars==0)
### Remove genes with 0 variance.
MAM.LA.old <- MAM.LA
MAM.LA <- MAM.LA[vars>0,]
### Showing frequency of minimum value in each gene.
min.val <- min(la.vec)
counts.of.min <- apply(MAM.LA, 1, function(i){sum(i==min.val)})
sum(counts.of.min==(ncol(MAM.LA)))

# Simulating data.
n.a <- dim(MAM.LA)[2]
n.b <- dim(MAM.LB)[2]
p <- dim(MAM.LA)[1]
sim.data.a <- rnorm(n.a*p) %>% matrix(nrow=p, ncol=n.a)
sim.data.b <- rnorm(n.b*p) %>% matrix(nrow=p, ncol=n.b)

## Create output directory.
out.dir <- filePath(getwd(), gsub("-", "_", Sys.time()))
dir.create(out.dir)
setwd(out.dir)

p.vals.l <- list()
n.iter <- 10
for(i in 1:n.iter){
  sim.data.a <- rnorm(n.a*p) %>% matrix(nrow=p, ncol=n.a)
  sim.data.b <- rnorm(n.b*p) %>% matrix(nrow=p, ncol=n.b)
  
  DCM.Null.SIM <- DCM_Kevin(sim.data.a,
                            sim.data.b, 
                            max.iter   = 10,
                            max.groups = 1,
                            max.time   = 100,
                            est.size   = 50,
                            alpha      = .05,  
                            strict     = 'low',
                            echo       = TRUE, 
                            QN         = FALSE,
                            initialize = FALSE,
                            resid.full = TRUE)
  load(file.path(out.dir, "Debug_Output/DCM_1.RData"))
  p.vals.l[[i]] <- DCM$"pvals"
}
save(p.vals.l, file=file.path(out.dir, "Null_P_Vals.RData"))
load(file.path(out.dir, "Null_P_Vals.RData"))
pdf(file=file.path(out.dir, "Null_P_Vals.pdf"))
for(i in 1:n.iter){
  hist(p.vals.l[[i]], main=paste0("P-Values from Random Sample ", i), xlab="P-Values")
}
dev.off()