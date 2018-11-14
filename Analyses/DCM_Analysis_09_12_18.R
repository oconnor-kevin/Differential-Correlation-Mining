# TITLE: DCM_Analysis_09_12_18.R
# AUTHOR: Kevin O'Connor
# DATE MODIFIED: 9/12/18

# In this script, we run DCM on row-normalized RNAseq data.

# Switches.
is.test <- FALSE

# Libraries and directories.
library(R.utils)
library(dplyr)
library(limma)

sourceDirectory("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/DCM/R")
sourceDirectory("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/Utils")
if(is.test){
  data.dir <- file.path("/Users/kevinoconnor/Documents/Research/DCM/Test")  
} else {
  data.dir <- file.path("/Users/kevinoconnor/Documents/Research/DCM/Data")
}
out.dir <- filePath(data.dir, gsub("-", "_", Sys.time()))
dir.create(out.dir)
setwd(out.dir)

# Reading data.
if(!(exists("SeqData") & exists("sampleTab"))){
  SeqData   <- read.delim(file.path(data.dir, "BRCA_1201_FINAL_Cluster.txt") , 
                          header=F)
  sampleTab <- read.delim(file.path(data.dir, "BRCA_1218_pam50scores_FINAL.txt"))
}
data   <- SeqData[-c(1:5),-(1:3)] 
SamBar <- SeqData[,-c(1:3)]
samTab <- filter(sampleTab, Use=="YES")
gene.names <- SeqData[-c(1:5), 1]
rownames(data) <- gene.names

# Match barcodes in both datasets.
SamBar0 <- unlist(lapply(SamBar[1,], as.character))
m       <- match(SamBar0, as.character(samTab$Barcode))
samTab0 <- samTab[m, ]
dataM0  <- apply(data, 2, as.numeric)
rownames(dataM0) <- as.vector(gene.names)

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

# Standardize row within each group.
stdize <- function(gene){
  # Make zero mean.
  gene <- gene - mean(gene)
  
  # Make sum of squares one.
  gene <- gene/sqrt(sum(gene^2))
  
  return(gene)
}

MAM.LA.std <- MAM.LA
for(r in 1:nrow(MAM.LA)){
  MAM.LA.std[r,] <- stdize(MAM.LA[r,])
}

MAM.LB.std <- MAM.LB
for(r in 1:nrow(MAM.LB)){
  MAM.LB.std[r,] <- stdize(MAM.LB[r,])
}

MAM.basal.std <- MAM.basal
for(r in 1:nrow(MAM.basal)){
  MAM.basal.std[r,] <- stdize(MAM.basal[r,])
}

# Running DCM
lumA.lumB.50 <- DCM_KO(MAM.LA.std,
                       MAM.LB.std, 
                       max.iter   = 10, 
                       max.time   = 100, 
                       max.groups = 1,
                       alpha      = .05,
                       init.size  = 50)
save(lumA.lumB.50, file=filePath(out.dir, "DCM.lumA.vs.lumB.50.RData"))

# Load data.
load(file.path(out.dir, "debug_output/Initialization.RData"))
load(file.path(out.dir, "debug_output/DCM_1.RData"))

# Compute differential correlation.
A <- DCM$it_sets[[length(DCM$it_sets)-1]]
corA <- MAM.LA.std[A,] %*% t(MAM.LA.std[A,])
corB <- MAM.LB.std[A,] %*% t(MAM.LB.std[A,])
corA.tot <- MAM.LA.std %*% t(MAM.LA.std)
corB.tot <- MAM.LB.std %*% t(MAM.LB.std)

# Means
mean(corA-corB)
mean(corA.tot - corB.tot)

# Histograms
hist(corA-corB)
hist(corA.tot-corB.tot)

notA.sub <- sample(setdiff(1:ncol(corA.tot), A), 100)
A.sub    <- sample(A, 100)
corA.sub <- corA.tot[c(A.sub, notA.sub),c(A.sub, notA.sub)]
corB.sub <- corB.tot[c(A.sub, notA.sub),c(A.sub, notA.sub)]
pdf(file.path(out.dir, "Heatmap.pdf"))
heatmap((corA.sub-corB.sub)/(max(corA.sub-corB.sub)-min(corA.sub-corB.sub)), Rowv=NA, Colv=NA, scale="row")
dev.off()

mean(corA.tot[notA.sub, notA.sub]-corB.tot[notA.sub, notA.sub])
mean(corA.tot[A.sub, A.sub]-corB.tot[A.sub, A.sub])


# Check differential correlation throughout the procedure.
pdf(file.path(out.dir, "Mean_Differential_Correlation.pdf"))
plot(c(unlist(initialization$diff_cors), unlist(DCM$it_diff_cor)), main="Mean Differential Correlation by Iteration", ylab="Mean Differential Correlation")
abline(v=length(initialization$diff_cor))
dev.off()

# Plot set sizes.
pdf(file.path(out.dir, "Set_Size.pdf"))
plot(sapply(DCM$it_sets, function(s){length(s)}), main="Set sizes", ylab="Set size")
dev.off()

