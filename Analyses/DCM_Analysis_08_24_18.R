# TITLE: DCM_Analysis_08_24_18.R
# AUTHOR: Kevin O'Connor
# DATE MODIFIED: 8/24/18

# In this script, we are testing a new terminating condition in the DCM_Kevin.R
#  script. This condition halts computation once a search does not converge.

# Switches.
is.test <- TRUE

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
out.dir.1 <- file.path(out.dir, "No QR")
dir.create(out.dir.1)
setwd(out.dir.1)
lumA.lumB.50.noQR <- DCM_Kevin(MAM.LA.std,
                               MAM.LB.std, 
                               max.iter   = 10, 
                               max.time   = 100, 
                               max.groups = 3,
                               alpha      = .05,
                               est.size   = 50,
                               strict     = 'low', 
                               echo       = TRUE)
save(lumA.lumB.50.noQR, file=filePath(out.dir.1, "DCM.lumA.vs.lumB.noQR.50.RData"))

out.dir.2 <- file.path(out.dir, "QR")
dir.create(out.dir.2)
setwd(out.dir.2)
lumA.lumB.50.QR <- DCM_Kevin(MAM.LA.std,
                             MAM.LB.std, 
                             max.iter   = 10, 
                             max.time   = 100, 
                             max.groups = 3,
                             alpha      = .05,
                             est.size   = 50,
                             QN         = TRUE,
                             strict     = 'low', 
                             echo       = TRUE)
save(lumA.lumB.50.QR, file=filePath(out.dir.2, "DCM.lumA.vs.lumB.QR.50.RData"))


# Looking at overlaps between found sets and known gene sets.
gene.names.filtered <- rownames(MAM.LA)
found.genes.1 <- gene.names.filtered[lumA.lumB.50.noQR$DC_sets[[1]]]
found.genes.1 <- sapply(found.genes.1, function(gene.str){
  char.str <- strsplit(gene.str, '')[[1]]
  ind <- which(char.str=='|')[1]
  return(paste0(char.str[1:(ind-1)], collapse=""))
})

overlaps.1 <- sapply(Hs.gmtl.c2, function(gene.set){
    return(length(intersect(gene.set, found.genes.1)))
}) %>% sort(decreasing=TRUE)

gene.names.filtered <- rownames(MAM.LA)
found.genes.2 <- gene.names.filtered[lumA.lumB.50.noQR$DC_sets[[2]]]
found.genes.2 <- sapply(found.genes.2, function(gene.str){
  char.str <- strsplit(gene.str, '')[[1]]
  ind <- which(char.str=='|')[1]
  return(paste0(char.str[1:(ind-1)], collapse=""))
})

overlaps.2 <- sapply(Hs.gmtl.c2, function(gene.set){
  return(length(intersect(gene.set, found.genes.2)))
}) %>% sort(decreasing=TRUE)

genes.set.sizes <- sapply(Hs.gmtl.c2, length)
genes.set.sizes[names(overlaps.2)[1:10]]