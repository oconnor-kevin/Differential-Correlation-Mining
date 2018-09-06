# TITLE: DCM_Analysis_08_28_18.R
# AUTHOR: Kevin O'Connor
# DATE MODIFIED: 8/28/18

# Filtering based on mean vs standard deviation.

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

# Filter and transform the data.
dataM <- filter_data(dataM0, min.var=10, max.var=1e24)

# Investigate mean vs standard deviation.
mean.v <- apply(dataM, 1, mean)
mad.v  <- apply(dataM, 1, mad)
sd.v   <- apply(dataM, 1, sd)

plot(sd.v, mean.v)
hist(sd.v/mad.v)

# Find genes with mad/sd exceeding threshold.
th <- 2
bad.genes <- sd.v/mad.v > th
plot(sd.v, 
     mean.v, 
     col=ifelse(bad.genes, "red", "black"), 
     cex=0.3, 
     main="RNAseq Mean vs. SD")

#mat.1 <- sapply(1:ncol(dataM), function(c){c/mad.v[c]})
#scaled.data <- dataM / mad.v 
#bad.genes <- rowSums(scaled.data > 5) > 0
#plot(mean.v, sd.v, col=ifelse(bad.genes==1, "red", "black"), cex=0.3)
#mean(bad.genes)
#hist(dataM[!bad.genes, ])

# Filter out bad genes.
dataM.good <- dataM[!bad.genes,]

# Create matrix for each group.
tar <- samTab0$Call
wb  <- which(tar=="Basal")
wla <- which(tar=="LumA")
wlb <- which(tar=="LumB")
MAM.basal <- dataM.good[, wb]
MAM.LA    <- dataM.good[, wla]
MAM.LB    <- dataM.good[, wlb]

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

# Read microarray data.
load(file=file.path(data.dir, "microarray_filtered.RData"))
load(file=file.path(data.dir, "microarray_luma_filtered.RData"))
load(file=file.path(data.dir, "microarray_lumb_filtered.RData"))

MA    <- microarray.dat
MA.LA <- microarray.luma.dat
MA.LB <- microarray.lumb.dat

MA.LA.std <- MA.LA
for(r in 1:nrow(MA.LA)){
  MA.LA.std[r,] <- stdize(MA.LA[r,])
}

MA.LB.std <- MA.LB
for(r in 1:nrow(MA.LB)){
  MA.LB.std[r,] <- stdize(MA.LB[r,])
}

# Investigate mean vs standard deviation.
mean.ma.v <- apply(MA, 1, mean)
mad.ma.v  <- apply(MA, 1, mad)
sd.ma.v   <- apply(MA, 1, sd)

bad.genes.ma <- sd.ma.v/mad.ma.v > th

MA.LA.std.good <- MA.LA.std[!bad.genes.ma,]
MA.LB.std.good <- MA.LB.std[!bad.genes.ma,]


# Plot histogram with normal fit.
plot.test.stats <- FALSE
if(plot.test.stats){
  x <- DCM$it_test_stats[[4]]
  h<-hist(x, breaks=20)
  xfit<-seq(min(x), max(x), length.out=20) 
  yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
  yfit <- yfit*diff(h$mids[1:2])*length(x) 
  lines(xfit, yfit, col="blue", lwd=2)
}

#cor.LA    <- cor(MAM.LA.std)
#cor.LB    <- cor(MAM.LB.std)
#cor.basal <- cor(MAM.basal.std)


#epsilon <- 10e-05
#problem.idcs <- which(abs(DCM$it_test_stats[[3]]) <= epsilon)
#MAM.LA.std.old <- MAM.LA.std
#MAM.LB.std.old <- MAM.LB.std
#MAM.LA.std <- MAM.LA.std[-problem.idcs,]
#MAM.LB.std <- MAM.LB.std[-problem.idcs,]

# Running DCM on RNAseq
out.dir.1 <- file.path(out.dir, "No QR")
dir.create(out.dir.1)
setwd(out.dir.1)
lumA.lumB.50.noQR <- DCM_Kevin(MAM.LA.std,
                               MAM.LB.std, 
                               max.iter   = 10, 
                               max.time   = 10, 
                               max.groups = 3,
                               alpha      = .05,
                               est.size   = 50,
                               strict     = 'low', 
                               echo       = TRUE)
save(lumA.lumB.50.noQR, file=filePath(out.dir.1, "DCM.lumA.vs.lumB.noQR.50.RData"))

# Running DCM on microarray data
out.dir.2 <- file.path(out.dir, "Microarray No QR Filtered by SD over MAD")
dir.create(out.dir.2)
setwd(out.dir.2)
lumA.lumB.50.MA.noQR <- DCM_Kevin(MA.LA.std.good,
                                  MA.LB.std.good, 
                                  max.iter   = 10, 
                                  max.time   = 10, 
                                  max.groups = 3,
                                  alpha      = .05,
                                  est.size   = 50,
                                  strict     = 'low', 
                                  echo       = TRUE)
save(lumA.lumB.50.MA.noQR, file=filePath(out.dir.2, "DCM.lumA.vs.lumB.noQR.50.RData"))

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