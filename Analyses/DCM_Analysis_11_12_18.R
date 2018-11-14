# TITLE: DCM_Analysis_11_12_18.R
# AUTHOR: Kevin O'Connor
# DATE CREATED:  11/12/18
# DATE MODIFIED: 11/12/18

# The purpose of this script is to evaluate the performance of the DCM 
#  algorithm with original variance estimation scheme in a null setting.

# Libraries and directories.
library(dplyr)
library(knitr)
library(R.utils)
library(dplyr)
library(limma)

sourceDirectory("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/DCM/R")
sourceDirectory("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/Utils")
data.dir <- file.path("/Users/kevinoconnor/Documents/Research/DCM/Data")
output.dir <- file.path("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/Analyses/11_12_18")

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
wh <- which(tar=="Her2")
wn <- which(tar=="Normal")
MAM.basal  <- dataM[, wb]
MAM.LA     <- dataM[, wla]
MAM.LB     <- dataM[, wlb]
MAM.her2   <- dataM[, wh]
MAM.normal <- dataM[, wn]

# Standardize row within each group.
stdize <- function(gene){
  # Make zero mean.
  gene <- gene - mean(gene)
  
  # Make sum of squares one.
  gene <- gene/sqrt(sum(gene^2))
  
  return(gene)
}

stdize_KO <- function(mat){
  # Standardize matrix of data, by rows.
  
  # Make mean 0.
  mat <- mat - rowMeans(mat)
  
  # Find sd of each row.
  sds <- sqrt(rowSums(mat*mat))
  
  # Set rows with zero standard deviation to 1.
  sds[which(sds==0)] <- 1
  
  # Return stdized matrix.
  return(mat/sds)
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

MAM.her2.std <- MAM.her2
for(r in 1:nrow(MAM.her2)){
  MAM.her2.std[r,] <- stdize(MAM.her2[r,])
}

MAM.normal.std <- MAM.normal
for(r in 1:nrow(MAM.normal)){
  MAM.normal.std[r,] <- stdize(MAM.normal[r,])
}

dat.luma   <- MAM.LA.std
dat.lumb   <- MAM.LB.std
dat.basal  <- MAM.basal.std
dat.her2   <- MAM.her2.std
dat.normal <- MAM.normal.std

# Define functions for variance estimation.
makeVars <- function(Uis, M_A){
  # stdized to ss = 1
  n1 = ncol(Uis)
  
  r1s = M_A %*% t(Uis)
  k = length(r1s)
  
  Uis = Uis*sqrt(n1-1)
  W = colMeans(M_A)*sqrt(n1-1)
  Y = (t(r1s) %*% M_A^2)*(n1-1)/k
  rA = mean(r1s)
  
  mat = 1/4*rA^2*Uis^4 + rA*Y/2*Uis^2 + W^2*Uis^2 + Y^2/4 - W*Y*Uis - rA*W*Uis^3
  
  return(rowMeans(mat)/n1)
}

makeVar <- function(Ui, M_A){
  # stdized to ss = 1
  
  n1 = length(Ui)
  k = dim(M_A)[1]
  
  r1s = M_A %*% Ui
  
  Ui = Ui*sqrt(n1-1)
  W = colMeans(M_A)*sqrt(n1-1)
  Y = (t(r1s) %*% M_A^2)*(n1-1)/k
  rA = mean(r1s)
  
  mat = 1/4*rA^2*Ui^4 + rA*Y/2*Ui^2 + W^2*Ui^2 + Y^2/4 - W*Y*Ui - rA*W*Ui^3
  
  
  return(mean(mat)/n1)
}

n1 <- ncol(dat.luma)
n2 <- ncol(dat.lumb)
p <- nrow(dat.luma)

################################################################################
# Experiment 1: Distribution of set sizes after a single iteration.
set.lengths <- c()
for (iter in 1:10){
  A <- sample.int(p, 1000)
  M1 <- rnorm(n1*p) %>% matrix(ncol=n1) %>% stdize_KO()
  M2 <- rnorm(n2*p) %>% matrix(ncol=n2) %>% stdize_KO()

  xA1 = M1[A,]
  xA2 = M2[A,]
  y1 = M1[-A,]
  y2 = M2[-A,]
  
  # Find mean vectors
  mean1 = colMeans(xA1)
  mean2 = colMeans(xA2)
  
  # Find norms
  n_m1 = sqrt(sum(mean1^2))
  n_m2 = sqrt(sum(mean2^2))
  
  # Length of A
  k = length(A)
  
  # Find test quantities for all variables
  corsm1 = t(cor(mean1, t(y1)))
  corsm2 = t(cor(mean2, t(y2)))
  
  # Test stat and variance
  obs = corsm1*n_m1 - corsm2*n_m2
  sd = sqrt(makeVars(y1, xA1) + makeVars(y2, xA2))
  
  # p-values for rows not in A
  test_out = pt(-obs/sd, min(c(n1-1, n2-1)), 0)
  
  ## Calculate p-values for rows in A
  # Adjust means to not include row
  mean1s = -t(t(xA1) - mean1*k)/(k-1)
  mean2s = -t(t(xA2) - mean2*k)/(k-1)
  
  # Find new mean norms
  n_m1s = sqrt(rowSums(mean1s*mean1s))
  n_m2s = sqrt(rowSums(mean2s*mean2s))
  
  # Find cors of rows with means
  corsm1 = rowMeans(stdize(mean1s)*xA1)*n1
  corsm2 = rowMeans(stdize(mean2s)*xA2)*n2
  
  # Make test stat and variance
  obss = corsm1*n_m1s - corsm2*n_m2s
  sds = sqrt(sapply(1:k, function(x) makeVar(xA1[x,], xA1[-x,]) + makeVar(xA2[x,], xA2[-x,])))
  
  # Find pvals
  test_in = pt(-obss/sds, min(c(n1-1, n2-1)), 0)
  
  # Combine all p-values
  test = P
  test[-A] = test_out
  test[A] = test_in
  
  # Update A to include significant rows
  newA = bhy(test, alpha = alpha)
  
  set.lengths <- c(set.lengths, length(newA))  
}

################################################################################
# Experiment 2: Distribution of set sizes after a single iteration with a few 
#  bad genes.
# Conclusion: Bad genes of this type don't seem to affect the algorithm like 
#  we've seen.
set.lengths <- c()
for (iter in 1:100){
  A <- sample.int(p, 50)
  M1 <- rnorm(n1*p) %>% matrix(ncol=n1)
  M2 <- rnorm(n2*p) %>% matrix(ncol=n2)
    
  # Set a few of the genes values to 0.
  bad.genes <- sample.int(p, 3340)
  M1[bad.genes, sample.int(n1, 200)] <- min(min(M1), min(M2))
  M2[bad.genes, sample.int(n2, 100)] <- min(min(M1), min(M2)) + 10
  
  M1 <- stdize_KO(M1)
  M2 <- stdize_KO(M2)
  
  xA1 = M1[A,]
  xA2 = M2[A,]
  y1 = M1[-A,]
  y2 = M2[-A,]
  
  # Find mean vectors
  mean1 = colMeans(xA1)
  mean2 = colMeans(xA2)
  
  # Find norms
  n_m1 = sqrt(sum(mean1^2))
  n_m2 = sqrt(sum(mean2^2))
  
  # Length of A
  k = length(A)
  
  # Find test quantities for all variables
  corsm1 = t(cor(mean1, t(y1)))
  corsm2 = t(cor(mean2, t(y2)))
  
  # Test stat and variance
  obs = corsm1*n_m1 - corsm2*n_m2
  sd = sqrt(makeVars(y1, xA1) + makeVars(y2, xA2))
  
  # p-values for rows not in A
  test_out = pt(-obs/sd, min(c(n1-1, n2-1)), 0)
  
  ## Calculate p-values for rows in A
  # Adjust means to not include row
  mean1s = -t(t(xA1) - mean1*k)/(k-1)
  mean2s = -t(t(xA2) - mean2*k)/(k-1)
  
  # Find new mean norms
  n_m1s = sqrt(rowSums(mean1s*mean1s))
  n_m2s = sqrt(rowSums(mean2s*mean2s))
  
  # Find cors of rows with means
  corsm1 = rowMeans(stdize(mean1s)*xA1)*n1
  corsm2 = rowMeans(stdize(mean2s)*xA2)*n2
  
  # Make test stat and variance
  obss = corsm1*n_m1s - corsm2*n_m2s
  sds = sqrt(sapply(1:k, function(x) makeVar(xA1[x,], xA1[-x,]) + makeVar(xA2[x,], xA2[-x,])))
  
  # Find pvals
  test_in = pt(-obss/sds, min(c(n1-1, n2-1)), 0)
  
  # Combine all p-values
  test = P
  test[-A] = test_out
  test[A] = test_in
  
  # Update A to include significant rows
  newA = bhy(test, alpha = alpha)
  
  set.lengths <- c(set.lengths, length(newA))  
}
hist(set.lengths)

################################################################################
# Count number of multimodal genes.
n.multmod <- dat.luma %>% apply(1, function(r){
  Mclust(r)$G > 1
}) %>% sum() 
# Approximately 3300 / 4100 = 80%

################################################################################
# Experiment 3: Making some of the individuals correlated.
# Conclusion: Correlation between individuals (columns) can give bowl-shaped
#  p-value distributions and large sets in a null setting.
A <- sample.int(p, 50)
M1 <- rnorm(n1*p) %>% matrix(ncol=n1)
M2 <- rnorm(n2*p) %>% matrix(ncol=n2)

# Add correlation between columns.
for (cor.iter in 1:50){
  i <- sample.int(n1)
  j <- sample.int(n1)
  k <- sample.int(n1)
  M1[,i] <- M1[,i] + 0.05*M1[,j] + 0.1*M1[,k]
  
  i <- sample.int(n2)
  j <- sample.int(n2)
  k <- sample.int(n2)
  M2[,i] <- M2[,i] + 0.05*M2[,j] + 0.1*M2[,k]
}

M1 <- stdize_KO(M1)
M2 <- stdize_KO(M2)

xA1 = M1[A,]
xA2 = M2[A,]
y1 = M1[-A,]
y2 = M2[-A,]

# Find mean vectors
mean1 = colMeans(xA1)
mean2 = colMeans(xA2)

# Find norms
n_m1 = sqrt(sum(mean1^2))
n_m2 = sqrt(sum(mean2^2))

# Length of A
k = length(A)

# Find test quantities for all variables
corsm1 = t(cor(mean1, t(y1)))
corsm2 = t(cor(mean2, t(y2)))

# Test stat and variance
obs = corsm1*n_m1 - corsm2*n_m2
sd = sqrt(makeVars(y1, xA1) + makeVars(y2, xA2))

# p-values for rows not in A
test_out = pt(-obs/sd, min(c(n1-1, n2-1)), 0)

## Calculate p-values for rows in A
# Adjust means to not include row
mean1s = -t(t(xA1) - mean1*k)/(k-1)
mean2s = -t(t(xA2) - mean2*k)/(k-1)

# Find new mean norms
n_m1s = sqrt(rowSums(mean1s*mean1s))
n_m2s = sqrt(rowSums(mean2s*mean2s))

# Find cors of rows with means
corsm1 = rowMeans(stdize(mean1s)*xA1)*n1
corsm2 = rowMeans(stdize(mean2s)*xA2)*n2

# Make test stat and variance
obss = corsm1*n_m1s - corsm2*n_m2s
sds = sqrt(sapply(1:k, function(x) makeVar(xA1[x,], xA1[-x,]) + makeVar(xA2[x,], xA2[-x,])))

# Find pvals
test_in = pt(-obss/sds, min(c(n1-1, n2-1)), 0)

# Combine all p-values
test = P
test[-A] = test_out
test[A] = test_in

# Update A to include significant rows
newA = bhy(test, alpha = alpha)

pdf(file=file.path(output.dir, "null_corr_pvals.pdf"))
hist(test, main="Null P-values with Correlated Columns")
dev.off()

################################################################################
# Experiment 4: Remove correlation between columns via PCA then run the alg.
pc.luma <- prcomp(dat.luma, center=TRUE, scale=TRUE)
pc.lumb <- prcomp(dat.lumb, center=TRUE, scale=TRUE)
pc.dat.luma <- pc.luma$x
pc.dat.lumb <- pc.lumb$x

A <- sample.int(p, 50)
M1 <- pc.dat.luma %>% stdize_KO()
M2 <- pc.dat.lumb %>% stdize_KO()

xA1 = M1[A,]
xA2 = M2[A,]
y1 = M1[-A,]
y2 = M2[-A,]

# Find mean vectors
mean1 = colMeans(xA1)
mean2 = colMeans(xA2)

# Find norms
n_m1 = sqrt(sum(mean1^2))
n_m2 = sqrt(sum(mean2^2))

# Length of A
k = length(A)

# Find test quantities for all variables
corsm1 = t(cor(mean1, t(y1)))
corsm2 = t(cor(mean2, t(y2)))

# Test stat and variance
obs = corsm1*n_m1 - corsm2*n_m2
sd = sqrt(makeVars(y1, xA1) + makeVars(y2, xA2))

# p-values for rows not in A
test_out = pt(-obs/sd, min(c(n1-1, n2-1)), 0)

## Calculate p-values for rows in A
# Adjust means to not include row
mean1s = -t(t(xA1) - mean1*k)/(k-1)
mean2s = -t(t(xA2) - mean2*k)/(k-1)

# Find new mean norms
n_m1s = sqrt(rowSums(mean1s*mean1s))
n_m2s = sqrt(rowSums(mean2s*mean2s))

# Find cors of rows with means
corsm1 = rowMeans(stdize(mean1s)*xA1)*n1
corsm2 = rowMeans(stdize(mean2s)*xA2)*n2

# Make test stat and variance
obss = corsm1*n_m1s - corsm2*n_m2s
sds = sqrt(sapply(1:k, function(x) makeVar(xA1[x,], xA1[-x,]) + makeVar(xA2[x,], xA2[-x,])))

# Find pvals
test_in = pt(-obss/sds, min(c(n1-1, n2-1)), 0)

# Combine all p-values
test = P
test[-A] = test_out
test[A] = test_in

# Update A to include significant rows
newA = bhy(test, alpha = alpha)

pdf(file=file.path(output.dir, "PC_Pvals.pdf"))
hist(test, main="P-values of RNAseq data when columns replaced by PC's")
dev.off()

################################################################################
# Experiment 5: Getting p-value distribution for RNAseq data.
A <- sample.int(p, 50)
M1 <- dat.luma %>% stdize_KO()
M2 <- dat.lumb %>% stdize_KO()

xA1 = M1[A,]
xA2 = M2[A,]
y1 = M1[-A,]
y2 = M2[-A,]

# Find mean vectors
mean1 = colMeans(xA1)
mean2 = colMeans(xA2)

# Find norms
n_m1 = sqrt(sum(mean1^2))
n_m2 = sqrt(sum(mean2^2))

# Length of A
k = length(A)

# Find test quantities for all variables
corsm1 = t(cor(mean1, t(y1)))
corsm2 = t(cor(mean2, t(y2)))

# Test stat and variance
obs = corsm1*n_m1 - corsm2*n_m2
sd = sqrt(makeVars(y1, xA1) + makeVars(y2, xA2))

# p-values for rows not in A
test_out = pt(-obs/sd, min(c(n1-1, n2-1)), 0)

## Calculate p-values for rows in A
# Adjust means to not include row
mean1s = -t(t(xA1) - mean1*k)/(k-1)
mean2s = -t(t(xA2) - mean2*k)/(k-1)

# Find new mean norms
n_m1s = sqrt(rowSums(mean1s*mean1s))
n_m2s = sqrt(rowSums(mean2s*mean2s))

# Find cors of rows with means
corsm1 = rowMeans(stdize(mean1s)*xA1)*n1
corsm2 = rowMeans(stdize(mean2s)*xA2)*n2

# Make test stat and variance
obss = corsm1*n_m1s - corsm2*n_m2s
sds = sqrt(sapply(1:k, function(x) makeVar(xA1[x,], xA1[-x,]) + makeVar(xA2[x,], xA2[-x,])))

# Find pvals
test_in = pt(-obss/sds, min(c(n1-1, n2-1)), 0)

# Combine all p-values
test = P
test[-A] = test_out
test[A] = test_in

# Update A to include significant rows
newA = bhy(test, alpha = alpha)

pdf(file=file.path(output.dir, "RNAseq_Pvals.pdf"))
hist(test, main="P-values of RNAseq data")
dev.off()



################################################################################
# Experiment 6: Null distribution of p-values.
A <- sample.int(p, 50)
M1 <- rnorm(n1*p) %>% matrix(ncol=n1) %>% stdize_KO()
M2 <- rnorm(n2*p) %>% matrix(ncol=n2) %>% stdize_KO()

xA1 = M1[A,]
xA2 = M2[A,]
y1 = M1[-A,]
y2 = M2[-A,]

# Find mean vectors
mean1 = colMeans(xA1)
mean2 = colMeans(xA2)

# Find norms
n_m1 = sqrt(sum(mean1^2))
n_m2 = sqrt(sum(mean2^2))

# Length of A
k = length(A)

# Find test quantities for all variables
corsm1 = t(cor(mean1, t(y1)))
corsm2 = t(cor(mean2, t(y2)))

# Test stat and variance
obs = corsm1*n_m1 - corsm2*n_m2
sd = sqrt(makeVars(y1, xA1) + makeVars(y2, xA2))

# p-values for rows not in A
test_out = pt(-obs/sd, min(c(n1-1, n2-1)), 0)

## Calculate p-values for rows in A
# Adjust means to not include row
mean1s = -t(t(xA1) - mean1*k)/(k-1)
mean2s = -t(t(xA2) - mean2*k)/(k-1)

# Find new mean norms
n_m1s = sqrt(rowSums(mean1s*mean1s))
n_m2s = sqrt(rowSums(mean2s*mean2s))

# Find cors of rows with means
corsm1 = rowMeans(stdize(mean1s)*xA1)*n1
corsm2 = rowMeans(stdize(mean2s)*xA2)*n2

# Make test stat and variance
obss = corsm1*n_m1s - corsm2*n_m2s
sds = sqrt(sapply(1:k, function(x) makeVar(xA1[x,], xA1[-x,]) + makeVar(xA2[x,], xA2[-x,])))

# Find pvals
test_in = pt(-obss/sds, min(c(n1-1, n2-1)), 0)

# Combine all p-values
test = P
test[-A] = test_out
test[A] = test_in

# Update A to include significant rows
newA = bhy(test, alpha = alpha)

pdf(file=file.path(output.dir, "null_pvals.pdf"))
hist(test, main="P-values for null data")
dev.off()