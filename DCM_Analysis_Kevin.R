# TITLE: DCM_Analysis_Kevin.R
# AUTHOR: Kevin O'Connor(, Di Wu)
# DATE MODIFIED: 6/22/18

# Switches.
is.test <- FALSE
run.small.sample.test <- FALSE

# Libraries and directories.
library(R.utils)
library(dplyr)
library(limma)

sourceDirectory("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/DCM/R")
if(is.test){
  setwd("/Users/kevinoconnor/Documents/Research/DCM/Test")  
} else {
  setwd("/Users/kevinoconnor/Documents/Research/DCM")
}

# Reading data.
if(!(exists("SeqData") & exists("sampleTab"))){
  SeqData   <- read.delim("BRCA.1201.FINAL.Cluster.txt" , header=F)
  sampleTab <- read.delim("BRCA.1218_pam50scores.FINAL.txt")
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


# Running DCM
## Create output directory.
out.dir <- filePath(getwd(), gsub("-", "_", Sys.time()))
dir.create(out.dir)
setwd(out.dir)

## Test run for small number of genes.
if(run.small.sample.test){
  a <- DCM_Kevin(MAM.basal[1:500, ], 
                 MAM.LA[1:500, ], 
                 max.iter = 10, 
                 max.time = 100, 
                 alpha = .05,  
                 strict='low', 
                 echo=TRUE)
  ### max.iter = 10 can be increased from 10 to 20, 30…. but it will take longer.
  save(a, file=filePath(out.dir,"SmallSampleTest.RData"))
}
  
## LumA vs LumB, No QR (Quantile Normalization ?)
lumA.lumB.50.noQR <- DCM_Kevin(MAM.LA,
                               MAM.LB, 
                               max.iter = 10, 
                               max.time = 100, 
                               alpha    = .05,
                               est.size = 50,
                               strict   = 'low', 
                               echo     = TRUE)
save(lumA.lumB.50.noQR, file=filePath(out.dir, "DCM.lumA.vs.lumB.noQR.50.RData"))

lumA.lumB.200.noQR <- DCM_Kevin(MAM.LA,
                                MAM.LB, 
                                max.iter = 10, 
                                max.time = 100, 
                                alpha    = .05,
                                est.size = 200,
                                strict   = 'low', 
                                echo     = TRUE)
save(lumA.lumB.200.noQR, file=filePath(out.dir, "DCM.lumA.vs.lumB.noQR.200.RData"))

lumA.lumB.125.noQR <- DCM_Kevin(MAM.LA,
                                MAM.LB, 
                                max.iter = 10, 
                                max.time = 100, 
                                alpha    = .05,
                                est.size = 125,
                                strict   = 'low', 
                                echo     = TRUE)
save(lumA.lumB.125.noQR, file=filePath(out.dir, "DCM.lumA.vs.lumB.noQR.125.RData"))

#[1] "Removing 4697 rows - too many zero expressions."
#[1] "Warning: Data is non-normal; consider quantile normalizing.  (Set QN = TRUE)"
#[1] "Done checking data."
#[1] "Overall cor, group 1: 0.014379"
#[1] "Overall cor, group 2: 0.010060"

DCM.lumA.vs.basal.noQR=DCM(MAM.LA,
                           MAM.basal, 
                           max.iter = 10, 
                           max.time = 100, 
                           alpha    = .05,
                           est.size = 50,
                           strict   = 'low', 
                           echo     = TRUE)
save(DCM.lumA.vs.basal.noQR, file=filePath(out.dir, "DCM.lumA.vs.basal.noQR.RData"))
#save(DCM.lumA.vs.basal.noQR, file="/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/RData/DCM.lumA.vs.basal.noQR.RData" )

DCM.lumB.vs.lumA.noQR=DCM(MAM.LB, 
                          MAM.LA, 
                          max.iter = 10, 
                          max.time = 100, 
                          alpha = .05,  
                          strict='low', 
                          echo = TRUE)
save(DCM.lumB.vs.lumA.noQR, file=filePath(out.dir, "DCM.lumB.vs.lumA.noQR.RData"))
#save(DCM.lumB.vs.lumA.noQR, file="/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/RData/DCM.lumB.vs.lumA.noQR.RData" )

DCM.basal.vs.lumA=DCM(MAM.basal,
                      MAM.LA, 
                      max.iter = 10, 
                      max.time = 100, 
                      alpha = .05,  
                      strict='low', 
                      echo = T, 
                      QN = TRUE, 
                      resid.full = TRUE)
save(DCM.basal.vs.lumA, file=filePath(out.dir, "DCM.basal.vs.lumA.RData"))
#save(DCM.basal.vs.lumA, file="/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/RData/DCM.basal.vs.lumA.RData" )

DCM.basal.vs.lumA[[3]]
DCM.basal.vs.lumA[[5]]

DCM.lumA.vs.lumB=DCM(MAM.LA,
                     MAM.LB, 
                     max.iter = 10, 
                     max.time = 100, 
                     alpha = .05,  
                     strict='low', 
                     echo = TRUE, 
                     QN = TRUE , 
                     resid.full = TRUE)
save(DCM.lumA.vs.lumB, file=filePath(out.dir, "DCM.lumA.vs.lumB.RData"))
#save(DCM.lumA.vs.lumB, file="/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/DCM.lumA.vs.lumB.RData" )

#length(DCM.lumA.vs.lumB[[1]])
#length(DCM.lumA.vs.lumB[[1]][[1]])
#length(DCM.lumA.vs.lumB[[1]][[2]])

length(DCM.lumA.vs.lumB[[1]][[2]])
DCM.lumA.vs.lumB[[2]] 
#g32=fit.can[canF,] $genes$symbol [DCM.lumA.vs.lumB[[1]][[1]]]
# g32[order(g32)]

#g99=fit.can[canF,] $genes$symbol [DCM.lumA.vs.lumB[[1]][[2]]]
#g99[order(g99)]

###

DCM.lumA.vs.basal=DCM(MAM.LA,
                      MAM.basal, 
                      max.iter = 10, 
                      max.time = 100, 
                      alpha = .05,  
                      strict='low', 
                      echo=TRUE, 
                      QN=TRUE, 
                      resid.full=TRUE)
save(DCM.lumA.vs.basal, file=filePath(out.dir, "DCM.lumA.vs.basal.RData"))
#save(DCM.lumA.vs.basal, file="/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/DCM.lumA.vs.basal.RData" )


DCM.lumB.vs.lumA=DCM(MAM.LB, 
                     MAM.LA, 
                     max.iter = 10, 
                     max.time = 100, 
                     alpha = .05,  
                     strict='low', 
                     echo=TRUE, 
                     QN=TRUE, 
                     resid.full=TRUE)
save(DCM.lumB.vs.lumA, file=filePath(out.dir, "DCM.lumB.vs.lumA.RData"))
#save(DCM.lumB.vs.lumA, file="/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/DCM.lumB.vs.lumA.RData" )

#####


#setwd("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/RData")
list.files()

for (i in 1:8){load(list.files()[i])}
list.wQR  <- list(DCM.basal.vs.lumA, DCM.lumA.vs.basal, DCM.lumA.vs.lumB, DCM.lumB.vs.lumA)
list.noQR <- list(DCM.basal.vs.lumA.noQR, DCM.lumA.vs.basal.noQR, DCM.lumA.vs.lumB.noQR, DCM.lumB.vs.lumA.noQR)
# list.noQR<-list(  DCM.basal.vs.lumA.noQR , DCM.lumA.vs.basal.noQR, DCM.lumA.vs.lumB.noQR, DCM.lumB.vs.lumA.noQR)
list.full <- list(DCM.basal.vs.lumA ,DCM.basal.vs.lumA.noQR,DCM.lumA.vs.basal, DCM.lumA.vs.basal.noQR,DCM.lumA.vs.lumB, DCM.lumA.vs.lumB.noQR,DCM.lumB.vs.lumA , DCM.lumB.vs.lumA.noQR)
names(list.full)<-c("basal.vs.lumA", "basal.vs.lumA.noQR",
                    "lumA.vs.basal",  "lumA.vs.basal.noQR", 
                    "lumA.vs.lumB",  "lumA.vs.lumB.noQR",
                    "lumB.vs.lumA",   "lumB.vs.lumA.noQR")
lapply(list.full, function(x)    unlist(lapply(x[[1]], length)   )   )
lapply(list.full, function(x)    unlist( x[[4]])   )
lapply(list.full, function(x)    unlist(x[[5]]) )
lapply(list.full, function(x)    unlist(x[[2]]) ) # iteration

#table(duplicated (fit$gene))
#FALSE  TRUE 
#20504    32
#fit$gene[duplicated (fit$gene)]
#[1] "CLID"    "EWEIGHT" "EWEIGHT" "?"       "?"       "?"       "?"      
#[8] "?"       "?"       "?"       "?"       "?"       "?"       "?"      
#[15] "?"       "?"       "?"       "?"       "?"       "?"       "?"      
#[22] "?"       "?"       "?"       "SLC35E2" "?"       "?"       "?"      
#[29] "?"       "?"       "?"       "?"      

dcm.brAr<-lapply(  list.full, function(x)    lapply(x[[1]], function(y)   fit$gene[y]    )    )
totGF<-fit$gene
#load("/Users/diwu/Dropbox/MacPro2015Documents/Harvard/YiFang/MsigDBv2.5/human_c2_v2.5.rdata")
iset.msig2<-NULL
for (k in 1:length(Hs.gmtl.c2)){
  iset.msig2[[k]]<-which(totGF%in%as.character(Hs.gmtl.c2[[k]]))
}
names(iset.msig2)<- names(Hs.gmtl.c2)
#phyper(x,m,n,k, lower.tail=F)+dhyper(x,m,n,k)
# pTest<- array(c(vector1,vector2),dim = c(  length(Hs.gmtl.c2) ,1,28))

pL<-NULL
for(i in 1:length(dcm.brAr)){
  #print(i)
  pTest<-matrix(99, ncol= length(dcm.brAr[[i]]), nrow= length(Hs.gmtl.c2)  )
  for (j in 1:length(dcm.brAr[[i]])) {
    
    aset<- dcm.brAr[[i]][[j]]
    m=length(aset)
    n=length(totGF)-m
    
    #print(paste("j", j)) 
    
    for (kk in 1:length(iset.msig2))
    {
      
      #print(paste("kk", kk)) 
      
      k=length(iset.msig2[[kk]])
      bset<-totGF[(iset.msig2[[kk]])]
      ovset<-bset[bset %in% aset]
      # print(ovset)
      x<-length(ovset)
      
      pTest[kk,j]<-phyper(x,m,n,k, lower.tail=F)+dhyper(x,m,n,k)
    }
  }
  
  pL[[ i ]] <- pTest
  
  
} 

table(p.adjust(pL[[1]], method="BH")   <0.05)#60


names(Hs.gmtl.c2) [ x  <0.05 ]


dsL<-NULL
for (i in 1:length(dcm.brAr))
{
  
  ds<-NULL
  for (j in 1:length(dcm.brAr[[i]])) {
    
    #x= p.adjust( pL[[i]][,j]  , method="BH")
    
    x= pL[[i]][,j]
    ds [[j]]<-names(Hs.gmtl.c2) [ x  <0.05]
    
    
    
  }
  
  dsL[[i]] <- ds
}
names(dsL ) <- names(full.list)

lapply(dsL, length)
lapply(dsL,function(x) lapply(x, length) )






#library(gplots)
#heatmap.2(data ,  dendrogram ="both",col=greenred(32), cexRow=.01,trace="none", labCol=c(rep("A", dim(MAM.LA)[2]),rep("B", dim(MAM.LA)[2]) )  )

library(reshape2)
library(ggplot2)
library("cowplot")



#  

n.basal=dim(MAM.basal)[2]
n.LA=dim(MAM.LA)[2]
n.LB=dim(MAM.LB)[2]


corPlot<-function(data, n.basal=dim(MAM.basal)[2], n.LA=dim(MAM.LA)[2], n.LB=dim(MAM.LB)[2] ){
  
  if( i%in%c(1:4)) {
    CA.a<- cor(t(data[, 1:n.basal]))
    CA.b<- cor(t(data[, (n.basal+1):(n.basal+n.LA) ]  ))
  }
  
  if( i%in%c(5:8)) {
    CA.a<- cor(t(data[, 1:n.LA]))
    CA.b<- cor(t(data[, (n.LA+1):(n.LA+n.LB) ]  ))
  }
  
  
  nc<- dim(CA.a)[2]
  nc1=nc-1
  corA<-(sum(CA.a)-nc) /(nc1*nc1)
  #[1] 0.108802
  corB= (sum(CA.b)-nc) /(nc1*nc1)
  #[1] 0.4850927
  
  print(paste("ngenes",nc, "nsam", dim(data)[2]))
  
  
  if( i%in%c(1,2,5,6)) {
    print(paste("setA", corA, "setB", corB)) }
  
  
  if( i%in%c(3,4,7,8)) {
    print(paste("setA", corB, "setB", corA)) }
  
  
  melted_cormatA <- melt(CA.a)
  #head(melted_cormatA)
  p1=ggplot(data = melted_cormatA, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile()
  melted_cormatB <- melt(CA.b)
  #head(melted_cormatB)
  p2=ggplot(data = melted_cormatB, aes(x=Var1, y=Var2, fill=value), title=paste("setA", dim(data)[1], ".pdf", sep="")) +   geom_tile()
  return(list(p1,p2))
}





tD.basA=cbind( MAM.basal, MAM.LA)
tD.AB=cbind(MAM.LA, MAM.LB)

a<-nm<-NULL
for (i in 1:8){
  xx=list.full[[i]][[1]]
  
  if( i%in%c(1:4)) {
    dataS<-lapply(   xx, function(x)  tD.basA[ x , ]    )
  }  
  if( i%in%c(5,6,7,8)) {
    dataS<-lapply(   xx, function(x)  tD.AB[ x , ]    )
  } 
  
  
  a[[i]]= lapply( dataS, function(x) corPlot(x)  )
  nm[[i]]= paste (names(list.full)[i], ".pdf", sep="")
}





i=1
pdf (  paste("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/corHeat/", nm[[i]], sep=""),   height=6, width=6)
a[[i]]
dev.off()
i=2
pdf (  paste("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/corHeat/", nm[[i]], sep=""),   height=6, width=6)
a[[i]]
dev.off()
i=3
pdf (  paste("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/corHeat/", nm[[i]], sep=""),   height=6, width=6)
a[[i]]
dev.off()
i=4
pdf (  paste("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/corHeat/", nm[[i]], sep=""),   height=6, width=6)
a[[i]]
dev.off()
i=5
pdf (  paste("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/corHeat/", nm[[i]], sep=""),   height=6, width=6)
a[[i]]
dev.off()
i=6
pdf (  paste("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/corHeat/", nm[[i]], sep=""),   height=6, width=6)
a[[i]]
dev.off()
i=7
pdf (  paste("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/corHeat/", nm[[i]], sep=""),   height=6, width=6)
a[[i]]
dev.off()
i=8
pdf (  paste("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/corHeat/", nm[[i]], sep=""),   height=6, width=6)
a[[i]]
dev.off()


data.lumB.vs.lumA<-lapply(   DCM.lumB.vs.lumA[[1]], function(x)  tD.basA[ x , ]    )
a= lapply( data.lumB.vs.lumA, function(x) corPlot(x)  )
nm="lumB.vs.lumA.pdf"
pdf (  paste("/Users/diwu/Dropbox/UNCwork/projects/GTEX2018/br/BR_RNAseq/corHeat/", nm, sep=""),   height=6, width=6)
a
dev.off()
















