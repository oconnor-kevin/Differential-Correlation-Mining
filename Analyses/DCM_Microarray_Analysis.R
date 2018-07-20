# TITLE: DCM_Microarray_Analysis.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 7/16/18
# DATE MODIFIED: 7/16/18

# Cleanup.
rm(list=ls())

# Switches.
is.test <- FALSE

# Libraries and directories.
library(R.utils)
library(dplyr)
library(limma)

sourceDirectory("/Users/kevinoconnor/Documents/Research/DCM/Differential-Correlation-Mining/DCM/R", modifiedOnly=FALSE)
if(is.test){
  out.root <- file.path("/Users/kevinoconnor/Documents/Research/DCM/Test/Data")  
} else {
  out.root <- file.path("/Users/kevinoconnor/Documents/Research/DCM/Data")
}
data.dir <- file.path("/Users/kevinoconnor/Documents/Research/DCM/Data")
out.dir <- file.path(out.root, gsub("-", "_", Sys.time()))
dir.create(out.dir, showWarnings=FALSE)
setwd(out.dir)

# Read data.
load(file=file.path(data.dir, "microarray_luma_filtered.RData"))
load(file=file.path(data.dir, "microarray_lumb_filtered.RData"))

MA.LA <- microarray.luma.dat
MA.LB <- microarray.lumb.dat

# Running DCM.
DCM.LumA.LumB.MA <- DCM_Kevin(MA.LA,
                              MA.LB, 
                              max.iter   = 10, 
                              max.time   = 100, 
                              alpha      = .05,
                              est.size   = 1500,
                              strict     = 'low', 
                              echo       = TRUE,
                              QN         = TRUE, 
                              resid.full = TRUE)
save(DCM.LumA.LumB.MA, file=filePath(out.dir, "DCM_A_B_1500_wQR.RData"))
