# TITLE: DCM_Analysis_PVals_2018_07_26.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 7/20/18
# DATE MODIFIED: 7/26/18

# This script plots the p-values that were used in running the DCM algorithm on
#  the microarray data from Di Wu.

# Libraries and directories.
library(dplyr)
data.dir <- file.path("/Users/kevinoconnor/Documents/Research/DCM/Data/7_26_18/Data_2")

# Collect data from each search.
p_vals_l_l <- lapply(1:2, function(i){
  # Load data from i'th search.
  load(file.path(data.dir, paste0("DCM_", i, ".RData")))
  
  # Save p-values from each iteration.
  p_vals_l <- lapply(1:length(DCM$it_p_vals), function(j){
    return(DCM$it_p_vals[[j]])
  })
  
  return(p_vals_l)
})

# Plotting p-values.
pdf(file.path(data.dir, "Microarray_PVals.pdf"))
for(i in 1:2){
  for(j in 1:length(p_vals_l_l[[i]])){
    hist(p_vals_l_l[[i]][[j]], 
         main=paste0("Microarray Data Search ", i, ", Iteration ", j),
         xlab="p value")
  }
}
dev.off()