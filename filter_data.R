# TITLE: filter_data.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 7/8/18
# DATE MODIFIED: 7/9/18

filter_data <- function(data.mat,
                        min.var=0,
                        max.var=1e20,
                        min.mean=-1e20,
                        max.mean=1e20,
                        transform=TRUE,
                        add.error=TRUE){
  require(dplyr)
  
  # Transforms and filters rows of a gene expression data matrix.
  #
  # Arguments:
  ## data.mat: Numeric data matrix to be filtered.
  ## min.var: Minimum variance to be kept.
  ## max.var: Maximum variance to be kept.
  ## min.mean: Minimum mean to be kept.
  ## max.mean: Maximum mean to be kept.
  ## transform: Logical indicating whether a log2 transform should be applied to
  ##  the data. If true, will also set all 0 values to the minimum value.
  ## add.error: Logical indicating whether some error should be added when 
  ##  setting 0 values to minimum value.
  #
  # Returns:
  ## Filtered numeric data matrix.
  
  out.data.mat <- data.mat
  
  # Get mean and variance of each row.
  dat.mean <- apply(out.data.mat, 1, mean)
  dat.var  <- apply(out.data.mat, 1, var)
  
  # Remove rows with na's.
  out.data.mat <- out.data.mat[apply(out.data.mat, 1, function(x){
      !any(is.na(x))
    }),]
  
  # Filter by mean and variance.
  good.idcs <- (dat.mean >= min.mean) & (dat.mean <= max.mean) & 
    (dat.var >= min.var) & (dat.var <= max.var)
  out.data.mat <- out.data.mat[good.idcs,]
  
  # Transform.
  if(transform){
    # Get minimum value.
    min.val <- min(apply(out.data.mat, 2, function(x){min(x[x!=0])}))
    
    if(add.error){
      # Get number of values approximately equal to 0.
      n.near.zero <- sum(near(out.data.mat, 0))
      
      # Generate random noise above 0.
      noise <- runif(n.near.zero, 2*min.val/3, min.val+2*min.val/3)
      
      # Add noise to values that are 0.
      out.data.mat[near(out.data.mat, 0)] <- noise
    }
    
    # Transform data.
    out.data.mat <- log2(out.data.mat)
    
    # Filter out rows that have variance close to 0 after transforming.
    dat.var.tr <- apply(out.data.mat, 1, var)
    out.data.mat <- out.data.mat[(dat.var.tr >= min.var) & 
                                   (dat.var.tr <= max.var),]
  }
}