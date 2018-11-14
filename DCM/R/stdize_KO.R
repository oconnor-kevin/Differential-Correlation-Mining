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