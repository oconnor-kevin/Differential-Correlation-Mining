stdize_KO <- function(mat){
  # Standardize matrix of data, by rows.
  
  # Make mean 0.
  mat <- mat - rowMeans(mat)
  
  # Find sd of each row.
  sds <- sqrt(rowSums(mat*mat))
    
  # Return stdized matrix.
  return(mat/sds)
}