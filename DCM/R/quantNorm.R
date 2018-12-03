#' quantNorm
#' 
#' Quantile normalize rows of a given matrix.
#' 
#' @param Mat Matrix to be quantile normalized.
#' 
#' @return Matrix with rows quantile normalized.
#' 
#' @export
quantNorm <- function(Mat){
	# Find row length
	n <- dim(Mat)[2]
	# Rank each row
	ranks <- t(apply(Mat, 1, rank))
  # Store quantile-normalized data.
	Mat <- matrix(qnorm(ranks/(n+1)), ncol = n)
	
	return(Mat)
}
