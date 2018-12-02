#' stdize
#'
#' Standardize matrix to have mean=0 and sum of squares=1, by rows.
#'
#' @param mat Data matrix to be row-standardized.
#'
#' @return Row-standardized matrix.
#'
#' @export
stdize <- function(Mat){
	# Subtract mean of each row
	Mat <- Mat - rowMeans(Mat)
	# Find sd of each row
	sds <- sqrt(rowSums(Mat*Mat))

	return(Mat/sds)
}