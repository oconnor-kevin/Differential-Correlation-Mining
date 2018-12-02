#' sanitize_DCM
#'
#' Takes data matrices and returns rows with na's, too many zeroes, too many 
#' ones, and low variances.
#'
#' @param M1 p x n1 data matrix for group 1.
#' @param M2 p x n2 data matrix for group 2.
#' @param strict String indicating whether low variance rows should be removed
#' (\code{"high"} or \code{"low"})
#' @param validation Boolean indicating whether validation data should be 
#' stored and returned.
#'
#' @return If \code{validation == TRUE}, list containing,
#' \itemize{
#'    \item killrows - Vectors of rows to be removed from M1 and M2.
#'    \item na1 - Vector of rows in M1 with na's.
#'    \item na2 - Vector of rows in M2 with na's.
#'    \item lb1 - Vector of rows in M1 with more than \code{ncol(M1)/10} zeroes.
#'    \item lb2 - Vector of rows in M2 with more than \code{ncol(M2)/10} zeroes.
#'    \item ub1 - Vector of rows in M1 with more than \code{ncol(M1)/10} ones.
#'    \item ub1 - Vector of rows in M2 with more than \code{ncol(M2)/10} ones.
#'    \item var1 - Vector of rows in M1 with variance less than 
#'    \code{1/(100*ncol(M1))}.
#'    \item var2 - Vector of rows in M2 with variance less than
#'    \cose{1/(100*ncol(M2))}.
#' }
#' Else, vector of rows to be removed.
#'
#' @export
sanitize_DCM <- function(M1, M2, strict = "high", validation = FALSE){
	p  <- nrow(M1)
	n1 <- ncol(M1)
	n2 <- ncol(M2)
		
	# Find rows with na's.
	na1 <- apply(M1, 1, function(x){sum(is.na(x)) > 0})	
	na2 <- apply(M2, 1, function(x){sum(is.na(x)) > 0})
	# Find rows with too many zeroes.
	lb1 <- apply(M1, 1, function(x){sum(x == 0) > n1/10})	
	lb2 <- apply(M2, 1, function(x){sum(x == 0) > n2/10})	
	# Find rows with too many ones.
	ub1 <- apply(M1, 1, function(x){sum(x == 1) > n1/10})	
	ub2 <- apply(M2, 1, function(x){sum(x == 1) > n2/10})	
	# Find rows with low variance.
	var1 <- apply(M1, 1, function(x){var(x) < 1/(100*n1)})
	var2 <- apply(M2, 1, function(x){var(x) < 1/(100*n2)})
		
	if(strict == "high"){
		killrows <- (1:p)[na1 | na2 | lb1 | lb2 | ub1 | ub2 | var1 | var2]
	}else{
		killrows <- (1:p)[na1 | na2 | lb1 | lb2 | ub1 | ub2]
	}

	if (validation){
	  return(list("killrows"=killrows, 
	              "na1"=na1, 
	              "na2"=na2, 
	              "lb1"=lb1, 
	              "lb2"=lb2, 
	              "ub1"=ub1, 
	              "ub2"=ub2, 
	              "var1"=var1, 
	              "var2"=var2))
	} else {
	  return(killrows)  
	}
}