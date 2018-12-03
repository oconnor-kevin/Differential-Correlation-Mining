#' makeVars
#'
#' Estimates variance of the DCM test statistic for variables outside the set A
#' using the estimator of Steiger and Hakstian (1982). More details are provided
#' in section 3.2 of the DCM paper.
#'
#' @param Uis Row-normalized matrix containing rows outside of set A for a 
#' single group.
#' @param M_A Row-normalized matrix containing rows inside of A for a single 
#' group.
#'
#' @return Vector giving variance estimates for each row outside of set A.
#'
#' @export
makeVars <- function(Uis, M_A){
	n1 <- ncol(Uis)
	
	# Cross correlations between variables in A and those not in A.
	r1s <- M_A %*% t(Uis)
	k <- length(r1s)
	
	# Make Uis standardized (sd=1) rather than normalized (sum of squares = 1).
	Uis <- Uis*sqrt(n1-1)
	
	# Compute variables in variance formula.
	W <- colMeans(M_A)*sqrt(n1-1)
	Y <- (t(r1s) %*% M_A^2)*(n1-1)/k
	rA <- mean(r1s)
	mat <- 1/4*rA^2*Uis^4 + rA*Y/2*Uis^2 + W^2*Uis^2 + Y^2/4 - W*Y*Uis - rA*W*Uis^3
	
	return(rowMeans(mat)/n1)
}