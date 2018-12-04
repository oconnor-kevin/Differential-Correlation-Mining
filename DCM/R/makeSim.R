#' makeSim
#' 
#' Simulates two data matrices with pre-specified differentially correlated set.
#' 
#' @param p Number of variables.
#' @param k Size of differentially correlated set.
#' @param rho Correlation of set in first group.
#' @param n1 Number of columns in first data matrix.
#' @param n2 Number of columns in second data matrix.
#' @param bg String indicating what type of background should be used 
#' (\code{"n"}, \code{"pos"}, or \code{"noisy"}).
#' @param rho2 Correlation of set in second group.
#' @param varinf Boolean, if \code{varinf==TRUE}, sets diagonals of covariance
#' matrices to random uniform[0,5] values. Else, sets them to one.
#' @param large_k Integer describing how size of constant, positive submatrix to
#' add to covariance matrices.
#' 
#' @return List containing
#' \itemize{
#'     \item dat1 - Simulated p by n1 data matrix.
#'     \item dat2 - Simulated p by n2 data matrix.
#' }
#' 
#' @export
makeSim <- function(p, k, rho, n1, n2, bg = "n", rho2 = 0, varinf = FALSE, 
                    large_k = 0){
  require(MASS)
	# Initialize first covariance matrix.
	Sigma <- diag(p)
	# Add background to the covariance matrix.
	if(bg == "pos"){
		Sigma[] <- rho/3
	}else if(bg == "noisy"){
		temp <- array(rnorm(10*p), c(p, 10))
		Sigma <- cor(t(temp))
	}
	# Initialize second covariance matrix.
	Sigma2 <- Sigma
  # Add constant positive submatrix to both groups.
	if(large_k != 0){
		Sigma[1:large_k, 1:large_k] <- rho/2
		Sigma2[1:large_k, 1:large_k] <- rho/2
	}else{
		Sigma2[1:k, 1:k] <- rho2
	}
	# Create positively correlated subset in group 1.
	Sigma[1:k, 1:k] <- rho
  # Set diagonals of covariance matrix.
	if(varinf){
		vars <- runif(p, 0, 5)	
		Sigma[row(Sigma) == col(Sigma)] <- vars
		Sigma2[row(Sigma2) == col(Sigma2)] <- vars	
	}else{
		Sigma[row(Sigma) == col(Sigma)] <- 1
		Sigma2[row(Sigma2) == col(Sigma2)] <- 1
	}
	# Simulate data from specified covariance matrices.
	dat1 <- mvrnorm(n = n1, mu = rep(0, p), Sigma = Sigma)
	dat2 <- mvrnorm(n = n2, mu = rep(0, p), Sigma = Sigma2)
	
	return(list(dat1 = t(dat1), dat2 = t(dat2)))
}