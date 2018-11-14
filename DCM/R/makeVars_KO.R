# TITLE: makeVars_KO.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 9/10/18
# DATE MODIFIED: 9/10/18

makeVar_KO <- function(U_i, U_A){
  # Takes a vector U_i corresponding to the standardized values for a single 
  #  variable in a single group and a matrix U_A corresponding to the matrix of
  #  standardized values for the set of variables in the active set A and
  #  returns an estimate of the standard error of the test statistic.
  #
  # Args:
  # -U_i: Vector with mean=0 and sum of squares=1.
  # -U_A: Matrix with rowsums=0 and row sum of squares=1.
  # Returns:
  # -A numeric estimate of the standard error of the test statistic.
  
  n1 <- length(U_i)
  
  # Compute vector of correlations between variable i and the set A.
  r <- U_A %*% U_i
  r_A <- mean(r)
  
  W <- colMeans(U_A)
  Y <- as.vector(colMeans(t(r) %*% U_A^2))
  
  return(mean(r_A^2*U_i^4/4 - r_A*W*U_i^3 + (r_A*Y/2 + W^2)*U_i^2 - W*Y*U_i + Y^2/4))
}

makeVars_KO <- function(U_is, U_A){
  # Takes a matrix U_is corresponding to the standardized values for a set of 
  #  variables in a single group and a matrix U_A corresponding to the matrix of
  #  standardized values for the set of variables in the active set A and
  #  returns an estimate of the standard error of the test statistic.
  #
  # Args:
  # -U_is: Matrix with rowsums=0 and row sum of squares=1.
  # -U_A: Matrix with rowsums=0 and row sum of squares=1.
  # Returns:
  # -A vector of numeric estimates of the standard error of the test statistic.
  
  n1 <- ncol(U_is)
  
  # Compute matrix of correlations between variables in U_is and set A.
  r <- U_A %*% t(U_is)
  r_A <- colMeans(r)
  
  W <- colMeans(U_A)
  Y <- (t(r) %*% U_A^2)/nrow(U_A)
  
  r_A^2
}