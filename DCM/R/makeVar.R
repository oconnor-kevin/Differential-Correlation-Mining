#' makeVar
#'
#' Estimates variance of the DCM test statistic for a variable inside the set A
#' using the estimator of Steiger and Hakstian (1982). More details are provided
#' in section 3.2 of the DCM paper.
#'
#' @param Ui Normalized row in A for a single group.
#' @param M_A Row-normalized matrix containing rows inside of A excluding Uifor 
#' a single group.
#'
#' @return Number indicating variance of variable Ui.
#'
#' @export
makeVar <- function(Ui, M_A){
  n1 <- length(Ui)
  k <- dim(M_A)[1]
  
  # Get correlation of Ui with other variables in A.
  r1s <- M_A %*% Ui
  
  # Make Ui and M_Astandardized (sd=1) rather than normalized 
  # (sum of squares = 1).
  Ui <- Ui*sqrt(n1-1)
  W <- colMeans(M_A)*sqrt(n1-1)
  
  Y <- (t(r1s) %*% M_A^2)*(n1-1)/k
  rA <- mean(r1s)
  mat <- 1/4*rA^2*Ui^4 + rA*Y/2*Ui^2 + W^2*Ui^2 + Y^2/4 - W*Y*Ui - rA*W*Ui^3
  
  return(mean(mat)/n1)
}