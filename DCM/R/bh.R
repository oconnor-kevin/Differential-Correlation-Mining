#' bh
#'
#' Given a vector of p-values, selects a subset to reject based on the 
#' Benjamini-Hochberg FDR control procedure.
#'
#' @param pvals Vector of p-values to be tested.
#' @param alpha Number between 0 and 1 indicating level at which to control FDR.
#'
#' @return Indices of p-values to be rejected.
#'
#' @export
bh <- function(pvals, alpha = 0.05){
  # Sorted p-vals
  sp <- sort(pvals)
  
  # Save original order of p-vals
  ord <- order(pvals)
  
  # Find bh cutoff
  nums <- 1:length(pvals)
  
  # Find which p-vals are less than bh cutoff
  under <- sp < (nums/(length(pvals))*alpha)
  
  # Return indices of significant p-vals
  if(sum(under) == 0){
    return(c())
  }else{
    cutoff <- max(which(under))
    return(ord[1:cutoff])
  }
}