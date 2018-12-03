#' fisher_trunc
#'
#' Performs a Fisher transformation on a number between -1 and 1, truncating to 
#' -1 or 1 if the input lies outside of that range.
#'
#' @param r Number between -1 and 1 to be transformed.
#'
#' @return Fisher transformed input truncated at -1 or 1 if outside the interval
#' (-1, 1).
#'
#' @export
fisher_trunc <- function(r){
  # Return large number (approx infinity) for r outside range
  if(r <= -1){
    return(-10000000)
  } else if(r >= 1){
    return(10000000)
  } else {
    # Otherwise, fisher transform as usual
    return(.5*log((1+r)/(1-r)))
  }
}

# Vectorize function.
fisher_trunc <- Vectorize(fisher_trunc)