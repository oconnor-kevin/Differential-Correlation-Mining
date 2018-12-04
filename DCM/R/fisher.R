#' fisher
#'
#' Performs a Fisher transformation on a number between -1 and 1. If input is 
#' outside that range, returns 0.
#'
#' @param r Number between -1 and 1 to be transformed.
#'
#' @return If \code{-1 < r < 1}, Fisher transformed input. Else, 0.
#'
#' @export
fisher <- function(r){
	# Set r=0 if outside range.
	r <- r*(abs(r) < 1)
	
	# Return fisher transform
	return(.5*log((1+r)/(1-r)))
}