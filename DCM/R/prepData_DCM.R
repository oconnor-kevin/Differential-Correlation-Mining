#' prepData_DCM
#'
#' Takes two matrices and prepares them for use in the DCM procedure.
#'
#' @param M1 p x n1 data matrix for group 1.
#' @param M2 p x n2 data matrix for group 2.
#' @param QN Boolean indicating whether quantile normalization should be
#' performed.
#' @param resid.full Boolean indicating whether residualization of the full
#' matrix should be performed.
#' @param echo Boolean indicating whether information should be printed to the
#' console.
#' @param strict String indicating whether rows with low variance should be
#' filtered from the matrices. See sanitize_DCM.R for usage. (\code{"high"} or
#' \code{"low"}).
#' @param validation Boolean indicating whether validation data should be stored
#' and returned.
#'
#' @return If \code{validation == TRUE}, list containing,
#' \itemize{
#'    \item M1 - Prepared data matrix for group 1.
#'    \item M2 - Prepared data matrix for group 2.
#'    \item idcs - Vector of rows which have been included.
#'    \item preproc.steps - Vector containing strings that describe preprocessing
#'    steps that were performed.
#'    \item data.flags - Vector containing strings that describe any flags about
#'    the data that were observed.
#' }
#' Else, list containing,
#' \itemize{
#'    \item M1 - Prepared data matrix for group 1.
#'    \item M2 - Prepared data matrix for group 2.
#'    \item idcs - Vector of rows which have been included.
#' }
#'
#' @export
prepData_DCM <- function(M1, M2, QN = FALSE, resid.full = FALSE, echo = FALSE, 
                         strict = "low", validation = FALSE){
  # Initialize validation datasets.
  preproc.steps <- list()
  data.flags <- list()
  
	# Make sure data is numeric
	M1 <- t(apply(M1, 1, as.numeric))
	M2 <- t(apply(M2, 1, as.numeric))
	preproc.steps <- append(preproc.steps, "Cast data as numeric.")
	# Rescale.
	M1 <- (M1 - min(M1, na.rm = TRUE))/max(M1, na.rm = TRUE)
	M2 <- (M2 - min(M2, na.rm = TRUE))/max(M2, na.rm = TRUE)
	preproc.steps <- append(preproc.steps, "Subtract min and divide by max excluding na's.")	
	# Round to 5 decimal places.
	M1 <- round(M1, 5)
	M2 <- round(M2, 5)
	preproc.steps <- append(preproc.steps, "Round to 5 decimal places.")
	
	# Make sure same number of rows
	p1 <- nrow(M1)
	p2 <- nrow(M2)
	if(p1 != p2){
		stop("Datasets must have the same number of rows.")
	} else {
		p <- p1
	}
	preproc.steps <- append(preproc.steps, "Confirmed that datasets have same number of rows.")
	
	# Make sure 3 or more samples
	n1 <- ncol(M1)
	n2 <- ncol(M2)
	if(n1 < 3 | n2 < 3){
	  stop("Must have at least 3 samples in each group.")
	}
	preproc.steps <- append(preproc.steps, "Confirmed that both datasets have 3 or more samples.")	
			
	# Remove and report bad rows
	sanitize_result <- sanitize_DCM(M1, M2, strict, validation)
	if (validation){
	  killrows <- sanitize_result$killrows
	  na1 <- which(sanitize_result$na1)
	  na2 <- which(sanitize_result$na2)
	  lb1 <- which(sanitize_result$lb1)
	  lb2 <- which(sanitize_result$lb2)
	  ub1 <- which(sanitize_result$ub1)
	  ub2 <- which(sanitize_result$ub2)
	  var1 <- which(sanitize_result$var1)
	  var2 <- which(sanitize_result$var2)
	} else {
	  killrows <- sanitize_result
	}
	if(length(killrows) > 0){
		if(echo){print(sprintf("Removing %i rows - too many zero expressions.", length(killrows)))}
		
		# Record dimension
		idcs = (1:p)[-killrows]

		# Remove ignored rows
		M1 = M1[-killrows,]
		M2 = M2[-killrows,]
		
		# Store preprocessing steps.
		preproc.steps <- append(preproc.steps, "Checked for rows with na values.")
		preproc.steps <- append(preproc.steps, "Checked for rows with too many values equal to the smallest value.")
		preproc.steps <- append(preproc.steps, "Checked for rows with too many values equal to the largest value.")
		if(strict == "high"){
		  preproc.steps <- append(prepoc.steps, "Checked for rows with variance < 1/(100*n).")
		}
		
		# Store data flags.
		if(validation){
  		if(length(na1) > 0){
  		  data.flags <- append(data.flags, 
  		                       paste0("Found na values in the following rows of group 1: ", 
  		                              paste(na1, collapse=" ")))
  		}
  		if(length(na2) > 0){
  		  data.flags <- append(data.flags, 
  		                       paste0("Found na values in the following rows of group 2: ", 
  		                              paste(na2, collapse=" ")))
  		}
  		if(length(lb1) > 0){
  		  data.flags <- append(data.flags, 
  		                       paste0("Found too many small values in the following rows of group 1: ", 
  		                              paste(lb1, collapse=" ")))
  		}
  		if(length(lb2) > 0){
  		  data.flags <- append(data.flags, 
  		                       paste0("Found too many small values in the following rows of group 2: ", 
  		                              paste(lb2, collapse=" ")))
  		}
  		if(length(ub1) > 0){
  		  data.flags <- append(data.flags, 
  		                       paste0("Found too many large values in the following rows of group 1: ", 
  		                              paste(ub1, collapse=" ")))
  		}
  		if(length(ub2) > 0){
  		  data.flags <- append(data.flags, 
  		                       paste0("Found too many large values in the following rows of group 2: ", 
  		                              paste(ub2, collapse=" ")))
  		}
  		if(length(var1) > 0 & strict == "high"){
  		  data.flags <- append(data.flags, 
  		                       paste0("Found variance < 1/(100*n) in the following rows of group 1: ", 
  		                              paste(var1, collapse=" ")))
  		}
  		if(length(var2) > 0 & strict == "high"){
  		  data.flags <- append(data.flags, 
  		                       paste0("Found variance < 1/(100*n) in the following rows of group 2: ", 
  		                              paste(var2, collapse=" ")))
  		}
		}
	} else {
		idcs <- 1:p
		preproc.steps <- append(preproc.steps, "Checked for rows with too many values equal to smallest value and found none.")
		preproc.steps <- append(preproc.steps, "Checked for rows with too many values equal to largest value and found none.")
		if (strict == "high"){
		  preproc.steps <- append(preproc.steps, "Checked for rows with variance < 1/(100*n) and found none.")  
		}
	}
		
	p <- length(idcs)
	
	# If not quantile normalizing, check normality and warn if non-Gaussian
	if(!QN){
		rand <- sample((1:p), min(1000, p))
		pvs1 <- sapply(rand, function(x){shapiro.test(M1[x,])$p.value})
		pvs2 <- sapply(rand, function(x){shapiro.test(M2[x,])$p.value})
		
		bad1 <- length(bh(pvs1)) > 0
		bad2 <- length(bh(pvs2)) > 0
		
		preproc.steps <- append(preproc.steps, "Checked normality of data.")
		
		if((bad1 | bad2) & echo){
			print("Warning: Data is non-normal; consider quantile normalizing.  (Set QN = TRUE)")
		  data.flags <- append(data.flags, "Data is not normal. Consider quantile normalizing by setting QN=TRUE.")
		}
	}
		
	print("Done checking data.")
	
	################# Process data ###################
	
	# Quantile normalize
	if(QN){
		M1 <- quantNorm(M1)
		M2 <- quantNorm(M2)
		preproc.steps <- append(preproc.steps, "Quantile normalized each group (see quantNorm.R).")
	}
	
	# Standardize (for easier calculation of correlations).
	M1 <- stdize(M1)
	M2 <- stdize(M2)
	preproc.steps <- append(preproc.steps, "Row-standardized each group (see stdize.R).")
	
	m1 <- apply(M1, 2, mean)
	m2 <- apply(M2, 2, mean)
	c1 <- (sum(m1^2)*(n1^2) - n1)/(n1^2-n1)
	c2 <- (sum(m2^2)*(n2^2) - n2)/(n2^2-n2)
	
	if(echo){
		print(sprintf("Overall cor, group 1: %2f", c1))
		print(sprintf("Overall cor, group 2: %2f", c2))
	}
	
	# Residualize full matrix, if worried about systematic correlation within groups
	if(resid.full){
		M1 <- resid_DCM(M1, QN = QN)
		M2 <- resid_DCM(M2, QN = QN)
		preproc.steps <- append(preproc.steps, "Residualized full matrix (see resid_DCM.R).")
	}else if(abs(c1 - c2) > .1){
		if(echo){print("Warning: Large overall correlation difference.  Consider residualizing full matrix first.")}
	  preproc.steps <- append(preproc.steps, "Checked that difference in correlation between the two groups is not too big.")
	  data.flags <- append(data.flags, "Large overall correlation difference. Consider residualizing full matrix first (resid.full=TRUE).")
	}
	
	if(echo){print("Done preparing data.")}

	if(validation){
	  return(list("M1" = M1, 
	              "M2" = M2, 
	              "idcs" = idcs,
	              "preproc.steps" = unlist(preproc.steps),
	              "data.flags" = unlist(data.flags)))
	} else {
	  return(list("M1" = M1, 
	              "M2" = M2, 
	              "idcs" = idcs))
	}
}