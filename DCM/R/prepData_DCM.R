prepData_DCM <- function(mat.1, 
                         mat.2, 
                         QN              = FALSE, 
                         resid.full      = FALSE, 
                         echo            = FALSE, 
                         strict          = "low",
                         check.normality = FALSE){
  
  # Performs standard data preprocessing and checks necessary before running the
  #  DCM algorithm on input data matrices, mat.1 and mat.2.
  #
  # Args:
  #   -mat.1: Data matrix for group 1. Rows should correspond to variables and
  #     columns should correspond to samples.
  #   -mat.2: Data matrix for group 2. Orientation same as that of mat.1.
  #   -QN: If switched on, will perform quantile normalization.
  #   -resid.full: If switched on, will residualize full matrix.
  #   -echo: If switched on, will print information during computation.
  #   -strict: "low" or "high". If "high", will remove variables with variances
  #     <1/(100*n).
  # Returns:
  #   list containing updated data matrices for both groups as well as a vector 
  #     of indices corresponding to variables which made it through filtering.
  
	# Make sure data is numeric.
	mat.1 <- t(apply(mat.1, 1, as.numeric))
	mat.2 <- t(apply(mat.2, 1, as.numeric))
	
	# Rescale based on max and min.
	## Kevin: It seems this step is trying to rescale somehow? It doesn't seem to
	##  disturb the normality of the data but it does seem to shift if to be 
	##  centered at 1. Unclear if this is what we want. Maybe we should do
	##  (x - min)/(max-min) instead.
	mat.1 <- (mat.1 - min(mat.1, na.rm = TRUE))/max(mat.1, na.rm = TRUE)
	mat.2 <- (mat.2 - min(mat.2, na.rm = TRUE))/max(mat.2, na.rm = TRUE)
	
	# Round.	
	mat.1 <- round(mat.1, 5)
	mat.2 <- round(mat.2, 5)
	
	# Make sure same number of rows
	p1 <- nrow(mat.1)
	p2 <- nrow(mat.2)
	if(p1 != p2){
	  print("Datasets must have the same number of rows.")
		break
	}
	p <- p1
	
	# Make sure 3 or more samples.
	n1 <- ncol(mat.1)
	n2 <- ncol(mat.2)
	if(n1 < 3 | n2 < 3){
		print("Must have at least 3 samples in each group.")
		break
	}
			
	# Remove and report bad rows
	bad.rows <- sanitize_DCM(mat.1, mat.2, strict)
	if(length(bad.rows) > 0){
		if(echo){print(sprintf("Removing %i rows - too many zero expressions.", 
		                       length(bad.rows)))}
		
		# Record dimension
		good.idcs <- (1:p)[-bad.rows]

		# Remove ignored rows
		mat.1 <- mat.1[-bad.rows,]
		mat.2 <- mat.2[-bad.rows,]	
	}else{
		good.idcs <- 1:p
	}
		
	p <- length(good.idcs)
	
	# If not quantile normalizing, check normality and warn if non-Gaussian
	if(!QN & check.normality){
		
		rand = sample((1:p), min(1000, p))
		pvs1 = sapply(rand, function(x) shapiro.test(mat.1[x,])$p.value)
		pvs2 = sapply(rand, function(x) shapiro.test(mat.2[x,])$p.value)
		
		bad1 = length(bh(pvs1)) > 0
		bad2 = length(bh(pvs2)) > 0
		
		if((bad1 | bad2) & echo){
			print("Warning: Data is non-normal; consider quantile normalizing.  (Set QN = TRUE)")
		}

	}
		
	print("Done checking data.")
	
	################# Process data ###################
	
	# Quantile normalize
	if(QN){
		mat.1 <- quantNorm(mat.1)
		mat.2 <- quantNorm(mat.2)
	}
	
	# Standardize (for easier calculation of correlations)
	## Kevin: Not sure about the standardization that is happening here. After,
	##  rows have neither sd=1 nor sum of squares=1.
	#mat.1 <- stdize(mat.1)
	#mat.2 <- stdize(mat.2)
	mat.1 <- stdize_KO(mat.1)
	mat.2 <- stdize_KO(mat.2)
	
	col.mean.1 <- apply(mat.1, 2, mean)
	col.mean.2 <- apply(mat.2, 2, mean)
	
	# Compute overall correlations for each group.
	c1 <- (sum(col.mean.1^2)*(n1^2) - n1)/(n1^2-n1)
	c2 <- (sum(col.mean.2^2)*(n2^2) - n2)/(n2^2-n2)
	if(echo){
		print(sprintf("Overall cor, group 1: %2f", c1))
		print(sprintf("Overall cor, group 2: %2f", c2))
	}
	
	# Residualize full matrix, if worried about systematic correlation within 
	#  groups.
	if(resid.full){
		mat.1 <- resid_DCM(mat.1, QN = QN)
		mat.2 <- resid_DCM(mat.2, QN = QN)
	} else if(abs(c1 - c2) > .1){
		if(echo){print("Warning: Large overall correlation difference.  Consider 
		               residualizing full matrix first.")}
	}
	
	if(echo){print("Done preparing data.")}

	return(list(mat.1 = mat.1, mat.2 = mat.2, idcs = good.idcs))
}