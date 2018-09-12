# TITLE: run_DCM_KO.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 9/6/18
# DATE MODIFIED: 9/6/18

run_DCM_KO <- function(mat.1, 
                       mat.2, 
                       seed,
                       echo = FALSE, 
                       alpha = 0.05, 
                       max.iter = 50){
	
	# Calculate runtime
  starttime = Sys.time()
	
	# Initialize data lists.
	it_times      <- list()
	it_test_stats <- list()
	it_test_var   <- list()
	it_p_vals     <- list()
	it_sets       <- list(seed)
	it_diff_cor   <- list()

	# Find dimensions.	
	n1   <- ncol(mat.1) 
	n2   <- ncol(mat.2)
	p    <- nrow(mat.1)
	idcs <- 1:p

	# Initialize iteration variables.
	difference <- 1
	A          <- which(idcs %in% seed)
	prevA      <- p+1
	it         <- 0
	osc        <- FALSE
	diags      <- FALSE
	startdels  <- c()
	
	# Initialize list of differential correlations.
	# Store mean differential correlation.
	M1A.cor <- M1[idcs[A]] %*% t(M1[idcs[A]])
	M2A.cor <- M2[idcs[A]] %*% t(M2[idcs[A]])
	it_diff_cor[[1]] <- mean(M1A.cor-M2A.cor)
	
	# Continue searching until convergence on non-degenerate set
	while(difference > 0 & length(A) > 5){
		# Start iteration timer.
	  it_start <- Sys.time()
	  
	  # Store separate data matrices.
	  xA1 <- mat.1[A,]
	  xA2 <- mat.2[A,]
	  
	  y1 <- mat.1[-A,]
	  y2 <- mat.2[-A,]
	  
    # Calculate standard errors.
    #sd <- sqrt(sapply(1:nrow(y1), function(var){sqrt(makeVar_KO(y1[var,], xA1)^2/n1 + makeVar_KO(y2[var,], xA2)^2/n2)}))
    #sds <- sqrt(sapply(1:nrow(xA1), function(var){sqrt(makeVar_KO(xA1[var,], xA1[-var,])^2/n1 
    #                                                   + makeVar_KO(xA2[var,], xA2[-var,])^2/n2)}))
	  sd <- sqrt(makeVars(y1, xA1)^2/n1 + makeVars(y2, xA2)^2/n2)
	  sds <- sqrt(sapply(1:nrow(xA1), function(x){makeVar(xA1[x,], xA1[-x,])^2/n1 + makeVar(xA2[x,], xA2[-x,])^2/n2}))
    
    std.errs <- 1:p
    std.errs[-A] <- sd
    std.errs[A] <- sds
    
    ## Calculate p-values for rows in outside of A.
    # Find mean vectors
    mean1 <- colMeans(xA1)
    mean2 <- colMeans(xA2)
    
    # Find norms
    n_m1 <- sqrt(sum(mean1^2))
    n_m2 <- sqrt(sum(mean2^2))
    
    # Length of A
    k <- length(A)
    
    # Find test quantities for all variables
    corsm1 <- t(cor(mean1, t(y1)))
    corsm2 <- t(cor(mean2, t(y2)))
    
    # Test stat and variance
    obs <- corsm1*n_m1 - corsm2*n_m2
    
    # p-values for rows not in A
    test_out <- pt(-obs/sd, min(c(n1-1, n2-1)), 0)
    
    ## Calculate p-values for rows in A
    # Adjust means to not include row
    mean1s <- -t(t(xA1) - mean1*k)/(k-1)
    mean2s <- -t(t(xA2) - mean2*k)/(k-1)
    
    # Find new mean norms
    n_m1s <- sqrt(rowSums(mean1s*mean1s))
    n_m2s <- sqrt(rowSums(mean2s*mean2s))
    
    # Find cors of rows with means
    corsm1 <- rowMeans(stdize_KO(mean1s)*xA1)*n1
    corsm2 <- rowMeans(stdize_KO(mean2s)*xA2)*n2    
    
    # Compute test statistic.
    obss <- corsm1*n_m1s - corsm2*n_m2s
    
    # Find pvals
    test_in <- pt(-obss/sds, min(c(n1-1, n2-1)), 0)
    
    # Combine all p-values
    test = idcs
    test[-A] = test_out
    test[A] = test_in
  	
  	# Update A to include significant rows
  	newA <- bhy(test, alpha = alpha)
  
  	# Check for convergence
  	if(it >= max.iter){
  		# Check for max iterations reached
  		if(echo){print("Reached iteration limit.")}
  		difference = 0
  		newA = newA[newA %in% A]
  	}else{
  		# Calculate difference between current and updated sets
  		difference = sum(!(newA %in% A)) + sum(!(A %in% newA))
  	}
  		
  	print(newA)

  	# Check to see if algorithm is oscillating between two similar sets or jumping sideways between sets
  	if(sum(!(newA %in% prevA)) + sum(!(prevA %in% newA)) == 0){
  		if(echo){print("Oscillating...")}
  			
  		# If already oscillated, don't keep going.
  		if(osc){
  				
  			difference = 0
  			overlap = newA[newA %in% A]
  			
  			if(min(length(setdiff(overlap, newA))/length(newA), length(setdiff(overlap, A))/length(A)) < .05){
  				newA = overlap
  			}else{
  				newA = c()
  			}
  			startdels = c(startdels, newA, prevA[!(prevA %in% newA)])
  			osc = FALSE
  				
  		}else{
  			
  			osc = TRUE
  			newA = A[A %in% newA]
  				
  		}
  	}
  		
  	# Save old A, update current A
  	prevA = A
  	A = newA
  		
  	#print(prevA)
  	#print(newA)
  		
  	# Check for extended cycling. Assumes each set in it_sets is sorted.
  	A_sorted <- sort(A)
  	if(any(unlist(lapply(it_sets, function(s){
  	  if(length(s) == length(A_sorted)){
  	    return(all(s==A_sorted))
  	  } else {
  	    return(FALSE)
  	  }})))){
  	    if(echo){print("Cycle Detected")}
  	    difference <- 0
	  }		
  		
  	# Print progress if desired
  	if(echo){print(sprintf("Size = %i", length(A)))}	
  	
  	# Count iterations
  	it = it + 1
  		
  	# Stop iteration timer and save.
  	it_times[[it]] <- difftime(Sys.time(), it_start, units="secs")
  		
  	# Store iteration data.
  	it_sets[[it+1]] <- A_sorted
  	it_p_vals[[it]] <- test
	  it_test_stats[[it]] <- c(obs, obss)
	  it_test_var[[it]] <- c(sd^2, sds^2)
	  
	  # Store mean differential correlation.
	  M1A.cor <- M1[idcs[A]] %*% t(M1[idcs[A]])
	  M2A.cor <- M2[idcs[A]] %*% t(M2[idcs[A]])
	  it_diff_cor[[it+1]] <- mean(M1A.cor-M2A.cor)
  } #while(difference > 0 & length(A) > 10)
	
	# New length of A
	k = length(A)	

	# If algorithm didn't degenerate, find mean correlations
	if(k > 1){
		mA1 = apply(mat.1[A,], 2, mean)
		mA2 = apply(mat.2[A,], 2, mean)
		
		meanA1 = (sum(mA1^2)*(k^2) - k)/(k^2-k)
		meanA2 = (sum(mA2^2)*(k^2) - k)/(k^2-k)
  }else{
		meanA1 = NA
		meanA2 = NA
	}
	
	# Runtime
	time = as.double(difftime(Sys.time(), starttime, units = "secs"))
	
	# Deal with diag sitch
	if(diags){
		found = c(0, cut, idcs[A])
	}else{
		found = idcs[A]
	}
		
	# Save converged set and relevant data.
	return(list(found = found, 
	            mc1 = meanA1, 
	            mc2 = meanA2, 
	            its = it, 
	            time = time, 
	            pvals = test, 
	            startdels = startdels, 
	            it_sets = it_sets, 
	            it_p_vals = it_p_vals, 
	            it_test_stats = it_test_stats,
	            it_test_var = it_test_var,
	            it_times = it_times,
	            it_diff_cor = it_diff_cor))
}
