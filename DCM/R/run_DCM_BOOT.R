# TITLE: run_DCM_BOOT.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 9/15/18
# DATE MODIFIED: 9/15/18

run_DCM_BOOT <- function(mat.1, 
                         mat.2, 
                         seed,
                         echo = FALSE, 
                         alpha = 0.05, 
                         max.iter = 10,
                         boot.iter = 1000){
	
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
	#M1A.cor <- mat.1[idcs[A]] %*% t(mat.1[idcs[A]])
	#M2A.cor <- mat.2[idcs[A]] %*% t(mat.2[idcs[A]])
	#it_diff_cor[[1]] <- mean(M1A.cor-M2A.cor)
	
	# Continue searching until convergence on non-degenerate set
	while(difference > 0 & length(A) > 5){
		# Start iteration timer.
	  it_start <- Sys.time()
	  
	  # Store separate data matrices.
	  xA1 <- mat.1[A,]
	  xA2 <- mat.2[A,]
	  
	  y1 <- mat.1[-A,]
	  y2 <- mat.2[-A,]
	  
	  # Compute observed test statistics.
	  obs.test.stats <- 1:p
	  obs.test.stats[-A] <- colMeans(xA1 %*% t(y1) - xA2 %*% t(y2))
	  obs.test.stats[A] <- colMeans(xA1 %*% t(xA1) - xA2 %*% t(xA2))
	  
    start.time <- Sys.time()
    # Get bootstrapped distribution of test statistics.
    excedences <- rep(0, p)
	  for(b in 1:boot.iter){
	    # Shuffle observations between groups.
	    samp <- sample.int(n1+n2)
	    xA <- cbind(xA1, xA2)
	    A1.samp <- xA[, samp[1:n1]]
	    A2.samp <- xA[, samp[(n1+1):(n1+n2)]]
	     
	    # Restandardize to mean=0, ss=1.
	    A1.samp <- stdize_KO(A1.samp)
	    A2.samp <- stdize_KO(A2.samp)
	    
	    # Compute test statistics for genes outsides of A.
	    test.stats <- 1:p
	    test.stats[-A] <- colMeans(A1.samp %*% t(y1) - A2.samp %*% t(y2))
	    
	    # Compute test statistics for genes in A. This is biased but much faster.
	    test.stats[A] <- colMeans(A1.samp %*% t(A1.samp) - A2.samp %*% t(A2.samp))
	    
	    ## TOO SLOW
	    # Compute test statistics for genes in A.
	    #for(v in 1:length(A)){
	    #   test.stats[A[v]] <- mean(A1.samp[-v,] %*% A1.samp[v,] - A2.samp[-v,] %*% A2.samp[v,])
	    #}
	    
	    # Check whether each test stats is greater than the observed test stat.
	    excedences <- excedences + as.numeric(obs.test.stats < test.stats)	    
	  }
	  
	  # Compute p-values.
	  #p.vals <- sapply(1:p, function(i){mean(obs.test.stats[i] < test.stats.mat[i,])})
    p.vals <- excedences/boot.iter
	  
  	# Update A to include significant rows
  	newA <- bhy(p.vals, alpha = alpha)
  	
  	Sys.time() - start.time
  
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
  	it_p_vals[[it]] <- p.vals
  	#it_p_vals[[it]] <- test
	  #it_test_stats[[it]] <- c(obs, obss)
	  #it_test_var[[it]] <- c(sd^2, sds^2)
	  
	  # Store mean differential correlation.
	  #M1A.cor <- mat.1[idcs[A]] %*% t(mat.1[idcs[A]])
	  #M2A.cor <- mat.2[idcs[A]] %*% t(mat.2[idcs[A]])
	  #it_diff_cor[[it+1]] <- mean(M1A.cor-M2A.cor)
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
	            pvals = p.vals, 
	            startdels = startdels, 
	            it_sets = it_sets, 
	            it_p_vals = it_p_vals, 
	            it_test_stats = it_test_stats,
	            it_test_var = it_test_var,
	            it_times = it_times,
	            it_diff_cor = it_diff_cor))
}
