#' run_DCM
#'
#' Searches for a differentially correlated set of variables in matrices M1 and
#' M2 using an iterative testing procedure.
#'
#' @param M1 p x n1 data matrix for group 1.
#' @param M2 p x n2 data matrix for group 2.
#' @param seed Vector giving initial set of variables.
#' @param del Vector containing variables to be ignored.
#' @param echo Boolean indicating whether information should be printed to the
#' console.
#' @param alpha Number between 0 and 1 indicating the level at which to control
#' the false discovery rate at each iteration.
#' @param max.iter Maximum number of iterations to be performed.
#' @param max.size Maximum size of set to consider.
#' @param validation Boolean indicating whether validation data should be stored
#' and returned.
#'
#' @return If \code{validation==TRUE}, list containing,
#' \itemize{
#'    \item found - Vector containing the variables in the differentially 
#'    correlated set found by the algorithm.
#'    \item mc1 - Mean correlation of \code{found} in group 1.
#'    \item mc2 - Mean correlation of \code{found} in group 2.
#'    \item its - Number of iterations performed.
#'    \item time - Time elapsed during procedure.
#'    \item pvals - Vector of p-values produced at last iteration of testing 
#'    procedure.
#'    \item startdels - Vector of variables which were ignored.
#'    \item it.sets - List of vectors containing differentially correlated sets
#'    at each iteration.
#'    \item it.pvals - List of vectors containing p-values produced at each
#'    iteration.
#'    \item it.test.stats - List of vectors containing test statistics produced
#'    at each iteration.
#'    \item it.std.errs - List of vectors containing standard errors produced at
#'    each iteration.
#'    \item it.diff.cor - Vector of mean differential correlations at each 
#'    iteration.
#'    \item reason.for.terminating - String describing reason why the algorithm
#'    terminated.
#'    \item flags - Vector of strings containing any warnings that arose during
#'    the search.
#' }
#' Else, list containing,
#' \itemize{
#'    \item found - Vector containing the variables in the differentially 
#'    correlated set found by the algorithm.
#'    \item mc1 - Mean correlation of \code{found} in group 1.
#'    \item mc2 - Mean correlation of \code{found} in group 2.
#'    \item its - Number of iterations performed.
#'    \item time - Time elapsed during procedure.
#'    \item pvals - Vector of p-values produced at last iteration of testing 
#'    procedure.
#'    \item startdels - Vector of variables which were ignored.
#'    \item reason.for.terminating - String describing reason why the algorithm
#'    terminated.
#' }
#'
#' @export
run_DCM <- function(M1, M2, seed, del = c(), echo = FALSE, alpha = 0.05, 
                    max.iter = 50, max.size = 1000, validation = FALSE){
  # Initialize validation data.
  it.sets <- list(seed)
  it.pvals <- list()
  it.diff.cor <- list()
  it.test.stats <- list()
  it.std.errs <- list()
  flags <- list()
  # Calculate runtime.
	starttime <- Sys.time()
	# Find dimensions.
	n1 <- ncol(M1) 
	n2 <- ncol(M2)
	p <- nrow(M1)
	idcs <- 1:p
	
	# If we are ignoring certain rows, remove them.
	if(length(del) > 0){
		idcs <- idcs[-del]
		M1 <- M1[-del,]
		M2 <- M2[-del,]
	}
	# Find new dimensions.
	p <- nrow(M1)
	P <- 1:p

	# Initialize relevant parameters.
	difference <- 1
	A <- which(idcs %in% seed)
	prevA <- p+1
	it <- 1
	startdels <- c() # Variables that are to be avoided
	
	if(echo){print("Beginning search...")}
	
	# Continue searching until convergence on non-degenerate set
	while(difference > 0 & length(A) > 5 & length(A) < max.size){
	  # Length of A.
	  k <- length(A)
	  
	  # Create separate data matrices.
	  xA1 <- M1[A,]
		xA2 <- M2[A,]
		y1 <- M1[-A,]
		y2 <- M2[-A,]
		
		## Calculate p-values for rows not in A.
		# Find mean vectors
		mean1 <- colMeans(xA1)
		mean2 <- colMeans(xA2)
		# Find norms
		n_m1 <- sqrt(sum(mean1^2))
		n_m2 <- sqrt(sum(mean2^2))
		# Find test quantities for all variables
		corsm1 <- t(cor(mean1, t(y1)))
		corsm2 <- t(cor(mean2, t(y2)))
		# Test stat and variance.
		obs <- corsm1*n_m1 - corsm2*n_m2
		sd <- sqrt(makeVars(y1, xA1) + makeVars(y2, xA2))
		# P-values.
		test_out <- pt(-obs/sd, min(c(n1-1, n2-1)), 0)
		
		## Calculate p-values for rows in A
		# Adjust means to not include row.
		mean1s <- -t(t(xA1) - mean1*k)/(k-1)
		mean2s <- -t(t(xA2) - mean2*k)/(k-1)
		# Find new mean norms.
		n_m1s <- sqrt(rowSums(mean1s*mean1s))
		n_m2s <- sqrt(rowSums(mean2s*mean2s))
		# Find cors of rows with means.
		corsm1 <- rowMeans(stdize(mean1s)*xA1)*n1
		corsm2 <- rowMeans(stdize(mean2s)*xA2)*n2
		# Test stat and variance.
		obss <- corsm1*n_m1s - corsm2*n_m2s
		sds <- sqrt(sapply(1:k, function(x){
		  makeVar(xA1[x,], xA1[-x,]) + makeVar(xA2[x,], xA2[-x,])
		  }))
		# P-values.
		test_in <- pt(-obss/sds, min(c(n1-1, n2-1)), 0)
		
		# Combine all p-values.
		test <- P
		test[-A] <- test_out
		test[A] <- test_in
		
		# Combine all test statistics.
		test.stats <- P
		test.stats[-A] <- obs
		test.stats[A] <- obss
		# Combine all standard errors.
		std.errs <- P
		std.errs[-A] <- sd
		std.errs[A] <- sds
		
		# Update A to include significant rows.
		newA <- bhy(test, alpha=alpha)
		
		# Check for convergence
		if(it >= max.iter){
			# Check for max iterations reached
			if(echo){print("Reached iteration limit.")}
		  reason.for.terminating <- "Reached iteration limit"
			difference <- 0
			newA <- newA[newA %in% A]
		} else {
			# Calculate difference between current and updated sets
			difference <- sum(!(newA %in% A)) + sum(!(A %in% newA))
			if(difference == 0){
			  reason.for.terminating <- "Reached fixed point"
			}
		}
		
		# Check for cycling.
		cycle.depth.to.check <- 5
		cycle.depth <- 0
		cycling <- FALSE
		A.sorted <- sort(newA)
	  for(cyc.ind in 2:min(cycle.depth.to.check, it)){
	    if(it > cyc.ind && identical(A.sorted, sort(it.sets[[it-cyc.ind+1]]))){
	      cycling <- TRUE
	      cycle.depth <- cyc.ind
	      if(echo){print(paste0("Cycling observed of order ", cyc.ind, " observed."))}
	      flags <- append(flags, paste0("Cycling observed of order ", cyc.ind, " observed."))
	      # If cycles are long-range, then just start over avoiding these 
	      # variables.
	      if(cycle.depth > 2){
	        startdels <- c(startdels, newA)
	        newA <- c()
	      }
	      break
	    }
	  }
		
		# Check to see if algorithm is oscillating between two similar sets or 
		# jumping sideways between sets.
		if(cycle.depth == 2){
			overlap <- newA[newA %in% A]
			# If oscillating, i.e. consecutive sets very similar, set new A to the
			# overlap between current A and previous A.
			if(min(length(setdiff(overlap, newA))/length(newA), 
			       length(setdiff(overlap, A))/length(A)) < .05){
			  if(echo){print("Oscillating observed. Setting new set to the intersection of current set and previous set.")}
			  flags <- append(flags, "Oscillating observed. Setting new set to the intersection of current set and previous set.")
			  reason.for.terminating <- "Oscillating between similar sets"
				newA <- overlap
				difference <- 0
			} else {
        if(echo){print("Consecutive sets not similar enough, starting a new search.")}
			  flags <- append(flags, "Consecutive sets not similar enough, starting a new search.")
				newA <- c()
			}
			startdels <- c(startdels, newA, prevA[!(prevA %in% newA)])
		}

		# Save old A, update current A
		prevA <- A
		A <- newA
		
		# Print progress if desired
		if(echo){print(sprintf("Size = %i", length(A)))}
		
		# Check to make sure that length of A is not too small or too large.
		if(length(A) <= 5){
		  reason.for.terminating <- "A became too small"
		  startdels <- c(startdels, A)
		} else if(length(A) >= max.size){
		  reason.for.terminating <- "A became too large"
		  startdels <- c(startdels, A)
		}
		
		# Store iteration data.
		it.pvals[[it]] <- test
		it.sets[[it+1]] <- newA
		it.diff.cor[[it]] <- mean(obss)
		it.test.stats[[it]] <- test.stats
		it.std.errs[[it]] <- std.errs
		
		# Count iterations
		it <- it + 1
	} #while(difference > 0 & length(A) > 10)
	print(reason.for.terminating)
	
	# New length of A
	k <- length(A)	

	# If algorithm didn't degenerate, find mean correlations.
	if(k > 1){
		mA1 <- apply(M1[A,], 2, mean)
		mA2 <- apply(M2[A,], 2, mean)
		
		meanA1 <- (sum(mA1^2)*(k^2) - k)/(k^2-k)
		meanA2 <- (sum(mA2^2)*(k^2) - k)/(k^2-k)
	} else {
		meanA1 <- NA
		meanA2 <- NA
	}
	
	# Runtime
	time <- as.double(difftime(Sys.time(), starttime, units = "secs"))

	# Save found set.
	found <- idcs[A]
		
	# Save converged set and properties
	if(validation){
	  return(list(found = found, 
	              mc1 = meanA1, 
	              mc2 = meanA2, 
	              its = it, 
	              time = time, 
	              pvals = test, 
	              startdels = startdels,
	              it.sets = it.sets,
	              it.pvals = it.pvals,
	              it.test.stats = it.test.stats,
	              it.std.errs = it.std.errs,
	              it.diff.cor = unlist(it.diff.cor),
	              reason.for.terminating = reason.for.terminating,
	              flags = unlist(flags)))
	} else {
	  return(list(found = found, 
	              mc1 = meanA1, 
	              mc2 = meanA2, 
	              its = it, 
	              time = time, 
	              pvals = test, 
	              startdels = startdels,
	              reason.for.terminating = reason.for.terminating))
	}
}
