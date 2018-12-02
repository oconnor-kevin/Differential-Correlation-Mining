#' init_DCM
#'
#' Takes prepared Matrix (row-normalized and optionally Quantile Normalized).
#'  Looks for high correlation in M1 and low/negative correlation in M2, finds 
#'  group of size k to be used as initial set in iterative testing procedure.
#'
#' @param M1 p x n1 data matrix for group 1.
#' @param M2 p x n2 data matrix for group 2.
#' @param k Desired size of initial set.
#' @param start Desired initial set.
#' @param del Vector of variables to be excluded from initialization set.
#'
#' @return List containing,
#' \itemize{
#'    \item seed - Vector containing the seed set used.
#'    \item found - Vector containing the terminating set to be used as initial 
#'    set for iterative testing procedure.
#'    \item iterations - Number of iterations performed.
#'    \item time - Amount of time that the initialization procedure took.
#'    \item init.sets - List containing set at each iteration.
#'    \item it.scores - Vector containing the score for the set at each 
#'    iteration.
#'    \item init.steps - List containing strings that give information about 
#'    what steps are performed during initialization.
#' }
#'
#' @export
init_DCM <- function(M1, M2, k, start = c(), del = c()){
  # Initialize validation data.
  init.sets <- list()
  it.scores <- list()
  init.steps <- list()
  
  # Make sure we don't delete too many variables.
  if(length(del) + k > nrow(M1)){
    stop("Number of variables to delete plus set size larger than number of variables.")
  }
  
  # Remove unwanted variables.
  if(length(del) > 0){
    # Record dimension.
    real_p <- dim(M1)[1]
    idcs <- (1:real_p)[-del]
    
    # Remove ignored rows.
    M1 <- M1[-del,]
    M2 <- M2[-del,]
    
    init.steps <- paste0("Removed the following variables before initializing: ", paste(del, collapse=" "))
	} else {
		# Record dimension.
		p <- dim(M1)[1]
		idcs <- 1:p
	}
    
  # Lengths.
  n1 <- dim(M1)[2]
  n2 <- dim(M2)[2]
  p  <- dim(M1)[1]

	# Start with prespecified row set or with a random k rows.
	if(length(start) == k){
		A <- start
		init.steps <- append(init.steps, paste0("Starting from user-specified initial set: ", paste(start, collapse=" ")))
	} else {
		A <- sample(length(idcs), k)
		if(length(start) > 0){
		  init.steps <- append(init.steps, paste0("Starting set supplied by user does not match size specified by user, starting from random set of size ", k, "."))		  
		} else {
		  init.steps <- append(init.steps, paste0("Starting from random set of size ", k, "."))		  
		}
	}
	
	# Save starting point.
	orig_A <- idcs[A]

	# Make list of all indices, and of those not in seed set.
	P <- 1:p
	notA <- P[!(P %in% A)]

	# Find correlations of all genes with only genes in A. Assumes that M1 and M2
	#   are row-normalized. Should be taken care of by prepData_DCM.R.
	cross_1 <- round(M1 %*% t(M1[A,]), digits = 10)
	cross_2 <- round(M2 %*% t(M2[A,]), digits = 10)
	# Set 1's to 0 so that transformed value will be 0 and variance won't 
	#   contribute to sum
	cross_1[cross_1 == 1] <- 0
	cross_2[cross_2 == 1] <- 0
	# Fisher transform the correlations.
	cross_1 <- fisher(cross_1)*sqrt(n1 - 3)
	cross_2 <- fisher(cross_2)*sqrt(n2 - 3)
	# Note: resulting matrices have rows in order of data indices, columns 
	#   correspond to values of A in the order listed by A.

	# Rowsums represent total of pairwise correlations between [row] and A.
	rows_1 <- rowSums(cross_1)
	rows_2 <- rowSums(cross_2)

	# Initialize parameters for loop.
	done <- FALSE
	it <- 0
	d <- (p-k)*k # Number of in/out swap combos
	start.time <- Sys.time()
	
	# Iterate until local max of score is reached.
	while(!done){
    # Rowsum Differences.
    diffs12 <- rows_1 - rows_2
		
    # Save initial score.
    it.scores <- append(it.scores, mean(diffs12[A]))
    
		maxA <- function(i){
		  # Compute the change in score from swapping each variable outside of A with 
		  #   i'th variable in A. Return vector with index that maximizes the increase
		  #   and the corresponding increase in score.
			temp <- diffs12[notA] + c(cross_2[notA, i]) - c(cross_1[notA, i])
			idx <- which.max(temp)
			effect <- temp[idx]
			return(c(idx, effect))
		}

		# Find optimal swap of a gene in A for a gene outside of A.
		bestAs <- t(sapply(1:k, function(x){maxA(x)}))
		bestAs[,2] <- bestAs[,2] - diffs12[A]
		best_out <- which.max(bestAs[,2])
		best_in <- bestAs[best_out,1]
		
		if(bestAs[best_out,2] <= 0){ # If optimal swap decreases score
			done <- TRUE
		} else {
			# Find in and out data labels.
			out <- A[best_out]
			inn <- notA[best_in]

			# Switch "out" and "in" indices from their lists.
			A[best_out] <- inn
			notA[best_in] <- out
		
			# Find correlation vector for "inn".
			new_1 <- round(M1 %*% M1[inn,], digits = 10)
			new_1[inn] <- 0 # So Fisher isn't inf
			new_1 <- fisher(new_1)*sqrt(n1 - 3)
			# Edit rowsums to reflect inclusion of new index, exclusion of old.
			rows_1 <- rows_1 - cross_1[,best_out] + new_1
			# Replace relevant column of corr matrix with new correlations.
			cross_1[,best_out] <- new_1
		
			# Same thing for Group 2
			new_2 <- round(M2 %*% M2[inn,], digits = 10)
			new_2[inn] <- 0
			new_2 <- fisher(new_2)*sqrt(n2 - 3)
			rows_2 <- rows_2 - cross_2[,best_out] + new_2
			cross_2[,best_out] <- new_2
			
			# Update iteration count and store set.
			it <- it + 1
			init.sets[[it]] <- idcs[A]
		}
	}
	time <- difftime(Sys.time(), start.time, units="secs")
	
	# Translate back to real indices
	A <- idcs[A]
	
	init.steps <- append(init.steps, paste0("Performed ", it, " iterations during initialization."))
	
	return(list(seed=orig_A,
	            found=A, 
	            iterations=it, 
	            time=time,
	            init.sets=init.sets,
	            it.scores=unlist(it.scores),
	            init.steps=init.steps))	
}
