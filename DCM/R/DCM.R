#' DCM
#' 
#' Given two data matrices, runs the DCM algorithm, searching for sets of
#' variables with high correlation in the first group and no/less correlation in 
#' the second group.
#' 
#' @param M1 p x n1 data matrix for first group.
#' @param M2 p x n2 data matrix for second group.
#' @param max.groups Maximum number of DC sets to search for.
#' @param max.iter Maximum number of iterations to perform at each search.
#' @param max.time Maximum number of hours to run algorithm for.
#' @param est.size Integer giving size of initial set.
#' @param max.size Integer giving maximum set size to consider.
#' @param alpha Number between 0 and 1 giving level at which to control false
#' discovery rate at each iteration.
#' @param start Vector of variables to use as initial set.
#' @param QN Boolean indicating whether data matrices should be quantile 
#' normalized.
#' @param resid.full Boolean indicating whether full data matrices should be
#' residualized. Set to true if worried about systematic correlation differences
#' between the two groups.
#' @param echo Boolean indicating whether information should be printed during
#' runtime.
#' @param strict String indicating whether variables with low variances should 
#' be removed (\code{"high"} or \code{"low"}).
#' @param validation Boolean indicating whether validation data should be 
#' collected and compiled into pdf and txt documents.
#' @param validation.dir Directory in which validation documents should be 
#' saved.
#' 
#' @return If \code{validation==TRUE}, list containing
#' \itemize{
#'    \item DC_sets - List of vectors containing differentially correlated sets
#'    found from each search.
#'    \item iterations - List containing number of iterations from each search.
#'    \item time - List containing duration in seconds of each search.
#'    \item meanCor1 - List containing mean correlation in group 1 of the found
#'    set.
#'    \item meanCor2 - List containing mean correlation in group 2 of the found
#'    set.
#'    \item it.pvals - List of lists containing p-values from each iteration of 
#'    each search.
#'    \item it.test.stats - List of lists containing test statistics from each
#'    iteration of each search.
#'    \item it.std.errs - List of lists containing standard errors from each 
#'    iteration of each search.
#'    \item it.sets - List of lists containing sets from each iteration of each
#'    search.
#'    \item inits - List of lists containing initialization data from each 
#'    search.
#'    \item indices - Vector of indices indicating rows of the data that were 
#'    used.
#' }
#' Else, list containing
#' \itemize{
#'    \item DC_sets - List of vectors containing differentially correlated sets
#'    found from each search.
#'    \item iterations - List containing number of iterations from each search.
#'    \item time - List containing duration in seconds of each search.
#'    \item meanCor1 - List containing mean correlation in group 1 of the found
#'    set.
#'    \item meanCor2 - List containing mean correlation in group 2 of the found
#'    set.
#'    \item indices - Vector of indices indicating rows of the data that were 
#'    used.
#' }
#' 
#' @export
DCM <- function(M1, M2, max.groups = 5, max.iter = 50, max.time = 1, 
                est.size = nrow(M1)/10, max.size = 1000, alpha = 0.05, 
                start = c(), QN = FALSE, resid.full = FALSE, echo = FALSE, 
                strict = "low", validation = FALSE, validation.dir = ''){
  # Make sure validation directory is valid if validation==TRUE.
  if(validation & validation.dir == ""){
    validation.dir <- getwd()
  }
  
	################### Prepare Data #############################################
	prep <- prepData_DCM(M1, M2, QN, resid.full, echo, strict, validation)
	M1 <- prep$M1
	M2 <- prep$M2
	idcs <- prep$idcs
	if(validation){
	  preproc.steps <- prep$preproc.steps
	  data.flags <- prep$data.flags
	}
	p <- nrow(M1)
	################### Run Algorithm ############################################
	# Prepare to store results
	DC_sets <- list()
	iterations <- list()
	time <- list()
	meanCor1 <- list()
	meanCor2 <- list()
	it.pvals <- list()
	it.test.stats <- list()
	it.std.errs <- list()
	it.sets <- list()
	it.diff.cor <- list()
	inits <- list()

	# Initialize parameters.
	starttime <- Sys.time()
	tottime <- 0
	n.sets <- 1
	startdels <- c() # Vars to delete.
	
	# Keep searching until all rows or exhausted, or max groups or time is reached
	while(length(startdels) < p-est.size & n.sets <= max.groups & tottime < max.time){
	  # Find initial starting set, ignoring starting points already tried.
		if(length(start) < 5){
		  if(validation & length(start) > 0){
		    data.flags <- append(data.flags, "Desired starting set too small (< 5) so it is being ignored.")
		  }
			inits[[n.sets]] <- init_DCM(M1, M2, k=est.size, del=startdels)
			tmp <- inits[[n.sets]]$found
		}else{
		  inits[[n.sets]] <- list("found"=start, "it.scores"=c())
			tmp <- start
			start <- c()
		}
    
	  # Print the initialized set.
		if(echo){
			print("Initialized")
			print(tmp)
		}
		
		diag <- FALSE
		
		# Run search procedure with specified values until we reach a fixed point or
		#  approximate fixed point (oscillating).
		valid.termination <- FALSE
		reinitialize <- FALSE
		while(!valid.termination & difftime(Sys.time(), starttime, units = "hours") < max.time){
		  # If we've already called run_DCM on this initial set and got an invalid
		  #  termination, we need to reinitialize to avoid getting the same result.
		  if(reinitialize){
		    inits[[n.sets]] <- init_DCM(M1, M2, k=est.size, del=startdels)
		    tmp <- inits[[n.sets]]$found
		    if(echo){
		      print("Reinitialized...")
		      print(tmp)
		    }
		  }
		  # Run the search procedure.
		  DCM <- run_DCM(M1, M2, seed=tmp, echo=echo, max.iter=max.iter, 
		                 max.size=max.size, validation=validation, alpha=alpha)
		  startdels <- c(startdels, DCM$startdels)
		  res <- DCM$found
		  k <- length(res)
      # Decide whether the search procedure terminated in a valid way.
		  valid.termination <- DCM$reason.for.terminating %in% 
		    c("Reached fixed point", "Oscillating between similar sets")
		  # If we've made it this far through the loop and need to repeat the loop,
		  #  we have to reinitialize to make sure we don't get the same result.
		  reinitialize <- TRUE
		}
		
		# Create a validation document if desired.
		if(validation & k > 0){
		  # Make folder for this search.
		  validation.subdir <- file.path(validation.dir, paste0("DCM_Validation_Search", n.sets))
		  dir.create(validation.subdir, showWarnings=FALSE)
		  create_validation_doc(DCM, inits[[n.sets]], preproc.steps, data.flags, 
		                        validation.subdir)
		} else if(validation){
		  data.flags <- append(data.flags, "No set was found by iterative testing procedure.")
		}
		
		# Ignore groups that are too small to be significant
		if(k > 10){
			# Check for diag blocks
			if(res[1] == 0){
				diag <- TRUE
				cut <- res[2]
				res <- res[-c(1,2)]
			}
			if(diag){
  			res1 <- res[1:cut]
  			res2 <- res[-(1:cut)]
  			if(echo){
  				print(sprintf("Found an off-diagonal block group of size %i", k))
  			}
			} else {
				# Announce result
				if(echo){
					print(sprintf("Found a group of size %i", k))
					print(sprintf("cor 1 = %f", DCM$mc1))
					print(sprintf("cor 2 = %f", DCM$mc2))
				}
			}
			
			# Residualize remaining data for continued search
			M1[res,] <- resid_DCM(M1[res,], QN=QN)
			M2[res,] <- resid_DCM(M2[res,], QN=QN)
			
			# Save results
			if(diag){
				DC_sets[[n.sets]] <- list("Block Diag", idcs[res1], idcs[res2])
				meanCor1[[n.sets]] <- 0
				meanCor2[[n.sets]] <- 0
			}else{
				DC_sets[[n.sets]] <- idcs[res]
				meanCor1[[n.sets]] <- DCM$mc1
				meanCor2[[n.sets]] <- DCM$mc2
			}
			iterations[[n.sets]] <- DCM$its
			time[[n.sets]] <- DCM$time
			if(validation){
			  it.sets[[n.sets]] <- DCM$it.sets
			  it.pvals[[n.sets]] <- DCM$it.pvals
			  it.test.stats[[n.sets]] <- DCM$it.test.stats
			  it.std.errs[[n.sets]] <- DCM$it.std.errs
			  it.diff.cor[[n.sets]] <- DCM$it.diff.cor
			}

			# Count number of sets found
			n.sets <- n.sets + 1
		} else {
			# If group is degenerate, do not try it again
			startdels <- c(startdels, tmp)
		}# if(k > 10)
		
		# Check total runtime
		tottime <- difftime(Sys.time(), starttime, units = "hours")
	} # while(length(dels) < p-est.size & n.sets < max.groups)
	
	# Print reason for algorithm halting
	if(echo){
		if(tottime >= max.time){
			print("Timed out.")
		}else if(length(startdels) > p-est.size){
			print("Exhausted all searches.")
		}else if(n.sets > max.groups){
			print(sprintf("%i groups found.", max.groups))
		}
	}

	# Return results and properties each result
	if(validation){
	  return(list(DC_sets = DC_sets, 
	              iterations = iterations, 
	              time = time, 
	              meanCor1 = meanCor1, 
	              meanCor2 = meanCor2,
	              it.pvals = it.pvals,
	              it.test.stats = it.test.stats,
	              it.std.errs = it.std.errs,
	              it.sets = it.sets,
	              inits = inits,
	              indices = idcs))
	} else {
	  return(list(DC_sets = DC_sets, 
	              iterations = iterations, 
	              time = time, 
	              meanCor1 = meanCor1, 
	              meanCor2 = meanCor2, 
	              indices = idcs))
	}
}
