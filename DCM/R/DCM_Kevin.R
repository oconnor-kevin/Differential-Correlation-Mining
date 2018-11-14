# TITLE: DCM_Kevin.R
# AUTHOR: Kevin O'Connor(, Kelly Bodwin)
# DATE MODIFIED: 8/22/18

# DCM_Kevin contains added debugging functionality including storing
#  data from each iteration.

DCM_Kevin <- function(M1, 
                      M2, 
                      max.groups = 5, 
                      max.iter = 50, 
                      max.time = 1, 
                      est.size = nrow(M1)/10, 
                      alpha = 0.05, 
                      start = c(), 
                      QN = FALSE, 
                      resid.full = FALSE, 
                      echo = FALSE, 
                      strict = "low",
                      initialize = TRUE,
                      debug = TRUE,
                      debug.dir = "",
                      new.p.val.method = FALSE){
	
  ################### Setup Debugging ##################
  
  if(debug){
    # Create debug output directory.
    if(debug.dir == ""){
      debug.dir <- file.path(getwd(), "Debug_Output")
    }
    dir.create(debug.dir, showWarnings=FALSE)
    
    sink(file=file.path(debug.dir, "Console_Output.txt"))
  }
  
  
	################### Prepare Data #####################
	
	prep = prepData_DCM(M1, M2, QN, resid.full, echo, strict)
	M1 = prep$M1
	M2 = prep$M2
	idcs = prep$idcs
	
	p = nrow(M1)
	
	################### Run Algorithm ####################

	# Prepare to store results
	DC_sets = list()
	iterations = list()
	time = list()
	meanCor1 = list()
	meanCor2 = list()
	if(debug){
  	init = list()
	}

	# Initialize
	starttime = Sys.time()
	tottime = 0
	ans = 1
	kval = est.size
	startdels = c()
	converged = TRUE
	
	# Keep searching until all rows or exhausted, or max groups or time is reached
	while(length(startdels) < p-kval & ans <= max.groups & tottime < max.time & converged){
		
		# Find initial starting set, ignoring already tried starting points
		if(length(start) < 5){
			if(initialize){
  		  initialization <- init_DCM(M1, M2, k=kval, del=startdels)
	  	  if(debug){
		      save(initialization, file=filePath(debug.dir, "Initialization.RData"))
		    }
		    tmp <- initialization$found
			} else {
			  tmp <- sample(1:nrow(M1), kval)
			  initialization <- c()
			}
		}else{
			tmp = start
			start = c()
		}
	
		if(echo){
			print("Initialized")
			print(tmp)
		  #print(initialization)
		}
	  
		diag = FALSE
		
		# Run search procedure with specified values
		DCM = run_DCM(M1, M2, seed = tmp, echo = echo, max.iter = max.iter)
		if(debug){
		  save(DCM, file=filePath(debug.dir, paste0("DCM_", ans, ".RData")))
		}
		print(DCM)
		
		startdels <- c(startdels, DCM$startdels)
		res       <- DCM$found
		k         <- length(res)

		# Check whether algorithm converged.
		converged <- (DCM$its < max.iter)
				
		# Ignore groups that are too small to be significant
		if(k > 10){
			# Check for diag blocks
			if(res[1] == 0){
				diag = TRUE
				cut = res[2]
				res = res[-c(1,2)]
			}
			
			if(diag){
				res1 = res[1:cut]
				res2 = res[-(1:cut)]

				if(echo){
					print(sprintf("Found an off-diagonal block group of size %i", k))
				}
			}else{
				# Announce result
				if(echo){
					print(sprintf("Found a group of size %i", k))
					print(sprintf("cor 1 = %f", DCM$mc1))
					print(sprintf("cor 2 = %f", DCM$mc2))
				}
			}
			
			# Residualize remaining data for continued search
		  if(converged){
		    M1[res,] <- resid_DCM(M1[res,], QN=QN)
		    M2[res,] <- resid_DCM(M2[res,], QN=QN)
		  }
			
			# Save results
			if(diag){
				DC_sets[[ans]]  <- list("Block Diag", idcs[res1], idcs[res2])
				meanCor1[[ans]] <- 0
				meanCor2[[ans]] <- 0
			}else{
				DC_sets[[ans]]  <- idcs[res]
				meanCor1[[ans]] <- DCM$mc1
				meanCor2[[ans]] <- DCM$mc2
			}
			iterations[[ans]] <- DCM$its
			time[[ans]]       <- DCM$time
			if(debug){
  			init[[ans]]     <- initialization
			}
			
			# Count number of sets found
			ans <- ans+1
		}else{ # Else if found set smaller than 10 variables
			# If group is degenerate, do not try it again
			startdels <- c(startdels, tmp)
		}
		
		# Check total runtime
		tottime <- difftime(Sys.time(), starttime, units = "hours")
	} # while(length(dels) < p-kval & ans < max.groups)
	
	# Print reason for algorithm halting
	if(echo){
		if(tottime >= max.time){
			print("Timed out.")
		}else if(length(startdels) > p-kval){
			print("Exhausted all searches.")
		}else if(ans > max.groups){
			print(sprintf("%i groups found.", max.groups))
		}
	}
	
	# Remove sink.
	if(debug){sink()}
	
	# Return results and properties each result
	if(debug){
  	return(list(DC_sets = DC_sets, 
	              iterations = iterations, 
	              time = time, 
	              meanCor1 = meanCor1, 
	              meanCor2 = meanCor2, 
	              indices = idcs,
	              initialization = init))
	} else {
	  return(list(DC_sets = DC_sets, 
	              iterations = iterations, 
	              time = time, 
	              meanCor1 = meanCor1, 
	              meanCor2 = meanCor2, 
	              indices = idcs))
	}
}