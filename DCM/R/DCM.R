DCM <- function(M1, 
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
                validation = FALSE,
                validation.dir = ''){
	
	################### Prepare Data #############################################
	
	prep = prepData_DCM(M1, M2, QN, resid.full, echo, strict, validation)
	M1 = prep$M1
	M2 = prep$M2
	idcs = prep$idcs
	if(validation){
	  preproc.steps <- prep$preproc.steps
	  data.flags <- prep$data.flags
	}
	
	p = nrow(M1)
	
	################### Run Algorithm ############################################

	# Prepare to store results
	DC_sets = list()
	iterations = list()
	time = list()
	meanCor1 = list()
	meanCor2 = list()
	it.p.vals = list()
	it.test.stats = list()
	it.sets = list()
	inits = list()

	# Initialize parameters.
	starttime = Sys.time()
	tottime = 0
	n.sets = 1
	startdels = c() # Vars to delete.
	
	# Keep searching until all rows or exhausted, or max groups or time is reached
	while(length(startdels) < p-est.size & n.sets <= max.groups & tottime < max.time){
	  # Find initial starting set, ignoring already tried starting points.
		if(length(start) < 5){
		  if(validation & length(start) > 0){
		    data.flags <- append(data.flags, "Desired starting set too small (< 5) so it is being ignored.")
		  }
			inits[[n.sets]] <- init_DCM(M1, M2, k=est.size, del=startdels)
			tmp = inits[[n.sets]]$found
		}else{
			tmp = start
			start = c()
		}
	
		if(echo){
			print("Initialized")
			print(tmp)
		}
		
		diag = FALSE
		
		# Run search procedure with specified values
		DCM = run_DCM(M1, M2, seed = tmp, echo = echo, max.iter = max.iter)
		startdels = c(startdels, DCM$startdels)
		res = DCM$found
		k = length(res)
		
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
			M1[res,] = resid_DCM(M1[res,], QN = QN)
			M2[res,] = resid_DCM(M2[res,], QN = QN)
			
			# Save results
			if(diag){
				DC_sets[[n.sets]] = list("Block Diag", idcs[res1], idcs[res2])
				meanCor1[[n.sets]] = 0
				meanCor2[[n.sets]] = 0
			}else{
				DC_sets[[n.sets]] = idcs[res]
				meanCor1[[n.sets]] = DCM$mc1
				meanCor2[[n.sets]] = DCM$mc2
			}
			iterations[[n.sets]] = DCM$its
			time[[n.sets]] = DCM$time

			# Count number of sets found
			n.sets = n.sets + 1
		}else if(!is.null(res)){
			# If group is degenerate, do not try it again
			startdels = c(startdels, tmp)
		}# if(k > 10)
		
		# Check total runtime
		tottime = difftime(Sys.time(), starttime, units = "hours")
		
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
	
	# Create a validation document if desired.
	if(validation){
	  create_validation_doc(p.vals, 
	                        test.stats, 
	                        preproc.steps, 
	                        data.flags, 
	                        validation.dir)
	}
	
	# Return results and properties each result
	return(list(DC_sets = DC_sets, 
	            iterations = iterations, 
	            time = time, 
	            meanCor1 = meanCor1, 
	            meanCor2 = meanCor2, 
	            indices = idcs))
}
