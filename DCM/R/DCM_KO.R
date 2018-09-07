# TITLE: DCM_KO.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 9/6/18
# DATE MODIFIED: 9/6/18

DCM_KO <- function(mat.1,
                   mat.2,
                   alpha = 0.05,
                   init.size = round(nrow(mat.1)/10),
                   max.iter = 50,
                   max.groups = 1,
                   max.time = 1,
                   QN = FALSE,
                   debug = TRUE){
  
  # DCM_KO runs the DCM algorithm on the input data matrices, mat.1 and mat.2. 
  #
  # Args:
  #   -mat.1: Data matrix for group 1. Rows should correspond to variables and
  #     columns should correspond to samples.
  #   -mat.2: Data matrix for group 2. Orientation same as that of mat.1.
  #   -alpha: Significance level to use for iterative testing.
  #   -init.size: The initial size of the DC set to be used during
  #     initialization.
  #   -max.iter: Maximum number of iterations to perform during iterative
  #     testing.
  #   -max.groups: Maximum number of DC sets to search for.
  #   -max.time: Maximum amount of time in hours for the algorithm to run.
  #   -QN: When switched on, will quantile normalize the data.
  #   -debug: When switched on, will output iteration-level data.
  #
  # Returns:
  #   List containing the DC sets found as well as relevant data from the run.
  
  # Setup debugging directory.
  if(debug){
    debug.dir <- file.path(getwd(), "Debug_Output")
    dir.create(debug.dir, showWarnings=FALSE)
    sink(file=file.path(debug.dir, "Console_Output.txt"))
  }
  
  # Prepare data.
  prepared.data <- prepData_DCM(mat.1, 
                                mat.2, 
                                QN=FALSE, 
                                resid.full=FALSE,
                                echo=TRUE,
                                strict="low")
  mat.1     <- prepared.data$mat.1
  mat.2     <- prepared.data$mat.2
  good.idcs <- prepared.data$good.idcs

  p <- nrow(mat.1)
  
  # Prepare variables for storing data.
  DC_sets    <- list()
  iterations <- list()
  time       <- list()
  meanCor1   <- list()
  meanCor2   <- list()
  if(debug){
    init <- list()
  }
  
  # Initialize variables.
  start.time <- Sys.time()
  tot.time   <- 0
  search.num <- 1
  startdels  <- c()
  converged  <- TRUE
  
  # Search loop: repeatedly run DCM, residualize and run again as long as max
  #  number of groups and time has not been reached.
  while(search.num <= max.groups & tot.time < max.time & converged){
    # Find initial set.
    initialization <- init_DCM(mat.1, mat.2, k=init.size)
    if(debug){
      save(initialization, file=filePath(debug.dir, "Initialization.RData"))
    }
    init.set <- initialization$found
    
    # Run iterative search procedure.
    DCM <- run_DCM_KO(mat.1, mat.2, seed=init.set, max.iter=max.iter)
    if(debug){
      save(DCM, file=filePath(debug.dir, paste0("DCM_", search.num, ".RData")))
    }
    DC_set <- DCM$found
    
    # Check whether algorithm converged.
    converged <- (DCM$its < max.iter) & (length(DC_set) > 1)
    
    # If converged, residualize data.
    if(converged){
      mat.1[DC_set,] <- resid_DCM(mat.1[DC_set,], QN=QN)
      mat.2[DC_set,] <- resid_DCM(mat.2[DC_set,], QN=QN)
    }
    
    # Store iteration data.
    DC_sets[[search.num]]    <- good.idcs[DC_set]
    meanCor1[[search.num]]   <- DCM$mc1
    meanCor2[[search.num]]   <- DCM$mc2
    iterations[[search.num]] <- DCM$its
    time[[search.num]]       <- DCM$time
    if(debug){
      init[[search.num]] <- initialization
    }
    
    # Increment search number.
    search.num <- search.num + 1
    
    # Check total runtime.
    tot.time <- difftime(Sys.time(), start.time, units="hours")
    
    # Print reason for halting.
    if(tot.time >= max.time){
      print("Timed out.")
    } else if(search.num > max.groups){
      print(sprintf("%i groups found.", max.groups))
    }
  }
  
  # Remove sink.
  if(debug){sink()}  
  
  # Return results and relevant data.
  if(debug){
    return(list(DC_sets        = DC_sets, 
                iterations     = iterations, 
                time           = time, 
                meanCor1       = meanCor1, 
                meanCor2       = meanCor2, 
                indices        = good.idcs,
                initialization = init))
  } else {
    return(list(DC_sets        = DC_sets, 
                iterations     = iterations, 
                time           = time, 
                meanCor1       = meanCor1, 
                meanCor2       = meanCor2, 
                indices        = good.idcs))
  }  
}