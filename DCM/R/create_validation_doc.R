#' create_validation_doc
#'
#' Creates pdf and txt documents containing validation diagnostics for the DCM 
#' algorithm.
#'
#' @param validation.data Output of \code{run_DCM} with \code{validation==TRUE}.
#' A list containing lists, \code{it.pvals}, \code{it.test.stats}, 
#' \code{it.std.errs}, \code{it.sets}, and \code{it.diff.cor}.
#' @param init.data Output of \code{init_DCM}. A list containing a list, 
#' \code{it.scores}.
#' @param preproc.steps Vector containing strings with all preprocessing steps 
#' performed.
#' @param data.flags Vector containing strings with all data flags thrown by 
#' the algorithm.
#' @param validation.dir String giving the directory in which the pdf and txt
#' docs should be saved.
#'
#' @return None
#' 
#' @export
create_validation_doc <- function(validation.data, init.data,
                                  preproc.steps = list(), data.flags = list(),
                                  validation.dir = ''){
  sink(file.path(validation.dir, "DCM_Diagnostics_Text.txt"))
  ##############################################################################
  # List preprocessing steps.
  print("The following steps were performed in preprocessing:")
  for(step in preproc.steps){
    print(step)
  }
  cat("\n\n\n")
  
  ##############################################################################
  # List data flags.
  print("The following data flags were thrown during preprocessing:")
  for(flag in data.flags){
    print(flag)
  }
  cat("\n\n\n")
  sink()
  ##############################################################################
  # Extract data from validation.data.
  pvals <- validation.data$it.pvals
  test.stats <- validation.data$it.test.stats
  std.errs <- validation.data$it.std.errs
  sets <- validation.data$it.sets
  diff.cor <- validation.data$it.diff.cor
  init.scores <- init.data$it.scores
  
  # Plot diagnostics.
  n.iter <- length(pvals)
  pdf(file.path(validation.dir, "DCM_Diagnostic_Plots.pdf"))
  
  # Plot scores at each iteration of initialization.
  plot(init.scores, main="Score across iterations during initialization", 
       xlab="Iteration", ylab="Score")
  
  # Plot p-values at each iteration.
  for(i in 1:n.iter){
    hist(pvals[[i]], main=paste0("P-values from iteration ", i), 
         xlab="P-value")
  }
  
  # Plot test statistics at each iteration.
  for(i in 1:n.iter){
    hist(test.stats[[i]], main=paste0("Test statistics from iteration ", i), 
         xlab="Test statistic")
  }
  
  # Plot standard errors at each iteration.
  for(i in 1:n.iter){
    hist(std.errs[[i]], main=paste0("Standard errors from iteration ", i),
         xlab="Standard error")
  }
  
  # Plot set sizes at each iteration.
  plot(sapply(sets, function(x){length(x)}), main="Set size across iterations", 
       xlab="Iteration", ylab="Set size")
  
  # Plot Jaccard index at each iteration.
  plot(sapply(2:length(sets), function(ind){
    return(length(intersect(sets[[ind-1]], sets[[ind]]))/length(union(sets[[ind-1]], sets[[ind]])))
  }), main="Jaccard index (#(A intersect B)/#(A union B)) across iterations", 
  xlab="Iteration", ylab="Jaccard index")
  
  # Plot differential correlation at each iteration.
  plot(diff.cor, main="Mean differential correlation of active set across iterations", xlab="Iteration", ylab="Mean differential correlation of active set")
  
  dev.off()
}