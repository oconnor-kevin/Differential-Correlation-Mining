# TITLE: makeVars_boot.R
# AUTHOR: Kevin O'Connor
# DATE CREATED: 10/18/18
# DATE MODIFIED: 10/18/18

makeVars_boot <- function(X1, X2, A, boot.iter=1000){
  # Estimate approximate standard error of the test statistic via bootstrap.
  
  n1 <- ncol(X1)
  n2 <- ncol(X2)
  p <- nrow(X1)
  dat <- cbind(X1, X2)
  
  test.stats <- c()
  rand.rows <- sample(1:p, boot.iter)
  for(rand.r in rand.rows){
    perm <- sample.int(n1+n2)
    X1.perm <- dat[c(A,rand.r),perm[1:n1]]
    X2.perm <- dat[c(A,rand.r),perm[(n1+1):n2]]
    
    #test.stats <- c(test.stats, mean(cor(t(X1.perm)) - cor(t(X2.perm))[A, rand.r]))
    test.stats <- c(test.stats, mean(cor(t(X1.perm))[length(A)+1, 1:length(A)] - cor(t(X2.perm))[length(A)+1, 1:length(A)]))
  }
  
  return(var(test.stats))
}