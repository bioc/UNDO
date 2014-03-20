gene_expression_input <-
function(X) {
  
  
  ## make sure the input is not empty
  if (is.null(X)) 
    stop(" Missing the mixture dataset")
  
  if (dim(X)[2]==1) 
    stop('ERROR: Please input more than one samples.')
  
  if (dim(X)[2]>2){
    X <- X[,1:2]
    warning('Only the first two samples will be deconvoluted.')
  }
    
  
  ## make sure the input is nonnegative.
  if (min(X)<0) {
    stop('ERROR: The value of gene expression should be positive.')
  }
  
  
  ## delete missing value
  if (nrow(X)!=nrow(na.omit(X))) {
    X <- na.omit(X)
    warning('There is missing value in the expression matrix.')
    warning('NA in X is removed.')
  }
  
  ## The correlation of two samples is 1.
  if ((ncol(X)==2)&&(corr(X)==1)) {
    stop('Correlation between the two samples is 1.')
  }

  return(X)
}
