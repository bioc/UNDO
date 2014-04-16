marker_gene_selection <-
function(X,lowper,highper,epsilon1,epsilon2) {
   
  data_size <- dim(X)[1]
  
  ## If X contains more than 2 samples,  the dimension reduction funcion is called.
  if(dim(X)[2]>2) {
    dimenreduc <- dimension_reduction(X)
    X <- unlist(dimenreduc[[1]])
    dimenMatrix <- unlist(dimenreduc[[2]])
    
  } else dimenMatrix <- NULL
  
  
  sigNorm <- matrix(0,1,data_size)
  for (i in 1:data_size) {
    sigNorm[i] <- norm(matrix(X[i,]),"F") 
  }
  
  sigOrder <- sort(sigNorm,index.return=TRUE)
  Xfilter <- X[sigOrder$ix[(lowper*data_size+1):(data_size-highper*data_size)],]
  Xfnorm <- sigNorm[sigOrder$ix[(lowper*data_size+1):(data_size-highper*data_size)]]
  
  k <- Xfilter[,2]/Xfilter[,1]
  eps1 <- epsilon1
  eps2 <- epsilon2
  
  MG1 <- which((k>=max(k)-abs(eps1*max(k)))&k<=max(k))
  MG2 <- which((k<=min(k)+abs(eps2*min(k)))&k>=min(k))

  
  ## handle the situation that one or two samples are pure
  if (max(k)==Inf) {
    a1 <- c(1,0)
  } else {
    
    if (length(MG1)==1) {
      a1 <- Xfilter[MG1,]/Xfnorm[MG1]
    } else {
      a1 <- colMeans(Xfilter[MG1,]/Xfnorm[MG1])
    }
    
  }
  
  if (min(k)==0){
    a2 <- c(0,1)
  } else {
    if (length(MG2)==1) {
      a2 <- Xfilter[MG2,]/Xfnorm[MG2]
    } else {
      a2 <- colMeans(Xfilter[MG2,]/Xfnorm[MG2])
    }
  }
  
  
  markergene <- list(a1,a2,MG1,MG2,dimenMatrix)
  markergene
    
}
