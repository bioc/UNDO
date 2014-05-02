mixing_matrix_computation <-
function(X,a1,a2,dimenMatrix){

  
  Aest <- cbind(a1,a2)
  
  
  if(!is.null(dimenMatrix)){
    Atrans <- ginv(t(dimenMatrix)%*%dimenMatrix)%*%t(dimenMatrix)%*%Aest
    scale <- ginv(Atrans)%*%matrix(1,dim(X)[2],1)
    scale <- as.vector(scale)
    Aest <- Atrans%*%diag(scale)
  } else {
    
    scale <- ginv(Aest)%*%matrix(1,2,1)
    scale <- as.vector(scale)
    Aest <- Aest%*%diag(scale)
  }
  
  Aest <- Aest/matrix(rep(rowSums(Aest),2),nrow(Aest),2)

  Sest <- matrix(0,ncol=dim(Aest)[2],nrow=dim(X)[1])
  for (i in 1:nrow(X)){
    Sest[i,] <- coef(nnls(Aest,X[i,]))
  }

   
  ## create a new folder for results
  folder <- substr(date(),1,10)
  folder <- paste("Result",folder)
  folder <- paste(folder,floor(runif(1,1,100000)))
  dir.create(folder)
  
  write.table(Aest,file=paste(folder,'/Aest.txt',sep=''),col.names=FALSE,row.names=FALSE)
  write.table(Sest,file=paste(folder,'/Sest.txt',sep=''),col.names=FALSE,row.names=FALSE)

  result <- list(Aest,Sest)
  return(result)
  
  
  
}
