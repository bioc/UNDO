dimension_reduction <-
function(X) {
    ##eigenvalue decomposition is used to reduce dimension
    
    X <- t(X)
    Xmean <- rowMeans(X)
    Xmeanrm <- X-Xmean  
    r <- eigen(Xmeanrm%*%t(Xmeanrm))
    dimenMatrix <- t(r$vectors[,1:2])
    Xt <- dimenMatrix%*%X
    X <- t(Xt)
    
    dimenreduc <- list(X,dimenMatrix)
    dimenreduc
  }
