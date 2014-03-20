calc_E1 <-
function(A,Aest){
  
  P <- A%*%ginv(Aest)
  
  m <- dim(P)[1]
  n <- dim(P)[2]
  
  if (m!=n) {
    stop('ERROR: column number should be equal to row number')
    return
  }
  
  E11 <- 0
  
  for (i in 1:m){
    max_pi <- -1000
    for (k in 1:m){
      if (abs(P[i,k])>max_pi){
        max_pi <- abs(P[i,k])
      }
    }
    sum1 <- 0
    for (j in 1:m){
      sum1 <- sum1+abs(P[i,j])/max_pi
    }
    E11 <- E11+(sum1-1)
  }
  
  E12 <- 0
  
  for (j in 1:m){
    max_pj <- -1000
    for (k in 1:m){
      if (abs(P[k,j])>max_pj){
        max_pj <- abs(P[k,j])
      }
    }
    sum2 <- 0
    for (i in 1:m){
      sum2 <- sum2+abs(P[i,j])/max_pj
    }
    E12 <- E12+(sum2-1)
  }
  E1 <- E11+E12 
  E1
  
}
