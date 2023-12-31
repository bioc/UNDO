two_source_deconv <-
function(ExpressionData,lowper=0.4,highper=0.1,epsilon1=0.01,epsilon2=0.01,
                              A=NULL,S1=NULL,S2=NULL,return=0){
  
  ## check the input
  if (lowper<0||lowper>1||highper<0||highper>1){
    stop("ERROR: the gene removing percentage should be between 0 and 1.")
  }
  
  if (lowper+highper>1) stop('ERROR: the sum of lowper and highper should be smaller than 1.')
  
  if (epsilon1<0||epsilon2<0) stop('ERROR: the eps1/eps2 shoule be positive.')
  
  if (return!=0&&return!=1) stop('ERROR: return shoule be 0 or 1.')
  
  
  if (class(ExpressionData)=="matrix"){
    X <- ExpressionData
  } else if (class(ExpressionData)=="ExpressionSet"){
    X <- exprs(ExpressionData)
  } else {
    stop("ERROR: Please input expression data matrix or ExpressionSet data.")
  }
  
  
  ## select marker genes
  X <- unlist(gene_expression_input(X))
  markergene <- marker_gene_selection(X,lowper,highper,epsilon1,epsilon2)
  
  a1 <- unlist(markergene[[1]])
  a2 <- unlist(markergene[[2]])
  dimenMatrix <- unlist(markergene[[5]])
  
  deconvresult <- mixing_matrix_computation(X,a1,a2,dimenMatrix)
  Aest <- unlist(deconvresult[[1]])
  Sest <- unlist(deconvresult[[2]])
  
  
  ## check whether the input A is valid
  if(!is.null(A)) {
    if(dim(A)[1]!=dim(Aest)[1]||dim(A)[2]!=dim(Aest)[2]){
      warning('The size of input mixing matrix is different from estimated mixing matrix.')
      E1 <- NULL      
    }
    
    if(dim(A)[1]==dim(Aest)[1]&&dim(A)[2]==dim(Aest)[2]) {
      E1<- min(calc_E1(A,Aest),calc_E1(A,Aest[c(2,1),]))
    } else {
      E1<-NULL
    }
  } else E1 <- NULL
  
  
  ## calculate the correlation between estimated and measured expression
  if(!is.null(S1)){
    S_corr1 <- max(corr(cbind(S1,Sest[,1])),corr(cbind(S1,Sest[,2])))
  } else {
    S_corr1 <-NULL
  }
  
  if(!is.null(S2)){
    S_corr2 <- max(corr(cbind(S2,Sest[,1])),corr(cbind(S2,Sest[,2])))
  } else {
    S_corr2 <- NULL
  }
  
  ## output
  if (return==1) list(Estimated_Mixing_Matrix=Aest,E1=E1,S1_correlation=S_corr1,
                      S2_correlation=S_corr2,Sest)
  else list(Estimated_Mixing_Matrix=Aest,E1=E1,S1_correlation=S_corr1,S2_correlation=S_corr2)
  
}
