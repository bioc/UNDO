### R code from vignette source 'DeconRNASeq.Rnw'

###################################################
### code chunk number 1: set_width
###################################################
options( width = 60 )


###################################################
### code chunk number 2: packageLoad
###################################################
library(UNDO)
#load tumor stroma mixing tissue samples
data(NumericalMixMCF7HS27)
X <- NumericalMixMCF7HS27
#load mixing matrix for comparison
data(NumericalMixingMatrix)
A <- NumericalMixingMatrix


###################################################
### code chunk number 3: mixing matrix
###################################################
A


###################################################
### code chunk number 4: Deconvolution
###################################################
two_source_deconv(X,lowper=0.4,highper=0.1,epsilon1=0.01,epsilon2=0.01,
                              A,return=0)


###################################################
### code chunk number 5: Scatterplot
###################################################			
# compute the estimated pure source expressions
result <- two_source_deconv(X,lowper=0.4,highper=0.1,epsilon1=0.01,epsilon2=0.01,
                              A,return=1)
Sest <- result[[3]]
#load pure tumor stroma expressions 
data(PureMCF7HS27)
S <- exprs(PureMCF7HS27)
#draw the scatter plots between pure and estimated expressions of MCF7 and HS27
plot(S[,1],Sest[,1],main="MCF7" ,xlab="Estimated expression", ylab="Measured expression", pch=1, col="turquoise", cex=0.5)


###################################################
### code chunk number 5: Scatterplot2
###################################################		
plot(S[,2],Sest[,2],main="HS27" ,xlab="Estimated expression", ylab="Measured expression", pch=1, col="turquoise", cex=0.5)


											