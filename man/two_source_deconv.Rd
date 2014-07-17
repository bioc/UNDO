\name{two_source_deconv}
\alias{two_source_deconv}
\title{
Main function to call other subfunction to deconvolute the mixed expression data.
}
\description{
This is the main function that is to call all the other subfunctions and realize the deconvolution of mixed expression data. When the real mixing matrix exist, it will also compare the estimated mixing matrix and real mixing matrix and give the E1 measurement.
}
\usage{
two_source_deconv(ExpressionData, lowper = 0.4, highper = 0.1, epsilon1 = 0.01, epsilon2 = 0.01, A = NULL, S1=NULL, S2=NULL, return = 0)
}
\arguments{
  \item{ExpressionData}{
gene expression data matrix/ExpressionSet object
}
  \item{lowper}{
The percentage of genes the user wants to remove with lowest norm. The range should be between 0 and 1.
}
  \item{highper}{
The percentage of genes the user wants to remove with highest norm.The range should be between 0 and 1.
}
  \item{epsilon1}{
Influence the number of marker genes. With increasing of epsilon1, the number marker genes in source 1 will increase. The value should be positive.
}
  \item{epsilon2}{
Influence the number of marker genes. With increasing of epsilon1, the number marker genes in source 2 will increase. The value should be positive.
}
  \item{A}{
real mixing matrix if existing
}
  \item{S1}{
Pure expression profile of first source if existing
}
  \item{S2}{
Pure expression profile of second source if existing
}
  \item{return}{
if it is equal to 0, do not return estimated S; otherwise, return the estimated S.
}
}

\value{
\item{Aest}{estimated mixing matrix}
\item{E1}{E1 measurement between real and estimated mixing matrix}
}
\author{
Niya Wang (wangny@vt.edu)
}
\examples{

data(NumericalMixMCF7HS27)
X <- NumericalMixMCF7HS27
deconvResult <- two_source_deconv(X, lowper = 0.4, highper = 0.1, epsilon1 = 0.1, epsilon2 = 0.1, A = NULL, S1=NULL,S2=NULL, return = 0)

}

