\name{UNDO-package}
\alias{UNDO-package}
\alias{UNDO}
\docType{package}
\title{
Implementation of UNDO (unsupervised deconvolution of tumor-stromal mixed expressions)
}
\description{
This package contains main function "two_source_deconv" to implement the deconvolution of mixed tumor-stromal expressions in a completely unsupervised way. The prior knowledge of mixing matrix or pure expression is not needed. The package detects marker genes and calculate the mixing matrix and pure expressions automatically.
}
\details{
\tabular{ll}{
Package: \tab UNDO\cr
Type: \tab Package\cr
Version: \tab 1.5.2\cr
Date: \tab 2014-04-30\cr
License: \tab GPL version 2 or later\cr
}
two_source_deconv(ExpressionData,lowper=0.4,highper=0.1,epsilon1=0.01,epsilon2=0.01,
                              A=NULL,return=0)
}
\author{
Niya Wang <wangny@vt.edu>
}
\keyword{ package }
\keyword{methods}
\examples{
data(NumericalMixMCF7HS27)
X <- NumericalMixMCF7HS27
deconvResult <- two_source_deconv(X, lowper = 0.4, highper = 0.1, epsilon1 = 0.1, epsilon2 = 0.1, A = NULL, return = 0)
}
