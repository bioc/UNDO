\name{dimension_reduction}
\alias{dimension_reduction}
\title{
Dimension reduction function
}
\description{
When the number of input samples is larger than 2, this function is called to reduce the dimension to 2 by using PCA.
}
\usage{
dimension_reduction(X)
}
\arguments{
  \item{X}{
 gene expression data matrix
}
}
\value{

\item{X}{}
\item{dimenMatrix}{the dimension reduction matrix used to recover the mixing matrix for all the samples}
}
\author{
Niya Wang (wangny@vt.edu)
}

\examples{
X <- matrix(runif(5000),1000,5)
dimenResult <- dimension_reduction(X)
}

