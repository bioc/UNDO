\name{calc_E1}
\alias{calc_E1}
\title{
    function calculating the E1 measurement
}
\description{
   A function used to calculate the E1 measurement when the real mixing matrix is provided

}
\usage{
calc_E1(A, Aest)
}
\arguments{
  \item{A}{real mixing matrix
}
  \item{Aest}{
estimated mixing matrix
}
}

\value{
E1 measurement (numeric)
}

\author{
Niya Wang <wangny@vt.edu>
}

\examples{
A <- matrix(runif(4),2,2)
Aest <- matrix(runif(4),2,2)
E1 <- calc_E1(A,Aest)  # to calculate the similarity of two randowm 2*2 matrix

}

