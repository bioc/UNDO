\name{mixing_matrix_computation}
\alias{mixing_matrix_computation}
\title{
Calculate and scale the mixing matrix
}
\description{
Calculate the mixing matrix based on the output from marker_gene_selection(), and scale the mixing matrix to make the sum of proportions from tumor and stroma equal to 1.The pure expression levels of tumor and stroma are also computed. 
}
\usage{
mixing_matrix_computation(X, a1, a2, dimenMatrix)
}
\arguments{
  \item{X}{
Gene expression data matrix
}
  \item{a1}{
The slope of marker genes in source 1
}
  \item{a2}{
The slope of marker genes in source 2
}
  \item{dimenMatrix}{
The dimention reduction matrix used to recover mixing matrix for all the samples
}
}
\value{
\item{Aest}{estimated mixing matrix}
\item{Sest}{estimated pure gene expression of two sources}
}
\author{
Niya Wang (wangny@vt.edu)
}
\examples{
a1<- matrix(runif(2),2,1)
a2<- matrix(runif(2),2,1)
X <- 1000*matrix(runif(20000),10000,2)
dimenMatrix <- NULL
Deconv <- mixing_matrix_computation(X, a1, a2, dimenMatrix)
}

