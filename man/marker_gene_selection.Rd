\name{marker_gene_selection}
\alias{marker_gene_selection}
\title{
Select marker genes in two sources
}
\description{
Select the marker genes in tumor and stroma in an unsupervised way
}
\usage{
marker_gene_selection(X, lowper, highper, epsilon1, epsilon2)
}
\arguments{
  \item{X}{
gene expression data
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
}

\value{

 \item{a1}{The slope of marker genes in source 1}
 \item{a2}{The slope of marker genes in source 2}
 \item{MG1}{The gene list of marker genes in source 1}
 \item{MG2}{The gene list of marker genes in source 2}

}
\author{
Niya Wang (wangny@vt.edu)
}
\examples{
X <- matrix(runif(20000),10000,2)
MG_set <- marker_gene_selection(X, 0.4, 0.1, 1, 1)
}

