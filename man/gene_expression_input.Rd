\name{gene_expression_input}
\alias{gene_expression_input}
\title{
 Detect whether the input gene expression data are valid
}
\description{
Check the input gene expression data to see whether they are nonempty, nonnegative, etc.
}
\usage{
gene_expression_input(X)
}
\arguments{
  \item{X}{
gene expression data matrix with row representing genes/probe sets, and column representing samples.
}
}

\value{
If the input is valid, the output will be the same as the input; otherwise,
if the input contains NA, the corresponding rows will be deleted. if the input
contains negative value, the algorithm will stop and give error information.
}
\author{
 Niya Wang (wangny@vt.edu)
}

\examples{
gene_expression <- matrix(runif(2000),1000,2)
valid_gene_expression <- gene_expression_input(gene_expression)
}

