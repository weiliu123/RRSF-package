\name{mRNA_matrix}
\alias{mRNA_matrix}
\docType{data}
\title{
The expression data
}
\description{
An example of GBM expression data with 200 samples and 300 genes. We acknowledge the TCGA Research Network for generating the GBM datasets.
}
\usage{data("mRNA_matrix")}
\format{
  The format is:
 num [1:200, 1:300] -0.929 2.12 -1.59 0.112 1.86 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:200] "TCGA-12-1597" "TCGA-02-0033" "TCGA-76-4929" "TCGA-06-0648" ...
  ..$ : chr [1:300] "5806" "3588" "3627" "55240" ...
}
\examples{
data(mRNA_matrix)
}
\keyword{datasets}
