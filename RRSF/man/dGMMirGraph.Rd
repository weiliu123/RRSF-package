\name{dGMMirGraph}
\alias{dGMMirGraph}
\docType{data}
\title{
The global pathway graph
}
\description{
The directed global pathway graph constructed by the R package \code{iSubpathwayMiner}. 
}
\usage{data("dGMMirGraph")}
\format{
An \code{igraph} R object.
}
\details{
There are 7159 nodes and 39930 edges in \code{dGMMirGraph}. Each node in the graph represents a gene/miRNA/metabolite.
The global pathway graph is used to evaluate the topological importance of genes by directed random walk.
}
\examples{
data(dGMMirGraph)
}
\keyword{datasets}
