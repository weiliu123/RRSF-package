\name{DRW}
\alias{DRW}
\title{
Directed Random Walk
}
\description{
The directed random walk algorithm proposed by Liu et al (2013).
}
\usage{
DRW(igraphM, p0, EdgeWeight = FALSE, gamma = 0.3)
}
\arguments{
  \item{igraphM}{An igraph object containing the directed global pathway graph.}
  \item{p0}{A unit vector containing the initial weights of genes in the global pathway graph.}
  \item{EdgeWeight}{Logical. Should \code{igraphM} be converted to a weighted matrix or an un-weighted matrix (the default)?}
  \item{gamma}{A numeric value. The restart probability in directed random walk.}
}
\details{
  This function implements the directed random walk algorithm proposed by Liu et al (2013). It evaluates the topological
  weight of each gene according to its topological importance in the global pathway graph. The genes that close to many 
  other genes that have large initial weights will receive larger weights. The final weights reflect the topological importances 
  of genes in the global pathway graph. 
}
\value{
A numerical vector containing the topological weights of nodes in \code{igraphM}.
}
\references{
  Liu, W., et al., Topologically inferring risk-active pathways toward precise cancer classification by directed random walk. 
  Bioinformatics, 2013. 29(17): p. 2169-77.
}
\author{
Wei Liu <freelw@qq.com>
}

\examples{
  data(dGMMirGraph)
  vertexs <- V(dGMMirGraph)
  p0 <- runif(length(vertexs), min = 0, max = 1)
  names(p0) <- vertexs$name
  p0 <- p0/sum(p0)
  vertexWeight <- DRW(igraphM = dGMMirGraph, p0, EdgeWeight=FALSE, gamma = 0.3)
  names(vertexWeight) <- names(p0)
}
