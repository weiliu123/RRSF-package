\name{RRSF}
\alias{RRSF}
\title{
Reweighted random survival forest
}
\description{
A rewighted random survival forest algorithm for survival prediction and gene selection.
}
\usage{
RRSF(x, y, geneSel = FALSE, DEBUG = FALSE, standardize = TRUE, globalGraph = NULL, topo.Wt = NULL, Gamma = 0.3, ...)
}
\arguments{
  \item{x}{
a n x p matrix of gene expression measurements with n samples and p genes.
}
  \item{y}{
a n x 2 matrix of survival data. The two columns represent disease status 'status' and survival time 'time' respectively. 
}
  \item{geneSel}{
Logical. Should gene selection be performed.
}
  \item{DEBUG}{
Logical. Should debugging information be shown.
}
  \item{standardize}{
Logical flag for \code{x} standardization, prior to fitting the model. Default is \code{TRUE}.
}
  \item{globalGraph}{
An \code{igraph} R object containing the global pathway graph.
}
  \item{topo.Wt}{
The topological weights of genes.
}
  \item{Gamma}{
A numeric value. The restart probability in directed random walk. Default is \code{Gamma = 0.3}.
}
  \item{\dots}{
Arguments to be passed to \code{rfsrc} in R package \code{randomForestSRC}.
}
}
\details{
RRSF integrates gene interaction information into the random survival forest algorithm for accurate survial prediction and robust gene selection.
}
\value{
Fitted \code{"RRSF"} model object.
  \item{RRSFmodel }{An object of class "rfsrc"}
  \item{topo.Wt}{The topological weights of genes}
  \item{geneSelected}{If \code{geneSel = TRUE}, A matrix recording the genes selected at each level, else \code{NULL} }
}
\references{
Liu, W., et al., Integration of gene interaction information into a reweighted random survival forest approach for accurate survival prediction and survival biomarker discovery. 
  To be published.
}
\author{
Wei Liu
}

\seealso{
 \code{\link{predict.RRSF}}
}
\examples{
data(dGMMirGraph)
data(mRNA_matrix)
data(survData)

trainSmpl.Idx <- sample(1:dim(mRNA_matrix)[1], floor(1/2*dim(mRNA_matrix)[1]))
testSmpl.Idx <- setdiff(1:dim(mRNA_matrix)[1], trainSmpl.Idx)
trainSmpl <- mRNA_matrix[trainSmpl.Idx ,]
testSmpl <- mRNA_matrix[testSmpl.Idx ,]

res <- RRSF(x=trainSmpl, y=survData[trainSmpl.Idx ,], geneSel = TRUE, DEBUG=TRUE, standardize=TRUE, 
            globalGraph = dGMMirGraph, topo.Wt = NULL, Gamma=0.3, ntree=100)

}
