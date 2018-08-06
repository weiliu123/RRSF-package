\name{RRSF-package}
\alias{RRSF-package}
\alias{RRSF}
\docType{package}
\title{
\packageTitle{RRSF}
}
\description{
\packageDescription{RRSF}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{RRSF}
\packageIndices{RRSF}
Very simple to use. Accepts x,y data for RRSF models, and makes predictions for new samples. 
Only 2 functions:
\code{RRSF} A rewighted random survival forest algorithm for survival prediction and gene selection.
\code{predict.RRSF} predicts the risk of new patients using the model fitted by \code{RRSF}.
}
\author{
\packageAuthor{RRSF}

Maintainer: \packageMaintainer{RRSF}
}
\references{
A rewighted random survival forest algorithm for survival prediction and gene selection.
}
\keyword{ package }
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
lp <- predict.RRSF(res$RRSFmodel,testSmpl)

}
