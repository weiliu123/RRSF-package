\name{predict.RRSF}
\alias{predict.RRSF}
\title{
Make predictions from a "RRSF" object
}
\description{
This function predicts the risk of new samples from a fitted \code{RRSF} object.
}
\usage{
predict.RRSF(object, newx, ...)
}
\arguments{
  \item{object}{
Fitted \code{"rfsrc"} model object.
}
  \item{newx}{
A matrix with new samples to predict.
}
  \item{\dots}{
Arguments to be passed to \code{predict.rfsrc} in R package \code{randomForestSRC}.
}
}
\value{
Predicted results of new patients in \code{newx}.
}

\seealso{
\code{\link{RRSF}}
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
lp <- predict.RRSF(res$RRSFmodel,testSmpl)

}
