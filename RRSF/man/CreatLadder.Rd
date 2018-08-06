\name{CreatLadder}
\alias{CreatLadder}
\title{
Create a decreasing sequence
}
\description{
This function creates a decreasing sequence according to a specific decreasing ratio.
}
\usage{
CreatLadder(Ntotal, pRatio = 0.75, Nmin = 5)
}
\arguments{
  \item{Ntotal}{
Total number
}
  \item{pRatio}{
Decreasing ratio
}
  \item{Nmin}{
Minimum number
}
}
\value{
A decreasing sequence.
}
\examples{
ladder <- CreatLadder(Ntotal = 300, pRatio = 0.9, Nmin = 2)
}
