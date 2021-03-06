% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PM.R
\name{PM}
\alias{PM}
\title{ORFS method}
\usage{
PM(pik, EPS = 1e-11, comment = TRUE)
}
\arguments{
\item{pik}{a matrix of temporal inclusion probabilities.
Columns of \code{pik} correspond to the waves, and rows correspond to the units.
Inclusion probabilities can be totally unequal.}

\item{EPS}{a tolerance parameter. Default value is 1e-6.}

\item{comment}{a comment is written during the execution if \code{comment} is TRUE (default value).}
}
\value{
a matrix that contains temporal samples.
This is the update of \code{pik} and contains only 0s and 1s that indicates if a unit is selected or not at each wave.
}
\description{
Select temporal samples with fixed size at each wave.
It provides optimal time rotation of the selected samples using the systematic sampling method.
}
\references{
Quenouille, M. H. (1949). Approximate Tests of Correlation in time-Series. Royal Statistical Society, Series B Vol. 11, No. 1 (1949), pp. 68-84.
}
\seealso{
\code{\link{SystematicDesign}}
}
\author{
Esther Eustache \email{esther.eustache@unine.ch}
}
