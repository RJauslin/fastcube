% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ffphase.R
\name{ffphase}
\alias{ffphase}
\title{Cube method with reduction of the auxiliary variables matrix}
\usage{
ffphase(X, pik, Xcat)
}
\arguments{
\item{X}{a matrix of size (N x p) of auxiliary variables on which the sample must be balanced.}

\item{pik}{a vector of size N of inclusion probabilities.}
}
\value{
the updated vector of \code{pik} that contains only 0s and 1s that indicates if a unit is selected or not at each wave.
}
\description{
Modified cube method.
This function reduces considerably the execution time when the matrix of auxiliary variables \code{X} contains lot of 0s.
It is based on the function \code{\link[sampling:samplecube]{samplecube}} from the package \code{sampling}.
}
\details{
In case where the number of auxiliary variables is great (i.e. p very large), even if we use the fast implementation proposed by
(Chauvet and Tille 2005), the problem is time consuming.
This function reduces considerably the execution time when the matrix of auxiliary variables \code{X} contains lot of 0s.
It considers a reduced matrix \code{X} by removing columns and rows that sum to 0 in the flight phase of the method (see  \code{\link{ReducedMatrix}} and \code{\link{ReducedFlightphase}}).
}
\examples{
\dontrun{


########################################################


rm(list = ls())

library(sampling)

N <- 1000

Xcat <-as.matrix(data.frame(cat1 = rep(1:40,each = N/40),
                    cat2 = rep(1:50,each = N/50),
                    cat2 = rep(1:100,each = N/100)))


p <- 30
X <- matrix(rnorm(N*p),ncol = 30)


Xcat_tmp <- disjMatrix(Xcat)
Xcat_tmp <- do.call(cbind,apply(Xcat,MARGIN = 2,disjunctive))
Xred <- as.matrix(cbind(X,Xcat_tmp))

pik <- rep(300/N,N)
A <- Xred/pik


system.time(s1 <- fastcube(Xred,pik))
system.time(s1 <- fastcube(X,pik,Xcat))
system.time(s1 <- fastcube(X,pik))
as.vector(t(A)\%*\%s1)
as.vector(t(A)\%*\%pik)

system.time(s2 <- ReducedSamplecube(Xred,pik,redux = TRUE))
A <- Xred/pik
as.vector(t(A)\%*\%s2)
as.vector(t(A)\%*\%pik)

system.time(s3 <- fastflightcube(Xred,pik,order = 2,comment = FALSE))
as.vector(t(A)\%*\%s3)
as.vector(t(A)\%*\%pik)



system.time(s4 <- BalancedSampling::flightphase(pik,Xred))
as.vector(t(A)\%*\%s4)
as.vector(t(A)\%*\%pik)

}
}
\references{
Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
}
\seealso{
\code{\link[sampling:samplecube]{samplecube}}, \code{\link[sampling:landingcube]{landingcube}}, \code{\link{ReducedFlightphase}}, \code{\link{ReducedMatrix}}
}
\author{
Esther Eustache \email{esther.eustache@unine.ch}
}
