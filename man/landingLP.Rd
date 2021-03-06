% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/landing.R
\name{landingLP}
\alias{landingLP}
\title{Landing by linear programming}
\usage{
landingLP(X, pikstar, pik, Xcat)
}
\arguments{
\item{X}{matrix of auxiliary variables on which the sample must be balanced.}

\item{pikstar}{vector of updated inclusion probabilities by the flight phase.}

\item{pik}{vector of inclusion probabilities.}

\item{Xcat}{matrix of categorical variable (not in its disjunctive form).}
}
\description{
This function does the landing phase of the cube method using linear programming.
 It allows you to put some categorical variable using the \code{Xcat} variable.
}
\examples{
\dontrun{

rm(list = ls())
library(sampling)

N <- 100
n <- 10
Xcat <-as.matrix(data.frame(cat1 = rep(1:4,each = N/4),
                    cat2 = rep(1:5,each = N/5),
                    cat2 = rep(1:10,each = N/10)))

p <- 3
X <- cbind(rep(1,N),matrix(rnorm(N*p),ncol = 3))
Xcat_tmp <- disjMatrix(Xcat)
Xred <- as.matrix(cbind(X,Xcat_tmp))

pik <- rep(n/N,N)
A <- Xred/pik

pikstar <- fastcube(X,pik,Xcat)
s <- landingLP(X,pikstar,pik,Xcat)
sum(s)
}
}
\references{
Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
}
\author{
Raphaël Jauslin \email{raphael.jauslin@unine.ch}
}
