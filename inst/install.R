library(devtools)
install_github("Rjauslin/fastcube@master")
library(fastcube)

rm(list = ls())
# set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 500
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),50),
                sampling::inclusionprobabilities(runif(N),100),
                sampling::inclusionprobabilities(runif(N),150)),ncol = 3)
X <- PM(Pik)$PM
pik <- PM(Pik)$P
image(as(X,"sparseMatrix"))
# rrefArma(X)
dim(X)
order = 2
EPS = 1e-11
A=X/pik
system.time(test <- flightphase_arma(X,pik,redux = TRUE))
system.time(test1 <- ReducedSamplecube(X,pik,t = 3))
t(A)%*%test1
t(A)%*%pik
