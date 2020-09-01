#' Title
#'
#' @param Xaux
#' @param Xcat
#'
#' @return
#' @export
#'
#' @examples
#' rm(list = ls())
#'
#' library(sampling)
#'
#' N <- 1000
#'
#' Xcat <-data.frame(cat1 = rep(1:40,each = N/40),
#'                     cat2 = rep(1:50,each = N/50),
#'                     cat3 = rep(1:500,times = 2))
#'                     #cat3 = rep(1:100,each = N/100))
#'
#' pik <- inclusionprobastrata(Xcat[,1],rep(1,40))
#' p <- 30
#' X <- matrix(rnorm(N*p),ncol = 30)
findB <- function(X,
                  Xcat){
  eps <- 1e-9
  N <- nrow(X)
  pInit <- ncol(X)

  Xcat_tmp <- Xcat[1:(pInit+1),]
  ncat <- sum(apply(Xcat_tmp,MARGIN = 2,FUN <- function(x){nlevels(as.factor(x))}))
  ncat_tmp <- 0

  while(ncat != ncat_tmp){
    ncat_tmp = ncat
    p =  pInit  + ncat
    Xcat_tmp <- Xcat[1:(p+1),]
    ncat <- sum(apply(Xcat_tmp,MARGIN = 2,FUN <- function(x){nlevels(as.factor(x))}))
  }

  Xdev <- as.matrix(do.call(cbind,apply(Xcat_tmp,MARGIN = 2,FUN <- function(x){
    if(all(x == 1)){
     return(as.matrix(x))
    }else{
      return(as.matrix(model.matrix(~as.factor(x)-1)))
    }
    })))
  return(cbind(X[1:(p+1),],Xdev))


}
