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
  Xcat <- as.matrix(Xcat)
  X <- as.matrix(X)
  eps <- 1e-9
  N <- nrow(X)
  pInit <- ncol(X)

  Xcat_tmp <- Xcat[1:(pInit+1),]
  # ncat <- sum(apply(as.matrix(Xcat_tmp),MARGIN = 2,FUN <- function(x){nlevels(as.factor(x))}))
  n_all_cat <- sum(ncat(Xcat_tmp))
  n_all_cat_tmp <- 0

  while(n_all_cat != n_all_cat_tmp){
    n_all_cat_tmp = n_all_cat
    p =  pInit  + n_all_cat
    Xcat_tmp <- Xcat[1:(p+1),]
    # ncat <- sum(apply(as.matrix(Xcat_tmp),MARGIN = 2,FUN <- function(x){nlevels(as.factor(x))}))
    n_all_cat <- sum(ncat(Xcat_tmp))
  }

  Xcat_tmp <- as.matrix(Xcat_tmp)
  Xdev <- disjMatrix(Xcat_tmp)
  # if(ncol(Xcat_tmp) == 1){
  #   if(all(Xcat_tmp[,1] == Xcat_tmp[1,1])){
  #     Xdev <- Xcat_tmp
  #   }else{
  #     Xdev <- model.matrix(~as.factor(Xcat_tmp)-1)
  #   }
  # }else{
  #   Xdev <- as.matrix(do.call(cbind,apply(Xcat_tmp,MARGIN = 2,FUN <- function(x){
  #     if(all(x == 1)){
  #       return(as.matrix(x))
  #     }else{
  #       return(as.matrix(model.matrix(~as.factor(x)-1)))
  #     }
  #   })))
  # }

  return(cbind(X[1:(p+1),],Xdev))


}
