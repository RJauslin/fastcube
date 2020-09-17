#' find the sub-matrix B in the flight phase
#'
#' @description
#' This function is used in the \code{\link{ffphase}} function to find the smallest matrix B.
#'
#' @param X matrix of auxiliary variables.
#' @param Xcat matrix of categorical variables.
#'
#' @details
#'
#' The function find the smallest matrix B such that is contains only one additional row.
#'  It consecutively add the right number of row depending of the number of categroies that is added.
#'
#' @return matrix B
#'
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
  n_all_cat <- sum(ncat(Xcat_tmp))
  n_all_cat_tmp <- 0

  while(n_all_cat != n_all_cat_tmp){
    n_all_cat_tmp = n_all_cat
    p =  pInit  + n_all_cat
    Xcat_tmp <- Xcat[1:(p+1),]
    n_all_cat <- sum(ncat(Xcat_tmp))
  }

  Xcat_tmp <- as.matrix(Xcat_tmp)
  Xdev <- disjMatrix(Xcat_tmp)

  return(cbind(X[1:(p+1),],Xdev))

}
