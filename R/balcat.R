


#' Title
#'
#' @param X
#' @param Xcat
#' @param pik
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#' rm(list = ls())
#' N <- 100
#' pik <- rep(10/N,N)
#'
#' Xcat <-data.frame(cat1 = rep(1:50,each = 2))
#' Xcat_tmp <- disjMatrix(as.matrix(Xcat))
#' X <- cbind(pik,matrix(rnorm(N),ncol = 1))
#' A <- cbind(X,Xcat_tmp)/pik
#'
#' s <- balcat(X,Xcat,pik)
#' t(A)%*%s
#' t(A)%*%pik
#'
#'
#' rm(list = ls())
#' N <- 10000
#' pik <- rep(0.25,N)
#'
#' Xcat <- as.matrix(rep(1:1250,each = N/1250))
#' Xcat_tmp <- disjMatrix(Xcat)
#' X <- cbind(pik,matrix(rnorm(N),ncol = 1))
#' A <- cbind(X,Xcat_tmp)/pik
#'
#' system.time(s <- balcat(X,Xcat,pik))
#'
#'
#'
#' rm(list = ls())
#' N <- 1000
#' pik <- rep(0.25,N)
#'
#' Xcat <- as.matrix(rep(1:250,each = N/250))
#' Xcat_tmp <- disjMatrix(Xcat)
#' X <- cbind(pik,matrix(rnorm(N),ncol = 1))
#' A <- cbind(X,Xcat_tmp)/pik
#'
#'
#' system.time(s <- balcat(X,Xcat,pik))
#' system.time(s <- BalancedSampling::flightphase(pik,cbind(X,Xcat_tmp)))
#' system.time(s <- balancedstratification(X,Xcat,pik,comment=FALSE,method=2))
#'
#' microbenchmark( balcat(X,Xcat,pik),
#' BalancedSampling::flightphase(pik,cbind(X,Xcat_tmp)),
#' balancedstratification(X,Xcat,pik,comment=FALSE,method=2),times = 10)
#'
#'
balcat <- function(X,Xcat,pik){

  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------
  Xcat <- as.matrix(Xcat)
  EPS = 1e-8
  A = X/pik
  p = ncol(X)
  n_all_cat <- sum(ncat(Xcat))
  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------

  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)


  ##----------------------------------------------------------------
  ##            flight phase with categorical variable             -
  ##----------------------------------------------------------------


  while(i_size > p + n_all_cat){

    ##------ Remove unique category

    uCat <- i[duplicated(Xcat[i,]) | duplicated(Xcat[i,], fromLast = TRUE)]

    ##------ Find B

    B <- findB(as.matrix(A[uCat,]),as.matrix(Xcat[uCat,]))

    ##------ onestep and check if null

    tmp <-  onestep(B,pik[uCat[1:nrow(B)]],EPS)
    if(is.null(tmp)){
      break;
    }else{
      pik[uCat[1:nrow(B)]] <- tmp
    }

    ##------ update i

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)

  }

  # image(as(disjMatrix(as.matrix(Xcat[i,])),"sparseMatrix"))
  # head(t(AInit)%*%pik)
  # head(t(AInit)%*%pikInit)

  ##----------------------------------------------------------------
  ##            end of flight phase without categories             -
  ##----------------------------------------------------------------

  while(i_size > 0){


    ##------ Find B

    if(i_size <= p){
      B <- as.matrix(A[i,]) # if not enough row
    }else{
      B <- as.matrix(A[i[1:(p+1)],]) # A has p columns
    }

    ##------ onestep and check if null

    tmp <-  onestep(B,pik[i[1:nrow(B)]],EPS)
    if(is.null(tmp)){
      break;
    }else{
      pik[i[1:nrow(B)]] <- tmp
    }

    ##------ update i

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)

  }

  # head(t(AInit)%*%pik)
  # head(t(AInit)%*%pikInit)


  return(pik)
}
