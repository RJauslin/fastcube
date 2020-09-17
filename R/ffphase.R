#' @title Cube method with reduction of the auxiliary variables matrix
#'
#' @description Modified cube method.
#' This function reduces considerably the execution time when the matrix of auxiliary variables \code{X} contains lot of 0s.
#' It is based on the function \code{\link[sampling:samplecube]{samplecube}} from the package \code{sampling}.
#'
#'
#' @param X a matrix of size (N x p) of auxiliary variables on which the sample must be balanced.
#' @param pik a vector of size N of inclusion probabilities.
#'
#' @details In case where the number of auxiliary variables is great (i.e. p very large), even if we use the fast implementation proposed by
#' (Chauvet and Tille 2005), the problem is time consuming.
#' This function reduces considerably the execution time when the matrix of auxiliary variables \code{X} contains lot of 0s.
#' It considers a reduced matrix \code{X} by removing columns and rows that sum to 0 in the flight phase of the method (see  \code{\link{ReducedMatrix}} and \code{\link{ReducedFlightphase}}).
#'
#'
#' @return the updated vector of \code{pik} that contains only 0s and 1s that indicates if a unit is selected or not at each wave.
#'
#'
#' @author Esther Eustache \email{esther.eustache@@unine.ch}
#'
#'
#' @references
#' Chauvet, G. and Tille, Y. (2006). A fast algorithm of balanced sampling. Computational Statistics, 21/1:53-62
#'
#'
#' @seealso \code{\link[sampling:samplecube]{samplecube}}, \code{\link[sampling:landingcube]{landingcube}}, \code{\link{ReducedFlightphase}}, \code{\link{ReducedMatrix}}
#'
#'
#' @export
#' @examples
#' \dontrun{
#'
#'
#' ########################################################
#'
#'
#' rm(list = ls())
#'
#' library(sampling)
#'
#' N <- 1000
#'
#' Xcat <-as.matrix(data.frame(cat1 = rep(1:40,each = N/40),
#'                     cat2 = rep(1:50,each = N/50),
#'                     cat2 = rep(1:100,each = N/100)))
#'
#'
#' p <- 30
#' X <- matrix(rnorm(N*p),ncol = 30)
#'
#'
#' Xcat_tmp <- disjMatrix(Xcat)
#' Xcat_tmp <- do.call(cbind,apply(Xcat,MARGIN = 2,disjunctive))
#' Xred <- as.matrix(cbind(X,Xcat_tmp))
#'
#' pik <- rep(300/N,N)
#' A <- Xred/pik
#'
#'
#' system.time(s1 <- fastcube(Xred,pik))
#' system.time(s1 <- fastcube(X,pik,Xcat))
#' system.time(s1 <- fastcube(X,pik))
#' as.vector(t(A)%*%s1)
#' as.vector(t(A)%*%pik)
#'
#' system.time(s2 <- ReducedSamplecube(Xred,pik,redux = TRUE))
#' A <- Xred/pik
#' as.vector(t(A)%*%s2)
#' as.vector(t(A)%*%pik)
#'
#' system.time(s3 <- fastflightcube(Xred,pik,order = 2,comment = FALSE))
#' as.vector(t(A)%*%s3)
#' as.vector(t(A)%*%pik)
#'
#'
#'
#' system.time(s4 <- BalancedSampling::flightphase(pik,Xred))
#' as.vector(t(A)%*%s4)
#' as.vector(t(A)%*%pik)
#'
#' }
ffphase <- function(X, pik, Xcat){

  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------


  EPS = 1e-8
  A = X/pik
  p = ncol(X)
  if(!missing(Xcat)){
    CAT = TRUE
    n_all_cat <- sum(ncat(Xcat))
  }else{
    CAT = FALSE
    n_all_cat <- 0
  }



  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------

  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)

  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------

  while(i_size > 0){


    ##-----------------------------------
    ##  Find B
    ##-----------------------------------


    if(i_size <= p + n_all_cat){
      if(CAT == TRUE){
        Xdev <- disjMatrix(Xcat[i,])
        B <- cbind(A[i,],Xdev)
      }else{
        B <- A[i,]
      }
    }else{
      if(CAT == TRUE){
        B <- findB(A[i,],Xcat[i,])
      }else{
        B <- A[i[1:(p+1)],]
      }
    }
    tmp <-  onestep(B,pik[i[1:nrow(B)]],EPS)
    if(is.null(tmp)){
      break;
    }else{
      pik[i[1:nrow(B)]] <- tmp
    }


    ##------------
    ##  update i
    ##------------

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)


    # ##-----------------------------------
    # ##  Depending if we have enough row
    # ##-----------------------------------
    #
    # if(i_size <= p + n_all_cat){
    #
    #   if(CAT == TRUE){
    #     Xdev <- disjMatrix(Xcat[i,])
    #     B <- cbind(A[i,],Xdev)
    #   }else{
    #     B <- A[i,]
    #   }
    #
    #
    #
    #   # if(i_size > EPS){
    #   #   kern <- MASS::Null(B)
    #   #   if(length(kern)==0){
    #   #     break;
    #   #   }
    #   # }
    #   # pik[i] <- onestep(B,pik[i],EPS)
    #
    #   tmp <-  onestep(B,pik[i],EPS)
    #   if(is.null(tmp)){
    #     break;
    #   }else{
    #     pik[i] <- tmp
    #   }
    #
    # }else{
    #
    #   if(CAT == TRUE){
    #     B <- findB(A[i,],Xcat[i,])
    #   }else{
    #     B <- A[i[1:(p+1)],]
    #   }
    #
    #   tmp <-  onestep(B,pik[i[1:nrow(B)]],EPS)
    #   if(is.null(tmp)){
    #     break;
    #   }else{
    #     pik[i[1:nrow(B)]] <- tmp
    #   }
    #
    #   # pik[i[1:nrow(B)]] <- onestep(B,pik[i[1:nrow(B)]],EPS)
    # }
    #
    # ##------------
    # ##  update i
    # ##------------
    #
    # i <- which(pik > EPS & pik < (1-EPS))
    # i_size = length(i)

  }
  return(pik)
}








