




#' Internal function of ReducedSamplecube
#' @noRd
onestep <- function(B,pik,EPS){

  kern <- MASS::Null(B)
  N <- length(pik)
  u = kern[,1]

  l1=min(pmax((1-pik)/u,-pik/u))
  l2=min(pmax((pik-1)/u,pik/u))

  if(runif(1) < l2/(l1+l2)){
    l = l1;
  }else{
    l = -l2;
  }
  pik = pik + l*u

  return(pik);
}


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
#' pik <- rep(300/N,N)
#' system.time(s1 <- fastcube(X,Xcat,pik))
#' A <- X/pik
#' t(A)%*%s1
#' t(A)%*%pik
#'
#'
#' Xcat_tmp <- disjMatrix(Xcat)
#' Xred <- as.matrix(cbind(X,Xcat_tmp))
#' system.time(s2 <- ReducedSamplecube(Xred,pik,t = 1000))
#' system.time(s2 <- cube(Xred,pik,t = 1000))
#' A <- Xred/pik
#' t(A)%*%s2
#' t(A)%*%s1
#' t(A)%*%pik
#'
#'
#' }
fastcube <- function(X, Xcat, pik){



  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------


  EPS = 1e-8
  A = X/pik
  p = ncol(X)
  # ncat <- sum(apply(Xcat,MARGIN = 2,FUN <- function(x){nlevels(as.factor(x))}))
  n_all_cat <- sum(ncat(Xcat))


  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------

  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)
  B <- findB(A[i,],Xcat[i,])

  # i <- which(pik > EPS & pik < (1-EPS))[1:(J+1)]
  # i_size = length(i)



  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------

  while(i_size > 0){
    print(dim(B))
    # print(i_size)
    pik[i[1:nrow(B)]] <- onestep(B,pik[i[1:nrow(B)]],EPS)



    ##------------
    ##  update i
    ##------------

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)



    ##-----------------------------------
    ##  Depending if we have enough row
    ##-----------------------------------

    if(i_size <= p + n_all_cat){

      # if(ncol(as.matrix(Xcat[i,])) == 1){
      #   Xdev <- model.matrix(~as.factor(as.matrix(Xcat[i,]))-1)
      # }else{
      #   Xdev <- as.matrix(do.call(cbind,apply(as.matrix(Xcat[i,]),MARGIN = 2,FUN <- function(x){as.matrix(model.matrix(~as.factor(x)-1))})))
      # }

      Xdev <- disjMatrix(Xcat[i,])
      B <- cbind(A[i,],Xdev)
      if(i_size > EPS){
        kern <- MASS::Null(B)
        if(length(kern)==0){
          break;
        }
      }else{
        pik[i] <- onestep(B,pik[i],EPS)
      }



    }else{
      B <- findB(A[i,],Xcat[i,])
    }

  }


  return(pik)
}








