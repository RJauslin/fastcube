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
ReducedSamplecube <- function(X, pik, redux = TRUE, t){

  ##----------------------------------------------------------------
  ##                        Initialization                         -
  ##----------------------------------------------------------------


  EPS = 1e-8
  A = X/pik
  J = ncol(X)


  ##----------------------------------------------------------------
  ##                Number of non 0-1 inclusion prob               -
  ##----------------------------------------------------------------




  i <- which(pik > EPS & pik < (1-EPS))
  i_size = length(i)
  if(i_size >= (J+1)){
    i <- i[1:(J+1)]
    B = A[i,]
  }else{
    B = A[i,]
  }

  # i <- which(pik > EPS & pik < (1-EPS))[1:(J+1)]
  # i_size = length(i)



  ##---------------------------------------------------------------
  ##                          Main loop                           -
  ##---------------------------------------------------------------

  while(i_size > 0){


    ##-----------------------
    ##  if redux is desired
    ##-----------------------

    if(redux == TRUE){
      pik_tmp <- pik[i]
      tmp <- ReducedMatrix(B)

      B_tmp <- tmp$B
      pik_tmp[tmp$ind_row]<- onestep(B_tmp, pik_tmp[tmp$ind_row],EPS)
      pik[i] <- pik_tmp

    }else{
      pik[i] <- onestep(B,pik[i],EPS)
    }


    ##------------
    ##  update i
    ##------------

    i <- which(pik > EPS & pik < (1-EPS))
    i_size = length(i)



    ##-----------------------------------
    ##  Depending if we have enough row
    ##-----------------------------------


    if(i_size >= (J+1)){
      i <- i[1:(J+1)]
      B = A[i,]
    }else if(i_size > t){
      B = A[i,]
    }else{

      B = A[i,]

      if(i_size > EPS){
        kern <- MASS::Null(B)
        if(length(kern)==0){
          break;
        }
      }
    }
  }


  return(pik)
}


