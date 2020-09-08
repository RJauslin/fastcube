#include <RcppArmadillo.h>
#include "reduxArma.h"
#include "disj.h"

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
//' @title is Identiy matrix
//'
//' @param M matrix
//'
//' @return bool
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
bool isEye2(arma::mat& M){
  bool out = true;
  int p = M.n_cols;
  int N = M.n_rows;
  double eps = 1e-12;
  for(int i = 0;i< N; i++){
    for(int j = 0; j < p; j++){
      if(i == j && M(i,j) != 1){
        out = false;
        break;
      }
      if(i != j && M(i,j) != 0){
        out = false;
        break;
      }
    }
  }
  return out;
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @title reduced row echelon form arma implementation
//'
//'
//' @param M matrix
//'
//' @return NULL (transform matrix)
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
void rrefArma2(arma::mat& M){
  int lead = 0;
  int rowCount = M.n_rows;
  int columnCount = M.n_cols;
  double eps = 1e-11;
  int i,k;
  double temp;
  for(int r = 0; r < rowCount; r++){
    if(columnCount <= lead){
      return;
    }
    i = r;
    while(std::max(M(i,lead),-M(i,lead)) < eps ){
      M(i,lead) = 0.0;
      i = i + 1;
      if(i == rowCount){
        i = r;
        lead = lead + 1;
        if(columnCount == lead){
          return;
        }
      }
    }
    // swap rows i and r
    for(int k = 0; k < columnCount;k++){
      temp = M(i,k);
      M(i,k) = M(r,k);
      M(r,k) = temp;
    }
    // If M(r, lead) is not 0 divide row r by M(r, lead)
    if( M(r,lead) != 0.0 ){
      temp = M(r,lead);
      for(int k = 0; k < lead;k++){
        M(r,k) = 0.0;
      }
      for(int k = lead;k < columnCount;k++){
        M(r,k) = M(r,k)/temp;
      }
    }
    for(int i = 0;i < rowCount;i++){
      if( i != r ){
        // Subtract M(i, lead) multiplied by row r from row i
        temp = M(i,lead);
        for( k = 0;k < columnCount; k++){
          M(i,k) = M(i,k) - temp * M(r,k);
        }
      }
    }
    lead = lead + 1;
  }
  return;
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//' @description
//' description
//'
//'
//' @param prob inclusion probabilities
//' @param Bm matrix of auxiliary variables
//'
//' @details
//'
//' details
//'
//' @return a vector
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @seealso
//' func
//'
//' @examples
//'
//' @export
// [[Rcpp::export]]
arma::vec osffphase2(arma::vec prob, arma::mat Bm){
  int ncol = Bm.n_cols;
  int nrow = Bm.n_rows;

  arma::vec u(ncol,arma::fill::zeros);
  arma::uvec uset(ncol,arma::fill::zeros);

  double la1 = 1e+200;
  double la2 = 1e+200;
  double la, eps = 1e-9;
  int lead;
  double v, free = -1.0;
  // find nonzero vector u in Ker B (null space of B, i.e. Bu = 0)
  // with both positive and negative values
  // find reduced row echelon form of B
  // rrefArma(Bm);

  // std::cout << Bm << std::endl;
  for(int i = (nrow-1);i >= 0; i--){
    // find lead (first nonzero entry on row) if exists
    // if no lead, i.e lead = ncol, do nothing
    // if lead, the variables after are either set or free
    // free variables are alternately set to 1 or -1
    lead = 0;
    for(int j = 0; j < ncol; j++){
      if(Bm(i,j)==0.0){
        lead++;
      }else{
        break;
      }
    }
    // lead found
    if(lead<ncol){
      v = 0.0;
      for(int j = lead+1;j < ncol;j++){
        if( uset[j] == 0 ){
          uset[j] = 1;
          free *= -1.0;
          u[j] = free;
        }
        v -= u[j]*Bm(i,j);
      }
      u[lead] = v/Bm(i,lead);
      uset[lead] = 1;
    }
  }
  // unset u[i] are free and are set to 1 or -1, can only exist at beginning
  for(int i = 0;i < ncol;i++){
    if( uset[i] == 0 ){
      free *= -1.0;
      u[i] = free;
    }else{
      break;
    }
  }

  // find lambda1 and lambda2
  for(int i = 0;i < ncol;i++){
    if(u[i]>0){
      la1 = std::min(la1,(1-prob[i])/u[i]);
      la2 = std::min(la2,prob[i]/u[i]);
    }
    if(u[i]<0){
      la1 = std::min(la1,-prob[i]/u[i]);
      la2 = std::min(la2,(prob[i]-1)/u[i]);
    }
  }
  // random choice of p+lambda1*u and p-lambda2*u
  if(runif(1)[0]<la2/(la1+la2)){
    la = la1;
  }else{
    la = -la2;
  }
  // update prob
  for(int i = 0;i < ncol;i++){
    prob[i] = prob[i] + la * u[i];
    if(prob[i] < eps){ prob[i] = 0; }
    if(prob[i] > 1-eps){ prob[i] = 1; }
  }
  return prob;
}


// [[Rcpp::depends(RcppArmadillo)]]
//' @title Fast Flight phase
//'
//'
//' @description
//'
//' Modified version of \code{\link[BalancedSampling:flightphase]{flightphase}}
//'
//' @param prob vector of inclusion probabilities of size N.
//' @param Xbal Matrix of auxiliary variables of dimension N x p
//' @param order if reordering at first step, Default TRUE.
//' @param redux if the matrix should be reduced. Default FALSE.
//'
//' @details
//'
//' details
//'
//' @return a sample with at most p value not update to 0 or 1.
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::vec ffphase(arma::vec prob, arma::mat Xbal, bool redux = false){
  int N = prob.size();
  int naux = Xbal.n_cols;

  arma::uvec index(N);
  arma::vec p(N);

  int howmany;
  int k;
  for(int i = 0;i < N;i++){
    index[i]=i;
    p[i]=prob[i];
  }
  double eps = 1e-12;
  int tempInt, howlong;


  // put finished units at beginning of list
  int done = 0;
  for(int i = done;i < N;i++){
    if( p[index[i]]<eps || p[index[i]]>1-eps ){
      tempInt = index[done];
      index[done] = index[i];
      index[i] = tempInt;
      done = done + 1;
    }
  }


  // remaining are index from done to N-1
  while( done < N ){
    std::cout << done << std::endl;
    // find cluster of size howmany
    howmany = std::min(naux + 1,N-done);

    if( howmany > 1 ){

      arma::vec p_small(howmany);
      arma::vec dists(howmany); dists = 1e+20;
      arma::vec index_small(howmany);
      // arma::mat B(howmany-1,howmany);
      arma::mat B(naux,howmany);


      for(int i = 0;i < howmany; i++){
        index_small[i] = index[done+i];
        for(int j = 0;j < naux;j++){ // HERE WE CHANGED j < howmany by ----- > j < naux
          B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
        }
        p_small[i] = p[index_small[i]];
      }


      if(redux == true){
        Rcpp::List L = reduxArma(B.t());
        arma::mat B_tmp = L[0];
        B_tmp = B_tmp.t();
        arma::uvec ind_row = L[2];
        rrefArma2(B_tmp);

        if(howmany < naux + 1){

          bool test = isEye2(B_tmp);
          if(test == true){
            break;
          }else{
            p_small.elem(ind_row) = osffphase2(p_small.elem(ind_row),B_tmp);
          }
        }else{
          p_small.elem(ind_row) = osffphase2(p_small.elem(ind_row),B_tmp);
        }


      }else{
        rrefArma2(B);
        // std::cout << B << std::endl;
        if(howmany < naux + 1){

          bool test = isEye2(B);
          if(test == true){
            break;
          }else{
            p_small = osffphase2(p_small,B);
          }
        }else{

          std::cout << sum(p) << std::endl;
          p_small = osffphase2(p_small,B);
        }
      }



      // update prob
      for(int i = 0;i < howmany;i++){
        p[index_small[i]] = p_small[i];
      }
      // update done and index
      howlong = done + howmany;
      for(int i = done;i < howlong;i++){
        if( p[index[i]]<eps || p[index[i]]>1-eps ){
          tempInt = index[done];
          index[done] = index[i];
          index[i] = tempInt;
          done = done + 1;
        }
      }
    }else{
      // max one unit left
      if(runif(1)[0] < p[index[done]]){p[index[done]]=1;}else{p[index[done]]=0;}
      done = N;
    }

  }

  // round
  for(int i = 0;i < N;i++){
    if( p[index[i]] > 1-eps  ){
      p[index[i]] = 1;
    }
    if( p[index[i]] < eps  ){
      p[index[i]] = 0;
    }
  }
  return p;
}

/*** R
rm(list = ls())

library(sampling)
N <- 1000

Xcat <-as.matrix(data.frame(cat1 = rep(1:40,each = N/40),
                            cat2 = rep(1:50,each = N/50),
                            cat2 = rep(1:100,each = N/100)))

p <- 30
X <- matrix(rnorm(N*p),ncol = 30)



# rm(list = ls())
#
# library(sampling)
# N <- 100
#
# Xcat <-as.matrix(data.frame(cat1 = rep(1:4,each = N/4),
#                             cat2 = rep(1:5,each = N/5),
#                             cat2 = rep(1:10,each = N/10)))
#
# p <- 3
# X <- matrix(rnorm(N*p),ncol = 3)


Xcat_tmp <- disjMatrix(Xcat)
Xcat_tmp <- do.call(cbind,apply(Xcat,MARGIN = 2,disjunctive))
Xred <- as.matrix(cbind(X,Xcat_tmp))

pik <- rep(10/N,N)
A <- Xred/pik
system.time(s1 <- ffphase(pik,Xred,redux = FALSE))
system.time(s1 <- fastcube(X,Xcat,pik))
as.vector(t(A)%*%as.vector(s1))
as.vector(t(A)%*%pik)
*/




