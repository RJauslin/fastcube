#include <RcppArmadillo.h>



#include <RcppArmadillo.h>
#include "reduxArma.h"

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
bool isEye(arma::mat& M){
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


/*** R
isEye(diag(rep(1,10)))
isEye(diag(rep(1,10)) - diag(rep(1e-15,10)))

D <- diag(rep(1,10000))
D <- D[,1:5000]
system.time(isEye(D))
*/


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
void rrefArma(arma::mat& M){
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


/*** R
set.seed(1)
rm(list = ls())
N = 50
n = 30
p = 20
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A <- as.matrix(X/pik)

test <- t(A[1:(p+2),])
rrefArma(test)
test

*/


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
arma::vec osffphase(arma::vec prob, arma::mat Bm){
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
//' @title title
//'
//' @description
//' description
//'
//'
//' @param B matrix of auxiliary variables.
//' @param pik vector of inclusion probabilities
//' @param EPS tolerance
//'
//' @return updated pik
//'
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::vec onestep(arma::mat B,arma::vec pik,double EPS=0.0000001){

  arma::mat kern = arma::null(B);
  arma::uword N = pik.size();
  arma::vec u(N);
  u = kern.col(0);
  // int ncol = kern.n_cols;
  double l1 = 1e+200;
  double l2 = 1e+200;
  double l = 1e-9;

  for(arma::uword k = 0; k < N; k++){
    if(u[k]> 0){
      l1 = std::min(l1,(1.0 - pik[k])/u[k]);
      l2 = std::min(l2,pik[k]/u[k]);
    }
    if(u[k]< 0){
      l1 = std::min(l1,-pik[k]/u[k]);
      l2 = std::min(l2,(pik[k]-1.0)/u[k]);
    }
  }
  if(Rcpp::runif(1)[0]<l2/(l1+l2)){
    l = l1;
  }else{
    l = -l2;
  }
  for(arma::uword k = 0; k < N; k++){
    pik[k] = pik[k] + l*u[k];
    if(pik[k] < EPS){
      pik[k] = 0;
    }
    if(pik[k] > (1-EPS)){
      pik[k] = 1;
    }
  }
  return(pik);
}

// [[Rcpp::depends(RcppArmadillo)]]
//' @title title
//'
//'
//'
//' @param X matrix of auxiliary variables.
//' @param pik vector of inclusion probabilities
//' @param EPS tolerance
//'
//' @return a sample
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::vec flightphase_arma(arma::mat X,
                           arma::vec pik,
                           bool redux){
  double EPS = 1e-10;
  arma::mat D = arma::diagmat(1/pik);
  arma::mat A = D*X;
  unsigned int J = X.n_cols; // initial size of B

  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first"); // find first index of B
  arma::mat B = (A.rows(i)).t(); // extract B of A

  while(i.size() > 0){
    // std::cout << i.size() << std::endl;


    // if(i.size() < (J+1)){
    //   arma::mat kern = arma::null(B);
    //   if(kern.empty()){
    //     break;
    //   }
    // }

    // rrefArma(B);
    if(redux == true){
      Rcpp::List L = reduxArma(B.t());

      arma::mat B_tmp = L[0];
      arma::uvec ind_row = L[2];

      B_tmp = B_tmp.t();

      if(i.size() < (J+1)){
        arma::mat kern = arma::null(B_tmp);
        if(kern.empty()){
          break;
        }
      }
      rrefArma(B_tmp);
      // std::cout << B_tmp << std::endl;
      pik.elem(i.elem(ind_row)) = osffphase(pik.elem(i.elem(ind_row)),B_tmp);

    }else{
      if(i.size() < (J+1)){
        arma::mat kern = arma::null(B);
        if(kern.empty()){
          break;
        }
      }
      rrefArma(B);
      pik.elem(i) = osffphase(pik.elem(i),B);
    }


    i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
    B = (A.rows(i)).t();
    // std::cout << (A.t()*pik).t() << std::endl;
  }
  return(pik);
}



/*** R
library(sampling)
rm(list = ls())
N = 500
n = 20
p = 10
pik=inclusionprobabilities(runif(N),n)
X=cbind(pik,matrix(rnorm(N*p),c(N,p)))
A=X/pik
system.time(test <- flightphase_arma(X,pik,redux = TRUE))
t(A)%*%test
t(A)%*%pik






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





rm(list = ls())
set.seed(1)
eps <- 1e-13
library(Matrix)
N <- 250
Pik <- matrix(c(sampling::inclusionprobabilities(runif(N),70),
                sampling::inclusionprobabilities(runif(N),50),
                sampling::inclusionprobabilities(runif(N),30)),ncol = 3)
X <- PM(Pik)$PM
pik <- PM(Pik)$P
dim(X)
order = 2
method = 1
comment = TRUE
EPS = 1e-11

system.time(test <- flightphase_arma(X,pik))
# t(X)%*%test
# t(X)%*%pik
system.time(pikstar <- fastflightcubeSPOT(X, pik))
# t(X)%*%pikstar
# t(X)%*%pik
system.time(pikstar2 <- sampling::fastflightcube(X, pik))
# t(X)%*%pikstar2
# t(X)%*%pik
system.time(test <- BalancedSampling::flightphase(pik,X))
system.time(test <- fast.flight.cube(X,pik))
t(X)%*%test
t(X)%*%pik

*/
