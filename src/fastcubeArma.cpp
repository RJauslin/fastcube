#include <RcppArmadillo.h>
#include "reduxArma.h"
#include "disj.h"

using namespace Rcpp;


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
arma::vec onestep2(arma::mat B,arma::vec pik,double EPS=0.0000001){

  arma::mat kern = arma::null(B.t());
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
//' @title fastcubeArma
//'
//' @param X matrix
//' @param Xcat matrix
//' @param pik matrix
//'
//' @return sample
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::vec fastcubeArma(arma::mat X,
                       arma::umat Xcat,
                       arma::vec pik){

  double eps = 1e-9;
  arma::mat D = arma::diagmat(1/pik);
  arma::mat A = D*X;

  int p = X.n_cols;
  int n_all_cat = sum(ncat(Xcat));

  arma::mat B;
  // arma::uvec i = arma::find(pik > eps && pik < (1-eps), (p+n_all_cat)+1, "first");
  arma::uvec i = arma::find(pik > eps && pik < (1-eps));
  int i_size = i.size();
  // arma::mat B = findBarma(X.rows(i),Xcat.rows(i));

  while(i_size > 0){
    std::cout << i_size << std::endl;

    if(i_size <= p + n_all_cat){
      arma::umat Xdev = disjMatrix(Xcat.rows(i));
      arma::mat Xdev_tmp = arma::conv_to<arma::mat>::from(Xdev);
      B = arma::join_rows(A.rows(i),Xdev_tmp);
      if(i_size > eps){
        arma::mat kern = arma::null(B.t());
        if(kern.empty()){
          break;
        }
      }
      pik.elem(i) = onestep2(B,pik.elem(i));

    }else{
       B = findBarma(X.rows(i),Xcat.rows(i));
      pik.elem(i.head(B.n_rows)) = onestep2(B,pik.elem(i.head(B.n_rows)));
    }


    i = arma::find(pik > eps && pik < (1-eps));
    i_size = i.size();


  }
  return(pik);
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


Xcat_tmp <- disjMatrix(Xcat)
Xcat_tmp <- do.call(cbind,apply(Xcat,MARGIN = 2,disjunctive))
Xred <- as.matrix(cbind(X,Xcat_tmp))

pik <- rep(300/N,N)
A <- Xred/pik
system.time(s1 <- fastcubeArma(X,Xcat,pik))
system.time(s1 <- fastcube(X,Xcat,pik))
as.vector(t(A)%*%as.vector(s1))
as.vector(t(A)%*%pik)
*/

