#include <RcppArmadillo.h>


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
arma::vec onestep5(arma::mat B,arma::vec pik,double EPS=0.0000001){

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
arma::vec flightphase_arma5(arma::mat X,arma::vec pik,double EPS=0.0000001){

  arma::mat D = arma::diagmat(1/pik);
  arma::mat A = D*X;
  unsigned int J = X.n_cols; // initial size of B

  arma::uvec i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first"); // find first index of B
  arma::mat B = (A.rows(i)).t(); // extract B of A

  while(i.size() > 0){
    // std::cout << i.size() << std::endl;
    pik.elem(i) =  onestep5(B,pik.elem(i));
    i = arma::find(pik > EPS && pik < (1-EPS), J+1, "first");
    B = (A.rows(i)).t();
    // std::cout << B.n_cols << std::endl;
    // std::cout << B.n_rows << std::endl << std::endl;
    if(i.size() < (J+1)){
      arma::mat kern = arma::null(B);
      if(kern.empty()){
        break;
      }
    }
  }
  return(pik);
}

/*** R
rm(list = ls())
set.seed(1)
library(sampling)
N <- 100

Xcat <-as.matrix(data.frame(cat1 = rep(1:4,each = N/4),
                            cat2 = rep(1:5,each = N/5),
                            cat2 = rep(1:10,each = N/10)))

p <- 3
X <- matrix(rnorm(N*p),ncol = 3)


Xcat_tmp <- disjMatrix(Xcat)
# Xcat_tmp <- do.call(cbind,apply(Xcat,MARGIN = 2,disjunctive))
Xred <- as.matrix(cbind(X,Xcat_tmp))

pik <- rep(30/N,N)
A <- Xred/pik
system.time(s1 <- fastcubeArma(X,Xcat/pik,pik))
system.time(s1 <- fastcubeArma(X,Xcat/pik,pik))
sum(s1)
system.time(s1 <- fastcube(X,Xcat,pik))
as.vector(t(A)%*%as.vector(s1))
as.vector(t(A)%*%pik)
*/
