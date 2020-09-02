#include <RcppArmadillo.h>



// [[Rcpp::depends(RcppArmadillo)]]
//' @title disj
//'
//' @description
//'  disjunctive
//'
//' @param strata vector of integer
//'
//' @return vector
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::umat disj(arma::uvec strata,int ncat) {
  int N = strata.size();
  // arma::uvec cat = arma::unique(strata);
  // int ncat = cat.size();

  // arma::uvec val = arma::regspace<arma::uvec>(0,1,ncat-1); // 0 to unique(strata)
  // arma::uvec strata_tmp = strata;
  // for(arma::uword i = 0;i<ncat;i++){
  //   strata_tmp.replace(cat[i],val[i]);
  // }

  arma::umat m(N,ncat,arma::fill::zeros);
  for(arma::uword i = 0;i < N;i++){
    m(i,strata(i)-1) = 1;
  }
  return(m);
}






/***R
strata=c(-2,3,-2,3,4,4,4,-2,-2,3,4,0,0,0)
disj(strata)


rm(list = ls())
set.seed(1)
strata <- sample(x = 1:6, size = 50, replace = TRUE)
system.time(M <- disj(strata,6))
system.time(M <- model.matrix(~as.factor(strata)-1))


rm(list = ls())
N <- 100000
strata <- sample(x = 1:400, size = N, replace = TRUE)
system.time(M <- disj(strata,400))
system.time(M <- model.matrix(~as.factor(strata)-1))
system.time(M <- sampling::disjunctive(strata))



*/





// [[Rcpp::depends(RcppArmadillo)]]
//' @title ncat
//'
//' @description
//' number of cat in each column
//'
//' @param Xcat Matrix
//'
//' @return vector
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::rowvec ncat(arma::umat Xcat) {
  int p = Xcat.n_cols;
  arma::rowvec ncat(p,arma::fill::zeros);
  for(arma::uword i = 0; i < p;i++){
    arma::uvec cat = unique(Xcat.col(i));
    int tmp = cat.size();
    ncat(i) = tmp;
  }
  return(ncat);
}

/***R
rm(list = ls())
Xcat <-  sample(x = 1:6, size = 100, replace = TRUE)
for(i in 1:1000){
  Xcat <- cbind(Xcat, sample(x = 1:6, size = 100, replace = TRUE))
}

system.time(n <- apply(Xcat,MARGIN = 2,FUN <- function(x){nlevels(as.factor(x))}))
system.time(test <- ncat(Xcat))

*/



// [[Rcpp::depends(RcppArmadillo)]]
//' @title disjMatrix
//'
//' @description
//' disj on each column
//'
//' @param starta a Matrix
//'
//' @return a matrix
//'
//' @author Raphaël Jauslin \email{raphael.jauslin@@unine.ch}
//'
//' @export
// [[Rcpp::export]]
arma::umat disjMatrix(arma::umat strata) {

  int N = strata.n_rows;
  int p = strata.n_cols;
  arma::rowvec all_cat = ncat(strata);
  int n_all_cat = sum(all_cat);

  arma::umat m(N,n_all_cat,arma::fill::zeros);
  arma::rowvec subind = arma::round(cumsum(all_cat)-1);
  arma::uvec ind = arma::conv_to<arma::uvec>::from(subind);


  for(arma::uword i = 0;i < (p-1);i++){
    arma::uvec tmp = strata.col(i);
    arma::umat tmp_mat(N,all_cat(i),arma::fill::zeros);
    tmp_mat = disj(tmp,all_cat(i));

    if(i == 0){
      m.cols(0, ind(i)) = tmp_mat;
    }else{
      m.cols(ind(i-1)+1, ind(i)) = tmp_mat;
    }

  }

  return(m);
}


/***R
rm(list = ls())

N <- 1000
Xcat <- sample(x = 1:40, size = N, replace = TRUE)
for(i in 1:1000){
  Xcat <- cbind(Xcat,sample(x = 1:40, size = N, replace = TRUE))
}
# sum(ncat(Xcat))
# ncat(Xcat)


system.time(test <- apply(as.matrix(Xcat),MARGIN = 2,FUN <- function(x){as.matrix(model.matrix(~as.factor(x)-1))}))
system.time(test <-  disjMatrix(Xcat))

*/
