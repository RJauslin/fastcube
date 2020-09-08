#include <Rcpp.h>
#include "rref.h"
#include "onestepfastflightcube.h"

using namespace Rcpp;


// check is matrix is diagonal
bool isEye3(NumericMatrix M){
  bool out = true;

  int p = M.ncol();
  int N = M.nrow();
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



//' @title  flightphaseModif
//'
//' @param prob vector of inclusion probabilities of size N.
//' @param Xbal Matrix of auxiliary variables of dimension N x p
//'
//' @return a sample with at most p value not update to 0 or 1.
//'
//' @export
// [[Rcpp::export]]
NumericVector flightphaseModif(NumericVector prob, NumericMatrix Xbal){
  int N = prob.size();
  int naux = Xbal.ncol();

  IntegerVector index(N);
  NumericVector p(N);
  int i,j,k,howmany;
  for(i=0;i<N;i++){index[i]=i; p[i]=prob[i];}
  double eps = 1e-12;
  int done = 0, tempInt, howlong;
  // randomize order of index list
  NumericVector rnd = runif(N);
  for(i=0;i<N;i++){
    k = i + floor(rnd[i] * (N-i));
    tempInt = index[i];
    index[i] = index[k];
    index[k] = tempInt;
  }
  // put finished units at beginning of list
  for(i=done;i<N;i++){
    if( p[index[i]]<eps || p[index[i]]>1-eps ){
      tempInt = index[done];
      index[done] = index[i];
      index[i] = tempInt;
      done = done + 1;
    }
  }
  // remaining are index from done to N-1
  while( done < N ){

    std::cout << howmany << std::endl;
    // find cluster of size howmany
    howmany = std::min(naux+1,N-done);
    // stop if there are less than naux units left
    // if(howmany <= naux){done=N; break;} --------------------------------------------------- remove condition to break

    if( howmany > 1 ){
      NumericVector p_small(howmany);
      NumericVector dists(howmany,1e+200);
      IntegerVector index_small(howmany);
      // NumericMatrix B(howmany-1,howmany);
      NumericMatrix B(naux,howmany); //--------------------------------------------------- change to naux instead of howmany -1

      for(i=0;i<howmany;i++){
        index_small[i] = index[done+i];
        for(j=0; j < naux;j++){ //--------------------------------------------------- change to naux instead of howmany - 1
          B(j,i) = Xbal(index_small[i],j)/prob[index_small[i]];
        }
        p_small[i] = p[index_small[i]];
      }
      // p_small = onestepfastflightcube(p_small,B);

      rref(B); //--------------------------------------------------- change rref is done here instead of inside onestepfastflightcube
      if(howmany < naux + 1){
        bool test = isEye3(B);
        if(test == true){
          break;
        }else{
          p_small = onestepfastflightcube(p_small,B);
        }
      }else{
          p_small = onestepfastflightcube(p_small,B);
      }


      // update prob
      for(i=0;i<howmany;i++){
        p[index_small[i]] = p_small[i];
      }
      // update done and index
      howlong = done + howmany;
      for(i=done;i<howlong;i++){
        if( p[index[i]]<eps || p[index[i]]>1-eps ){
          tempInt = index[done];
          index[done] = index[i];
          index[i] = tempInt;
          done = done + 1;
        }
      }
    }else{
      // max one unit left
      if(runif(1)[0]<p[index[done]]){p[index[done]]=1;}else{p[index[done]]=0;}
      done = N;
    }
  }
  // round
  for(i=0;i<N;i++){
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
system.time(s1 <- flightphaseModif(pik,Xred))
system.time(s1 <- fastcube(X,Xcat,pik))
as.vector(t(A)%*%as.vector(s1))
as.vector(t(A)%*%pik)
*/
