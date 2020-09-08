#include <Rcpp.h>
using namespace Rcpp;


// one step fast flight cube
NumericVector onestepfastflightcube(NumericVector prob, NumericMatrix Bm){
  int ncol = Bm.ncol();
  int nrow = Bm.nrow();
  int i, j;
  NumericVector u(ncol,0.0);
  IntegerVector uset(ncol,0);
  double la1 = 1e+200;
  double la2 = 1e+200;
  double la, eps = 1e-9;
  int lead;
  double v, free = -1.0;
  // find nonzero vector u in Ker B (null space of B, i.e. Bu = 0)
  // with both positive and negative values
  // find reduced row echelon form of B
  // rref(Bm); -----------------------------------------------> remove here the rref (it is done in previous step)
  for(i=(nrow-1);i>=0;i--){
    // find lead (first nonzero entry on row) if exists
    // if no lead, i.e lead = ncol, do nothing
    // if lead, the variables after are either set or free
    // free variables are alternately set to 1 or -1 
    lead = 0;
    for(j=0;j<ncol;j++){if(Bm(i,j)==0.0){lead++;}else{break;}}
    // lead found
    if(lead<ncol){
      v = 0.0;
      for(j=lead+1;j<ncol;j++){
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
  for(i=0;i<ncol;i++){
    if( uset[i] == 0 ){
      free *= -1.0;
      u[i] = free;
    }else{break;}	
  }
  // find lambda1 and lambda2
  for(i=0;i<ncol;i++){
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
  for(i=0;i<ncol;i++){
    prob[i] = prob[i] + la * u[i];
    if(prob[i] < eps){ prob[i] = 0; }
    if(prob[i] > 1-eps){ prob[i] = 1; }
  }
  return prob;	
}