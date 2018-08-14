#ifndef SAMPLING_UTILITIES_H
#define SAMPLING_UTILITIES_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

double logit(double);
double invlogit(double);
double propose_positive(double, double);
double propose_prob(double, double);
double propose_norm(double, double);

double logit(double p){
  return((double) (log(p/(1-p))));
}

double invlogit(double x){
  return((double) (exp(x)/(1+exp(x))));
}

double propose_positive(double x, double mh_sd){
  double new_x = x + rnorm(1,0,mh_sd)[0];
  while(new_x <= 0){
    new_x = x + rnorm(1,0,mh_sd)[0];
  }
  return(new_x);
}

double propose_prob(double x, double mh_sd){
  double new_x = x + rnorm(1,0,mh_sd)[0];
  while(new_x < 0 || new_x > 1){
    if(new_x < 0)
      new_x = -new_x;
    if(new_x > 1)
      new_x = 2 - new_x;
  }
  return(new_x);
}

double propose_norm(double x, double mh_sd){
  double new_x = x + rnorm(1, 0, mh_sd)[0];
  return(new_x);
}

#endif
