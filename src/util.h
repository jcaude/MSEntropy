#ifndef UTIL__H
#define UTIL__H

#include <Rcpp.h>
using namespace Rcpp;

extern double StandardDeviation(NumericVector x, double mu);
extern NumericVector normcdf(NumericVector x, double mu, double sigma);
extern NumericVector MapMinMax(NumericVector x, double ymin, double ymax);
extern NumericMatrix hankel(NumericVector x, NumericVector y);

#endif // !def UTIL_H
