#include <cmath>
#include <algorithm>
#include <iostream>
#include <Rcpp.h>
using namespace std;
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#define MAP_LINEAR  1
#define MAP_NCDF    2
#define MAP_LOGSIG  3
#define MAP_TANSIG  4
#define MAP_SORT    5

#define EPSILON     1e-10

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' Signal Scaling (coarse-grained by average)
//'
//' @param x the signal as a numeric vector
//' @param scale the scaling factor
//'
//' @export
// [[Rcpp::export]]
NumericVector CoarseGraining(NumericVector x, int scale) {

  // scale == 1; no scaling
  if(scale == 1) return x;

  // init.
  int sx_size = x.size()/scale;
  NumericVector sx(sx_size);

  // avg. signal over scale windows
  for (int i=0; i < sx_size; i++) {
    sx[i] = 0;
    for(int k=0; k<scale; k++) {
      sx[i] += x[i*scale+k];
    }
    sx[i] /= scale;
  }

  // eop
  return sx;
}

// For internal use only
double StandardDeviation(NumericVector x, double mu) {

  // init.
  int xsize = x.size();
  double sum2 = 0.0;

  // loop
  for(int i=0; i<xsize; i++) {
    sum2 += (x[i] - mu)*(x[i] - mu);
  }

  // eop
  return(sqrt(sum2/(xsize - 1)));
}

// For internal use only
NumericVector normcdf(NumericVector x, double mu, double sigma) {

  // init.
  const double sqrt2=1.41421356237309504880;
  int xsize = x.size();
  NumericVector y(xsize);

  // loop over x..
  for(int i=0; i<xsize; i++)
    y[i] = erfc((mu-x[i])/(sigma*sqrt2))/2;

  // eop
  return y;
}

// For Internal use only
NumericVector MapMinMax(NumericVector x, double ymin, double ymax) {

  // init.
  int xsize = x.size();
  NumericVector y(xsize);
  double xmin = x[0];
  double xmax = x[0];

  // loop over x to compute xmin & xmax
  for(int i=0; i<xsize; i++) {
    if(x[i] < xmin) xmin = x[i];
    if(x[i] > xmax) xmax = x[i];
  }

  // loop agin to compute y ()
  if(xmin == xmax) {
    double yall = (xmin < ymin ? ymin : (xmax > ymax ? ymax : xmax) );
    for(int i=0; i<xsize; i++)
      y[i] = yall;
  }
  else {
    double ycoef = (ymax-ymin)/(xmax-xmin);
    for(int i=0; i<xsize; i++)
      y[i] = ycoef*(x[i]-xmin) + ymin;
  }

  // eop
  return y;
}

//' STEP-1
//'
//' @export
// [[Rcpp::export]]
NumericVector step1(NumericVector x, int ma, int nc, double mu, double sigma) {

  // init.
  int xsize = x.size();
  NumericVector z;
  NumericVector foo(xsize);

  // check/init mu_x
  double mu_x=(isnan(mu) ? accumulate(x.begin(), x.end(), 0.0)/ x.size() : mu);
  double sigma_x=(isnan(sigma) ? StandardDeviation(x,mu_x) : sigma);

  // switch according to the mapping mode
  switch(ma) {

  case MAP_LINEAR:
    z = MapMinMax(x,0,1);
    break;

  case MAP_NCDF:
    z = normcdf(x,mu_x,sigma_x);
    break;

  case MAP_LOGSIG:
    for(int i=0; i<xsize; i++)
      foo[i] = 1 / (1 + exp(-((x[i] - mu_x)/sigma_x)));
    z = MapMinMax(foo,0,1);
    break;

  case MAP_TANSIG:
    for(int i=0; i<xsize; i++)
      foo[i] = 2/(1+exp(-2*((x[i] - mu_x)/sigma_x)))-1;
    z = MapMinMax(foo,0,1);
    break;

  case MAP_SORT:

    // init.
    int N = (nc*floor(xsize/nc));
    NumericVector osx(N);
    NumericVector cx(N);
    NumericVector zN(N);

    // sort
    iota(osx.begin(),osx.end(),0);
    auto comparator = [&x](int a, int b){ return x[a] < x[b]; };
    sort(osx.begin(), osx.end(), comparator);
    for(int i=0; i<N; i++) cx[i] = ceil(i*nc/N)+1;
    for(int i=0; i<N; i++) zN[osx[i]] = cx[i];
    z = zN;
    break;

  }

  // final fix
  if(ma != MAP_SORT) {
    int zsize = z.size();
    for(int i=0; i<zsize; i++) {
      if(z[i] == 0) z[i] = EPSILON;
      if(z[i] == 1) z[i] = 1-EPSILON;
      z[i] = round(z[i]*nc + 0.5);
    }
  }

  // eop
  return z;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R

*/
