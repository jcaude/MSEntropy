#include <Rcpp.h>
using namespace Rcpp;

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

// for internal use only ... no longer needed (buggy fdisen)
NumericMatrix hankel(NumericVector x, NumericVector y) {

  // init.
  int xsize = x.size();
  int ysize = y.size();
  int dsize = xsize + ysize - 1;
  std::vector<int> data(dsize);
  NumericMatrix hmat(xsize,ysize);

  // fill data vector
  for(int i=0; i<dsize; i++) data[i] = (i >= xsize ? y[i-xsize+1] : x[i]);

  // fill matrix
  int delta = 0;
  for(int j=0;j<ysize; j++) {
    for(int i=0; i<xsize; i++) {
      hmat(i,j) = data[i + delta];
    }
    delta++;
  }

  // eop
  return hmat;
}
