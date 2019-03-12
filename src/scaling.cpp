#include <Rcpp.h>
using namespace Rcpp;

//' Signal Scaling (coarse-grained by average)
//'
//' Compute a coarse-grained copy of an input signal. The signal is average on
//' non-overlapping windows of size 'scale' (the scale factor).
//'
//' @param x the signal as a numeric vector
//' @param scale the scaling factor
//'
//' @return a numeric vector with the coarse-grained signal
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
