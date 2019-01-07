#include <Rcpp.h>
using namespace Rcpp;

// for internal use only
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
