// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>

#include <cmath>
#include <algorithm>
#include <iostream>

#include "util.h"

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


// [[Rcpp::export]]
NumericVector dispen_map(NumericVector x, int ma, int nc, double mu, double sigma) {

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

// [[Rcpp::export]]
NumericVector dispen_npdf(NumericVector z, int nc, int m, int tau) {

  // locals
  int nm = pow(nc,m);
  vector<int> patterns(nm*m);
  NumericVector keys(nm);
  int N = z.size();

  // forge patterns array
  int p;
  int v=0;
  int k;
  int c1=0;
  int c2=0;
  for(int j=0; j<m; j++) {
    k = pow(nc,j);
    for(int i=0; i<nm; i++) {
      //p = j*nm + i;
      p = i*m + j;
      v = c1;
      patterns[p] = v + 1;
      if(++c2 >= k) {
        if(++c1 == nc) c1 = 0;
        c2 = 0;
      }
    }
  }

  // forge keys vector
  vector<int> keys_coef(m);
  for(int i=m; i>0; i--) keys_coef[m-i] = pow(10,i-1);
  for(int i=0; i<nm; i++) {
    p = i*m;
    c1 = 0;
    for(int j=0; j<m; j++) {
      c1 = c1 + patterns[p+j] * keys_coef[j];
    }
    keys[i] = c1;
  }

  // embd2
  arma::Col<int> embd2(N - (m-1)*tau,arma::fill::zeros);
  for(int j=1;j<=m;j++) {
    c1 = pow(10,m-j);
    for(int i=(j-1)*tau; i<(N -(m-j)*tau); i++) {
      embd2[i-(j-1)*tau] += z[i]*c1;
    }
  }

  // npdf
  NumericVector npdf(nm);
  arma::uvec hits;
  for(int id=0; id<nm; id++) {
    c1 = keys[id];
    hits = find(embd2 == c1);
    npdf[id] = (hits.n_elem > 0 ? (double) hits.n_elem /(N-(m-1)*tau) : 0.0);
  }

  // eop
  return npdf;
}

// [[Rcpp::export]]
NumericVector fdispen_npdf(NumericVector z, int nc, int m, int tau) {

  // init. patterns
  int N = z.size();
  int pm = (m == 1 ? 1 : m-1);
  int nm = 2*nc-1;
  int pn = pow(nm,pm);
  vector<int> patterns(pn*pm);

  // patterns
  int p;
  int v=0;
  int k;
  int c1=0;
  int c2=0;
  for(int j=0; j<pm; j++) {
    k = (j == 0 ? 1 : pow(nm,j));
    for(int i=0; i<pn; i++) {
      //p = j*nm + i;
      p = i*pm + j;
      v = c1;
      patterns[p] = v + 1;
      if(++c2 >= k) {
        if(++c1 == nm) c1 = 0;
        c2 = 0;
      }
    }
  }

  Rcerr << "$$ Hello There 1\n";

  // key
  NumericVector key(pn);
  for(int i=0; i<pn;i++) {
    key[i] = 0;
    for(int ii=0; ii<pm; ii++) {
      key[i] = key[i]*100 + patterns[i*pm + ii];
    }
  }

  Rcerr << "$$ Hello There 2\n";

  // ind matrix
  int hn = m;
  int hm = N-pm*tau;
  NumericMatrix ind(hn,hm);
  for(int i=0; i<hn; i++) {
    c1 = i*tau+1;
    for(int j=0; j<hm; j++)
      ind(i,j) = c1++;
  }

  Rcerr << "$$ Hello There 3\n";

  // embd2 / dembd2
  arma::mat embd2(hn,hm);
  for(int i=0; i<hn; i++) {
    for(int j=0; j<hm; j++) {
      embd2(i,j)=z[ind(i,j)-1];
    }
  }
  arma::mat dembd2 = diff(embd2).t() + nc;

  Rcerr << "$$ Hello There 4\n";

  // emb
  arma::mat foo(dembd2.n_rows, dembd2.n_cols);
  for(int i=pm; i>0; i--)
    foo.col(i-1) = dembd2.col(i-1) * pow(100,i-1);
  arma::vec emb = sum(foo,1);

  Rcerr << "$$ Hello There 5\n";

  // npdf
  NumericVector npdf(pn);
  arma::uvec hits;
  for(int id=0; id<pn; id++) {
    c1 = key[id];
    hits = find(emb == c1);
    npdf[id] = (hits.n_elem > 0 ? (double) hits.n_elem /(N-(m-1)*tau) : 0.0);
  }

  Rcerr << "$$ Hello There 6\n";

  // eop
  return npdf;
}
