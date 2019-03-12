#include <Rcpp.h>
using namespace Rcpp;

//' Sample Entropy
//'
//' Calcul the sample entropy
//'
//' @param x the signal as a numeric vector
//' @param m the embedding dimension
//' @param r the tolerance factor for "similarity"
//' @param sd the standard-deviation of the signal
//'
//' @return the sample entropy value
//'
//' @export
// [[Rcpp::export]]
double SampEn(NumericVector x, int m, double r, double sd)
{

  // init.
  int N = x.size();
  int N_max = N-(m+1)+1;
  int A = 0, B = 0;
  double sum = 0.0;
  double err;
  bool flag;

  // init.
  sd = (std::isnan(sd) ? Rcpp::sd(x) : sd);
  err = sd * r;

  //
  for (unsigned int i = 0; i < N_max; i++) {
    for (unsigned int j = i + 1; j < N_max; j++) {
      flag = true;
      for (unsigned int k = 0; k < m; k++) {
        if (std::abs(x[i+k] - x[j+k]) > err) {
          flag = false;
          break;
        }
      }
      if (flag) A++;
      if (flag && std::abs(x[i+m] - x[j+m]) <= err)
        B++;
    }
  }

  // eop
  if (A > 0 && B > 0)
    return std::log((double) A / (double)B);
  else
    return 0.0;
}
