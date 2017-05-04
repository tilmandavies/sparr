#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(rng = false)]]
NumericVector kernel2d(NumericVector x, NumericVector y, NumericVector h) {
  NumericVector ans(Rcpp::no_init(x.length()));
  if (h.length() == 1) {
    double h2 = h[0]*h[0];
    for(int i = 0; i < x.length(); i++) {
      ans[i] = 1/(2*M_PI*h2)*exp(-0.5*(x[i]*x[i]+y[i]*y[i])/h2);
    }
  } else {
    for(int i = 0; i < x.length(); i++) {
      ans[i] = 1/(2*M_PI*h[i]*h[i])*exp(-0.5*(x[i]*x[i]+y[i]*y[i])/(h[i]*h[i]));
    }
  }
  return ans;
}
