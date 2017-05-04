#include <Rcpp.h>
using namespace Rcpp;

inline double K(double x, double y, double h) {
  return 1/(2*M_PI*h*h)*exp(-0.5*(x*x+y*y)/(h*h));
}

// [[Rcpp::export(rng = false)]]
NumericVector kernel2d(NumericVector x, NumericVector y, NumericVector h) {
  NumericVector ans(Rcpp::no_init(x.length()));
  if (h.length() == 1) {
    double h2 = h[0];
    for(int i = 0; i < x.length(); i++) {
      ans[i] = K(x[i],y[i],h2);
    }
  } else {
    for(int i = 0; i < x.length(); i++) {
      ans[i] = K(x[i],y[i],h[i]);
    }
  }
  return ans;
}
