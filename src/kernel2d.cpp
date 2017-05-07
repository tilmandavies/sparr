#include <Rcpp.h>
#include "kernel_lut.h"

using namespace Rcpp;
using namespace kernel_lut;

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

// [[Rcpp::export(rng=false)]]
NumericVector kernel2d_adapt_none(NumericVector x, NumericVector y, NumericVector h, NumericVector ex, NumericVector ey) {
  // this evaluates the adaptive kernel density estimate of a point pattern (x,y) with bandwidth h(x,y) at (ex,ey)
  // h is assumed a vector, and we have no edge correction
  NumericVector ans(ex.length());
  // run over the evaluation points
  for (int e = 0; e < ex.length(); e++) {
    // run over the point pattern
    for(int i = 0; i < x.length(); i++) {
      ans[e] += K(x[i]-ex[e], y[i]-ey[e], h[i]);
    }
  }
  return(ans);
}

// [[Rcpp::export(rng=false)]]
NumericVector ec_uniform(NumericVector ex, NumericVector ey, NumericVector h, double dA) {
  // this computes the uniform edge correction for the adaptive kernel density estimate of a point pattern with bandwidth h(ex,ey) at (ex,ey)
  // h is assumed a vector of same length as ex,ey
  NumericVector ans(ex.length());
  for (int i = 0; i < ex.length(); i++) {
    // integrate over the evaluation points
    for(int e = 0; e < ex.length(); e++) {
      ans[i] += K(ex[e]-ex[i], ey[e]-ey[i], h[i]);
    }
    ans[i] *= dA;
  }
  return ans;
}

// [[Rcpp::export(rng=false)]]
NumericVector ec_diggle(NumericVector x, NumericVector y, NumericVector h, NumericVector ex, NumericVector ey, double dA) {
  // this computes the diggle edge correction for the adaptive kernel density estimate of a point pattern (x,y) with bandwidth h(x,y) at (ex,ey)
  // h is assumed a vector.
  NumericVector ans(x.length());
  for (int i = 0; i < x.length(); i++) {
    // integrate over the evaluation points
    for(int e = 0; e < ex.length(); e++) {
      ans[i] += K(ex[e]-x[i], ey[e]-y[i], h[i]);
    }
    ans[i] *= dA;
  }
  return ans;
}

// [[Rcpp::export(rng=false)]]
NumericVector kernel2d_adapt_uniform(NumericVector x, NumericVector y, NumericVector h, NumericVector edge, NumericVector ex, NumericVector ey) {
  // this evaluates the adaptive kernel density estimate of a point pattern (x,y) with bandwidth h(x,y) at (ex,ey)
  // h is assumed a vector, and we have edge correction ala Uniform (edge, of same length as ex,ey)
  NumericVector ans(ex.length());
  // run over the evaluation points
  for (int e = 0; e < ex.length(); e++) {
    // run over the point pattern
    for(int i = 0; i < x.length(); i++) {
      ans[e] += K(x[i]-ex[e], y[i]-ey[e], h[i]);
    }
    ans[e] /= edge[e];
  }
  return(ans);
}

// [[Rcpp::export(rng=false)]]
NumericVector kernel2d_adapt_diggle(NumericVector x, NumericVector y, NumericVector h, NumericVector edge, NumericVector ex, NumericVector ey) {
  // this evaluates the adaptive kernel density estimate of a point pattern (x,y) with bandwidth h(x,y) at (ex,ey)
  // h is assumed a vector, and we have edge correction ala Diggle (edge, of same length as x,y)
  NumericVector ans(ex.length());
  // run over the evaluation points
  for (int e = 0; e < ex.length(); e++) {
    // run over the point pattern
    for(int i = 0; i < x.length(); i++) {
      ans[e] += K(x[i]-ex[e], y[i]-ey[e], h[i])/edge[i];
    }
  }
  return(ans);
}
