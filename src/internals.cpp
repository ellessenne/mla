#include <Rcpp.h>

// #' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector valfpa(Rcpp::NumericVector const &vw,
                           Rcpp::NumericVector const &b,
                           Rcpp::NumericVector const &delta,
                           Rcpp::Function const &funcpa) {
  Rcpp::NumericVector bk = b + (exp(vw) * delta);
  Rcpp::NumericVector out = funcpa(bk);
  return(-out);
}
