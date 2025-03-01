#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector RRI_cpp(NumericVector x, size_t ni, size_t nw) {
  NumericVector out(ni);
  assert(nw == 25);

  auto rri_kernel = [](const double* x) -> double {
      return (abs(-0.5 * x[0] - 0.5 * x[12] - 0.207106781186547 * x[1] -
      0.207106781186547 * x[5] - 0.207106781186547 * x[7] -
      0.207106781186547 * x[11] + 1.82842712474619 * x[6]) +
      abs(-x[10] + 2 * x[11] - x[12]) + abs(-0.5 * x[20] -
      0.5 * x[12] - 0.207106781186547 * x[15] - 0.207106781186547 *
      x[21] - 0.207106781186547 * x[11] - 0.207106781186547 *
      x[17] + 1.82842712474619 * x[16]) + abs(-x[22] + 2 *
      x[17] - x[12]) + abs(-0.5 * x[24] - 0.5 * x[12] - 0.207106781186547 *
      x[19] - 0.207106781186547 * x[23] - 0.207106781186547 *
      x[13] - 0.207106781186547 * x[17] + 1.82842712474619 *
      x[18]) + abs(-x[14] + 2 * x[13] - x[12]) + abs(-0.5 *
      x[4] - 0.5 * x[12] - 0.207106781186547 * x[3] - 0.207106781186547 *
      x[9] - 0.207106781186547 * x[7] - 0.207106781186547 *
      x[13] + 1.82842712474619 * x[8]) + abs(-x[2] + 2 * x[7] -
      x[12]) + abs(-0.5 * x[6] - 0.5 * x[18] - 0.207106781186547 *
      x[7] - 0.207106781186547 * x[11] - 0.207106781186547 *
      x[13] - 0.207106781186547 * x[17] + 1.82842712474619 *
      x[12]) + abs(-x[11] + 2 * x[12] - x[13]) + abs(-0.5 *
      x[16] - 0.5 * x[8] - 0.207106781186547 * x[7] - 0.207106781186547 *
      x[11] - 0.207106781186547 * x[13] - 0.207106781186547 *
      x[17] + 1.82842712474619 * x[12]) + abs(-x[7] + 2 * x[12] -
      x[17]))/12.0;
  };

  const double* xptr = x.begin();

  for (auto& val : out) {
    val = rri_kernel(xptr);
    xptr += nw;
  }

  return out;
};
