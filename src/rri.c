#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Visibility.h>


static size_t read_size(SEXP vec) {
  if (LENGTH(vec) != 1) Rf_error("scalar length expected");
  switch (TYPEOF(vec)) {
    case INTSXP: return INTEGER(vec)[0];
    case REALSXP: return (size_t) REAL(vec)[0];
    default: Rf_error("scalar length expected");
  };
}

// direct port of the RRI function
static double calculate_rri(const double* x) {
  return (
    fabs(-0.5 * x[0] - 0.5 * x[12] - 0.207106781186547 * x[1] -
    0.207106781186547 * x[5] - 0.207106781186547 * x[7] -
    0.207106781186547 * x[11] + 1.82842712474619 * x[6]) +
    fabs(-x[10] + 2 * x[11] - x[12]) + fabs(-0.5 * x[20] -
    0.5 * x[12] - 0.207106781186547 * x[15] - 0.207106781186547 *
    x[21] - 0.207106781186547 * x[11] - 0.207106781186547 *
    x[17] + 1.82842712474619 * x[16]) + fabs(-x[22] + 2 *
    x[17] - x[12]) + fabs(-0.5 * x[24] - 0.5 * x[12] - 0.207106781186547 *
    x[19] - 0.207106781186547 * x[23] - 0.207106781186547 *
    x[13] - 0.207106781186547 * x[17] + 1.82842712474619 *
    x[18]) + fabs(-x[14] + 2 * x[13] - x[12]) + fabs(-0.5 *
    x[4] - 0.5 * x[12] - 0.207106781186547 * x[3] - 0.207106781186547 *
    x[9] - 0.207106781186547 * x[7] - 0.207106781186547 *
    x[13] + 1.82842712474619 * x[8]) + fabs(-x[2] + 2 * x[7] -
    x[12]) + fabs(-0.5 * x[6] - 0.5 * x[18] - 0.207106781186547 *
    x[7] - 0.207106781186547 * x[11] - 0.207106781186547 *
    x[13] - 0.207106781186547 * x[17] + 1.82842712474619 *
    x[12]) + fabs(-x[11] + 2 * x[12] - x[13]) + fabs(-0.5 *
    x[16] - 0.5 * x[8] - 0.207106781186547 * x[7] - 0.207106781186547 *
    x[11] - 0.207106781186547 * x[13] - 0.207106781186547 *
    x[17] + 1.82842712474619 * x[12]) + fabs(-x[7] + 2 * x[12] -
    x[17])
  )/12.0;
}

// contract:
// - x must be a numeric() or integer()
// - ni and nw must be scalar lengths
//
// we do only minimal input validations here and rely on the terra package
SEXP batch_rri_impl(SEXP x, SEXP ni, SEXP nw) {
  size_t n = read_size(ni);
  if (read_size(nw) != 25) Rf_error("window size for RRI calculation must be 5");
  if (n*25 != XLENGTH(x)) Rf_error("input size mismatch");

  // output
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  double* outpr = REAL(out);

  // dispatch on the input type
  if (TYPEOF(x) == REALSXP) {
    const double* xptr = REAL(x);
    for (size_t i = 0; i < n; i++, xptr += 25) outpr[i] = calculate_rri(xptr);
  } else
  if (TYPEOF(x) == INTSXP) {
    const int* xptr = INTEGER(x);
    double buff[25];
    for (size_t i = 0; i < n; i++, xptr += 25) {
      // bulk-convert the window to double
      for (int j = 0; j < 25; j++) buff[j] = (xptr[j] == NA_INTEGER) ? NA_REAL : (double) xptr[j];
      outpr[i] = calculate_rri(buff);
    }
  } else {
    Rf_error("invalid input type (numeric expected)");
  }

  UNPROTECT(1);
  return out;
}

static const R_CallMethodDef CallEntries[] = {
    {"batch_rri_impl", (DL_FUNC) &batch_rri_impl, 3},
    {NULL, NULL, 0}
};

void R_init_SurfRough(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
