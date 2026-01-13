#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;
using namespace std;

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

// cpp version of the RRIk3 function (RRI with differences of order 3)
// [[Rcpp::export]]
NumericVector RRIK3_cpp(NumericVector x, size_t ni, size_t nw) {
  NumericVector out(ni);
  assert(nw == 25);

  auto rri_kernel = [](const double* x) -> double {
    return (
        abs(x[17]-3*x[12]+3*x[7]-x[2])
    +abs((0.207106781186548*x[11]+0.0857864376269049*x[12]+0.5*x[16]+0.207106781186547*x[17])
           -3*x[12]+3*(0.207106781186547*x[7]+0.5*x[8]+0.085786437626905*x[12]+0.207106781186547*x[13])
           -(0.242640687119285*x[3]+0.17157287525381*x[4]+0.34314575050762*x[8]+0.242640687119285*x[9]))
           +abs(x[11]-3*x[12]+3*x[13]-x[14])
           +abs((0.5*x[6]+0.207106781186548*x[7]+0.207106781186547*x[11]+0.085786437626905*x[12])
           -3*x[12]+3*(0.085786437626905*x[12]+0.207106781186547*x[13]+0.207106781186548*x[17]+0.5*x[18])
           -(0.34314575050762*x[18]+0.242640687119285*x[19]+0.242640687119285*x[23]+0.17157287525381*x[24]))
           +abs(x[7]-3*x[12]+3*x[17]-x[22])
           +abs((0.207106781186547*x[7]+0.5*x[8]+0.085786437626905*x[12]+0.207106781186547*x[13])
           -3*x[12]+3*(0.207106781186548*x[11]+0.0857864376269049*x[12]+0.5*x[16]+0.207106781186547*x[17])
           -(0.242640687119285*x[15]+0.34314575050762*x[16]+0.17157287525381*x[20]+0.242640687119285*x[21]))
           +abs(x[13]-3*x[12]+3*x[11]-x[10])
           +abs((0.085786437626905*x[12]+0.207106781186547*x[13]+0.207106781186548*x[17]+0.5*x[18])
           -3*x[12]+3*(0.5*x[6]+0.207106781186548*x[7]+0.207106781186547*x[11]+0.085786437626905*x[12])
           -(0.17157287525381*x[0]+0.242640687119285*x[1]+0.242640687119285*x[5]+0.34314575050762*x[6]))
    )/8;
  };

  const double* xptr = x.begin();

  for (auto& val : out) {
    val = rri_kernel(xptr);
    xptr += nw;
  }

  return out;
};

// cpp version of the RRIMin function (Minimum Radial Roughness index)
// [[Rcpp::export]]
NumericVector RRIMin_cpp(NumericVector x, size_t ni, size_t nw) {
  NumericVector out(ni);
  assert(nw == 25);

  auto rri_kernel = [](const double* x) -> double {
    array<double, 12> values = {
      abs(-0.5*x[0]-0.5*x[12]-0.207106781186547*x[1]-0.207106781186547*x[5]-0.207106781186547*x[7]-0.207106781186547*x[11]+1.82842712474619*x[6]),
      abs(-x[10]+2*x[11]-x[12]),
      abs(-0.5*x[20]-0.5*x[12]-0.207106781186547*x[15]-0.207106781186547*x[21]-0.207106781186547*x[11]-0.207106781186547*x[17]+1.82842712474619*x[16]),
      abs(-x[22]+2*x[17]-x[12]),
      abs(-0.5*x[24]-0.5*x[12]-0.207106781186547*x[19]-0.207106781186547*x[23]-0.207106781186547*x[13]-0.207106781186547*x[17]+1.82842712474619*x[18]),
      abs(-x[14]+2*x[13]-x[12]),
      abs(-0.5*x[4]-0.5*x[12]-0.207106781186547*x[3]-0.207106781186547*x[9]-0.207106781186547*x[7]-0.207106781186547*x[13]+1.82842712474619*x[8]),
      abs(-x[2]+2*x[7]-x[12]),
      abs(-0.5*x[6]-0.5*x[18]-0.207106781186547*x[7]-0.207106781186547*x[11]-0.207106781186547*x[13]-0.207106781186547*x[17]+1.82842712474619*x[12]),
      abs(-x[11]+2*x[12]-x[13]),
      abs(-0.5*x[16]-0.5*x[8]-0.207106781186547*x[7]-0.207106781186547*x[11]-0.207106781186547*x[13]-0.207106781186547*x[17]+1.82842712474619*x[12]),
      abs(-x[7]+2*x[12]-x[17])
    };

    // reduce the values
    return *min_element(begin(values), end(values));
  };

  const double* xptr = x.begin();

  for (auto& val : out) {
    val = rri_kernel(xptr);
    xptr += nw;
  }

  return out;
};


// cpp version of the RRIMax function (Maximum Radial Roughness index)
// [[Rcpp::export]]
NumericVector RRIMax_cpp(NumericVector x, size_t ni, size_t nw) {
  NumericVector out(ni);
  assert(nw == 25);

  auto rri_kernel = [](const double* x) -> double {
    array<double, 12> values = {
      abs(-0.5*x[0]-0.5*x[12]-0.207106781186547*x[1]-0.207106781186547*x[5]-0.207106781186547*x[7]-0.207106781186547*x[11]+1.82842712474619*x[6]),
      abs(-x[10]+2*x[11]-x[12]),
      abs(-0.5*x[20]-0.5*x[12]-0.207106781186547*x[15]-0.207106781186547*x[21]-0.207106781186547*x[11]-0.207106781186547*x[17]+1.82842712474619*x[16]),
      abs(-x[22]+2*x[17]-x[12]),
      abs(-0.5*x[24]-0.5*x[12]-0.207106781186547*x[19]-0.207106781186547*x[23]-0.207106781186547*x[13]-0.207106781186547*x[17]+1.82842712474619*x[18]),
      abs(-x[14]+2*x[13]-x[12]),
      abs(-0.5*x[4]-0.5*x[12]-0.207106781186547*x[3]-0.207106781186547*x[9]-0.207106781186547*x[7]-0.207106781186547*x[13]+1.82842712474619*x[8]),
      abs(-x[2]+2*x[7]-x[12]),
      abs(-0.5*x[6]-0.5*x[18]-0.207106781186547*x[7]-0.207106781186547*x[11]-0.207106781186547*x[13]-0.207106781186547*x[17]+1.82842712474619*x[12]),
      abs(-x[11]+2*x[12]-x[13]),
      abs(-0.5*x[16]-0.5*x[8]-0.207106781186547*x[7]-0.207106781186547*x[11]-0.207106781186547*x[13]-0.207106781186547*x[17]+1.82842712474619*x[12]),
      abs(-x[7]+2*x[12]-x[17])
    };

    // reduce the values
    return *max_element(begin(values), end(values));
  };

  const double* xptr = x.begin();

  for (auto& val : out) {
    val = rri_kernel(xptr);
    xptr += nw;
  }

  return out;
};


// cpp version of the TRIbi function (TRI corrected for diagonal distance)
// [[Rcpp::export]]
NumericVector TRIbi_cpp(NumericVector x, size_t ni, size_t nw) {
  NumericVector out(ni);
  assert(nw == 9);
  
  auto tri_kernel = [](const double* x) -> double {
    return (abs(x[4]-x[1])+
            abs(x[4]-(0.207106781186547*x[1]+0.5*x[2]+0.085786437626905*x[4]+0.207106781186547*x[5]))+
            abs(x[4]-x[5])+
            abs(x[4]-(0.085786437626905*x[4]+0.207106781186547*x[5]+0.207106781186548*x[7]+0.5*x[8]))+
            abs(x[4]-x[7])+
            abs(x[4]-(0.207106781186548*x[3]+0.0857864376269049*x[4]+0.5*x[6]+0.207106781186547*x[7]))+
            abs(x[4]-x[3])+
            abs(x[4]-(0.5*x[0]+0.207106781186548*x[1]+0.207106781186547*x[3]+0.085786437626905*x[4])))/8.0;
  };
  
  const double* xptr = x.begin();
  
  for (auto& val : out) {
    val = tri_kernel(xptr);
    xptr += nw;
  }
  
  return out;
};

// cpp version of the RRIcore function
// [[Rcpp::export]]
NumericVector RRIcore_cpp(NumericVector x, size_t ni, size_t nw) {
  NumericVector out(ni);
  assert(nw == 9);
  
  auto rricore_kernel = [](const double* x) -> double {
    return (abs(-0.5*x[0]-0.5*x[8]-0.207106781186547*x[1]-0.207106781186547*x[3]-0.207106781186547*x[5]-0.207106781186547*x[7]+1.82842712474619*x[4])+
            abs(-x[3]+2*x[4]-x[5])+
            abs(-0.5*x[6]-0.5*x[2]-0.207106781186547*x[1]-0.207106781186547*x[3]-0.207106781186547*x[5]-0.207106781186547*x[7]+1.82842712474619*x[4])+
            abs(-x[1]+2*x[4]-x[7]))/4.0;
  };
  
  const double* xptr = x.begin();
  
  for (auto& val : out) {
    val = rricore_kernel(xptr);
    xptr += nw;
  }
  
  return out;
};

// cpp version of the RRIk4 function
// [[Rcpp::export]]
NumericVector RRIk4_cpp(NumericVector x, size_t ni, size_t nw) {
  NumericVector out(ni);
  assert(nw == 25);
  
  auto rrik4_kernel = [](const double* x) -> double {
    return (abs(-x[22]+4*x[17]-6*x[12]+4*x[7]-x[2])+
            abs(-(0.242640687119285*x[15]+0.34314575050762*x[16]+0.17157287525381*x[20]+0.242640687119285*x[21])
                  +4*(0.207106781186548*x[11]+0.0857864376269049*x[12]+0.5*x[16]+0.207106781186547*x[17])
                  -6*x[12]+4*(0.207106781186547*x[7]+0.5*x[8]+0.085786437626905*x[12]+0.207106781186547*x[13])
                  -(0.242640687119285*x[3]+0.17157287525381*x[4]+0.34314575050762*x[8]+0.242640687119285*x[9])
            )+
            abs(-x[10]+4*x[11]-6*x[12]+4*x[13]-x[14])+
            abs(-(0.17157287525381*x[0]+0.242640687119285*x[1]+0.242640687119285*x[5]+0.34314575050762*x[6])
                  +4*(0.5*x[6]+0.207106781186548*x[7]+0.207106781186547*x[11]+0.085786437626905*x[12])
                  -6*x[12]+4*(0.085786437626905*x[12]+0.207106781186547*x[13]+0.207106781186548*x[17]+0.5*x[18])
                  -(0.34314575050762*x[18]+0.242640687119285*x[19]+0.242640687119285*x[23]+0.17157287525381*x[24]))
    )/4;
  };
  
  const double* xptr = x.begin();
  
  for (auto& val : out) {
    val = rrik4_kernel(xptr);
    xptr += nw;
  }
  
  return out;
};

// cpp version of the std function
// [[Rcpp::export]]
NumericVector std_cpp(NumericVector x, size_t ni, size_t nw) {
    NumericVector out(ni);
    
    auto sd_kernel = [nw](const double* window) -> double {
      // Create vector view of the window
      std::vector<double> values(window, window + nw);
      
      // Calculate mean
      double mean = std::accumulate(values.begin(), values.end(), 0.0) / nw;
      
      // Calculate variance using std::accumulate
      double variance = std::accumulate(values.begin(), values.end(), 0.0,
                                        [mean](double acc, double val) {
                                          return acc + (val - mean) * (val - mean);
                                        });
      
      // Sample standard deviation
      if (nw > 1) {
        //variance /= (nw - 1);
        variance /= nw;
      } else {
        variance = 0.0;
      }
      
      return std::sqrt(variance);
      
    };
    
    const double* xptr = x.begin();
    
    for (size_t i = 0; i < ni; ++i) {
      out[i] = sd_kernel(xptr);
      xptr += nw;
    }
    
    return out;
  };

// cpp version of the IQR function (method type=7)
// [[Rcpp::export]]
NumericVector iqr_cpp(NumericVector x, size_t ni, size_t nw) {
    NumericVector out(ni);
    
    auto iqr_kernel = [nw](const double* window) -> double {
      // Return NA if any value is NA (focalCpp default)
      for (size_t i = 0; i < nw; ++i) {
        //if (ISNA(window[i])) return NA_REAL;//gives issues in the border with NAs
        if (R_IsNA(window[i]) || R_IsNaN(window[i])) return NA_REAL;//this works properly
      }
      //if (nw < 2) return NA_REAL;
      
      std::vector<double> values(window, window + nw);
      std::sort(values.begin(), values.end());
      
      //Rs quantile method (type = 7)
      auto r_quantile = [&values, nw](double prob) -> double {
        double index = prob * (nw - 1);
        size_t lower = static_cast<size_t>(index);
        double fraction = index - lower;
        
        if (lower >= nw - 1) {
          return values.back();
        }
        
        return values[lower] + fraction * (values[lower + 1] - values[lower]);
      };
      
      double q1 = r_quantile(0.25);
      double q3 = r_quantile(0.75);
      
      return q3 - q1;
    };
    
    const double* xptr = x.begin();
    
    for (size_t i = 0; i < ni; ++i) {
      out[i] = iqr_kernel(xptr);
      xptr += nw;
    }
    
    return out;
  };
