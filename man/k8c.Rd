% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{k8c}
\alias{k8c}
\title{basic kernels}
\format{
just matrices.
}
\source{
Sebastiano Trevisani
}
\usage{
k8c
}
\description{
Kernels for computing directional differences for specific directions and lag distances.
These have been constructed using bilinear interpolation for directions out of main axes.
The kernels are intended to be used with "Terra" focal functions (i.e., convolution).
}
\examples{
#to see kernels (each one is a list with 4 kernels) of order 1
#These should be used with a detrended "surface"
#lag 1 pixel
k1c
#lag 2 pixels
k2c
#lag 4 pixels
k4c
#lag 6 pixels
k6c
#lag 8 pixels
k8c
#kernels of order 2 (differences of differences)
#these can be applied directly without detrending
#lag 05 pixel
k05ck2
#lag 1 pixel
k1ck2
#lag 2 pixels
k2ck2
}
\references{
\enumerate{
\item Trevisani, S. & Rocca, M. 2015. MAD: Robust image texture analysis for applications in high resolution geomorphometry.
Computers and Geosciences, vol. 81, pp. 78-92.
\item Trevisani, S. Teza, G., Guth, P., 2023. A simplified geostatistical approach for characterizing key aspects of short-range roughness.
CATENA,Volume 223, ISSN 0341-8162,https://doi.org/10.1016/j.catena.2023.106927
}
}
\keyword{datasets}
