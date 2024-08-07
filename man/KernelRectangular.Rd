% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeostTextureLibrarySharedRoxy.R
\name{KernelRectangular}
\alias{KernelRectangular}
\title{Build a rectangular kernel of size X x Y}
\usage{
KernelRectangular(lenx, leny)
}
\arguments{
\item{lenx}{The size in pixels along x}

\item{leny}{The size in pixels along y}
}
\value{
The matrix (square/rectangular) with the selected pixels
}
\description{
Build a rectangular kernel of size X x Y
}
\examples{
#A rectangular moving window 5x5 pixels
w=KernelRectangular(5,5)
w
}
