% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeostTextureLibrarySharedRoxy.R
\name{anisoDir}
\alias{anisoDir}
\title{Calculate the direction of maximum continuity considering 4 directions}
\usage{
anisoDir(N, NE, E, SE)
}
\arguments{
\item{N}{Spatial variability along N-S direction}

\item{NE}{Spatial variability along NE-SW direction}

\item{E}{Spatial variability along E-W direction}

\item{SE}{Spatial variability along SE-NW direction}
}
\value{
A raster with the direction (in degrees, geographical) of maximum continuity
}
\description{
The input is represented by four rasters with the spatial variability index (e.g., MAD, variogram, etc.) computed
in four directions (N-S, NE-SW, E-W, SE-NW)
}
