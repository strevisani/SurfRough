###update 3 March 2023###
#news
#1) Trik2 and RRI functions
#2) Madscan and Meanscan functions return the 3 layers raster with the names of the
# roughness indexes (names(result)=c("IsoRough","AnisoDir","AnisoR"))
#
###Date 1 September 2022 Venice###
#
##Begin of license
#Copyright (c) 2022 and 2023 Sebastiano Trevisani (strevisani@iuav.it)
#MIT type license
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#1)The above copyright notice and this permission notice shall be included in
#all copies or substantial portions of the Software.
#2)When using this software, in particular for publications, please cite the related papers:
#1) Trevisani, S. & Rocca, M. 2015. MAD: Robust image texture analysis for applications in high resolution geomorphometry.
#Computers and Geosciences, vol. 81, pp. 78-92. https://biblioproxy.cnr.it:2481/10.1016/j.cageo.2015.04.003
#2) Trevisani, S. Teza, G., Guth, P., 2023. A simplified geostatistical approach for characterizing key aspects of short-range roughness.
#CATENA,Volume 223, ISSN 0341-8162,https://doi.org/10.1016/j.catena.2023.106927
#3)Trevisani S., Teza G., Guth P.L., 2023. Hacking the topographic ruggedness index. Geomorphology
# https://doi.org/10.1016/j.geomorph.2023.108838
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#THE SOFTWARE.
##End of License
#
##Main functions implementing geostatistical-based
#surface/image texture indexes,using Terra package.
#With these functions you can
#compute classical geostatistical indexes
#(e.g., variogram and madogram) as well as the robust
#version MAD based on the median of absolute directional differences.
#
#You will find also some functions for computing roughness based
#on dispersion of normal vectors to surface (see paper for details).
#
#In this public version of the code, I did not included the functions for
#creating the kernels for computing the
#directional differences on the fly for any direction and lag.
#I provide a basic set of pre-computed kernels
#for calculating directional differences
#for basic lags (1,2,4,6,and 8 pixels) in four directions. I also provide the kernels for
#computing directional differences of order 2, permitting to calculate the new MADk2 based metrics,
#which do not require detrending.
#However, consider that in the development of ad-hoc kernels,
#not necessarily limited to bivariate indexes, there is a lot
#of potential for detecting interesting aspects of surface/image texture.
#Other potential relies in the various approaches for the decomposition of
#trend and residuals, and in multiscale smoothing approaches.
#Finally, when you need to study long-range distances you can resample
#the original DEM/image and the resampling should be function of the lag distance.
#
#These functions are developed as a starting point, for promoting
#creativity in surface/image texture analysis, and eventually implementing
#these in other software environments.
#

#load the required library Terra
#library(terra)
###Utilities###


#' Conversion from geographical directions to mathematical ones
#'
#' @param alpha The angle in degrees from one system of representation (e.g. geographical)
#'
#' @return The angle in degrees in the other system
#'
#' @noRd
matDeg<-function(alpha){
  eta=450-alpha
  if(eta>=360) eta=eta-360 else eta=eta
}

#conversion to radians
#' Conversion to radians
#'
#' @param degree Angle in degrees
#'
#' @return The angle in radians
#'
#' @noRd
rad<-function (degree) {
  (degree * pi)/180
}

#Load the kernels.
#I use square kernels, however in N-S and W-E directions
#you could use a vector (i.e., 1D kernel)
#The kernels are furnished with the source code.
#load("basicKernels.RData")
#load("kernels.RData")
###End Utilities###

######Surface texture functions######

#' Build a circular moving window
#'
#' @param radius The radius of the moving window
#'
#' @return A matrix with selected pixels
#' @export
#'
#' @examples
#' #A circular moving window with a radius of 3 pixels
#' w=KernelCircular(3)
#' w
KernelCircular=function(radius){
  #This function builds a simple circular kernel
  #with a given radius expressed in pixels
  size=radius*2+1
  center=c((size+1)/2,(size+1)/2)
  kernel=matrix(nrow=size,ncol=size,NA)
  for (i in 1:size) {
    for (j in 1:size){
      if (sqrt((j-center[1])**2+(i-center[1])**2)<=radius){kernel[i,j]=1}
    }
  }
  return(kernel)
}

#' Build a rectangular kernel of size X x Y
#'
#' @param lenx The size in pixels along x
#' @param leny The size in pixels along y
#'
#' @return The matrix (square/rectangular) with the selected pixels
#' @export
#'
#' @examples
#' #A rectangular moving window 5x5 pixels
#' w=KernelRectangular(5,5)
#' w
KernelRectangular=function(lenx,leny){
  #This function builds a simple rectangular kernel
  #with a given radius expressed in pixels
  kernel=matrix(nrow=leny,ncol=lenx,1)
  return(kernel)
}

#' Calculate the median of absolute values found in a search window for each raster in a list
#'
#' @param deltas A list of rasters with the values from which calculate the median of absolute values (e.g., directional differences of order K)
#' @param w The moving window used (e.g. w=KernelCircular(3))
#'
#' @return A list of rasters with the median of absolute values in the search window
#' @import terra
CalcMedians=function(deltas,w){
  #compute the medians of directional
  #absolute differences from a list
  #and returns a list of rasters.
  #Deltas-> list with rasters containing directional differences (of any order...)
  #w-> the search window (kernel, e.g., w=KernelCircular(3)).
  #For madogram copy the function and use the mean instead of the median in the
  # focal function, then divide by two (if you need it). If you need the variogram you should consider the
  #square of directional differences, calculate the mean, and divide by 2.

  medians=list()
  nlayers=nlyr(deltas)
  #Instead of a for cycle you could
  #use sapp...but generally we have
  #very few iterations
  for (i in 1:nlayers){
    medians[[i]]=focal(abs(deltas[[i]]),w,median,na.rm=FALSE,expand=F)
  }
  return(rast(medians))
}

#' Calculate the direction of maximum continuity considering 4 directions
#'
#' The input is represented by four rasters with the spatial variability index (e.g., MAD, variogram, etc.) computed
#' in four directions (N-S, NE-SW, E-W, SE-NW)
#' @param N Spatial variability along N-S direction
#' @param NE Spatial variability along NE-SW direction
#' @param E Spatial variability along E-W direction
#' @param SE Spatial variability along SE-NW direction
#'
#' @return A raster with the direction (in degrees, geographical) of maximum continuity
anisoDir=function(N,NE,E,SE){
  #AnisoDir calculates the direction of maximum continuity
  #using a circular statistics approach and using four directions
  #along N, NE, E and SE (in this precise order!).
  #Returns a raster.
  #So the input is MAD (or other spatial variability indexes)
  #calculated in the four directions
  180-(57.2957*0.5*atan2((NE-SE),(E-N)))
  #If you need more directions you should define a new function
  #taking as argument direction and modulus, quite easy.
}

#' Calculate the direction of maximum continuity considering 4 directions
#'
#' The input is represented by a list of rasters with the spatial variability index (e.g., MAD, variogram, etc.) computed
#' in four directions (N-S, NE-SW, E-W, SE-NW)
#' @param x A list of rasters with the spatial variability along 4 directions (see function anisoDir())
#' @importFrom terra lapp
#'
#' @return A raster with the direction (in degrees, geographical) of maximum continuity
anisoDirL=function(x){
  #AnisoDirL calculates the direction of maximum continuity
  #using a circular statistics approach and using a list of four
  #rasters of directional differences stored in the following order of
  #directions:N, NE, E and SE.
  ##Returns a raster.
  #See function anisoDir() for details.
  #anisoDir(x[[1]],x[[2]],x[[3]],x[[4]])
  lapp(x, anisoDir)
}

#' Calculate the index of anisotropy considering the spatial variability along 4 directions
#'
#' The input is represented by four rasters with the spatial variability index (e.g., MAD, variogram, etc.) computed
#' in four directions (N-S, NE-SW, E-W, SE-NW)
#' @param N Spatial vairability along N-S direction
#' @param NE Spatial vairability along NE-SW direction
#' @param E Spatial vairability along E-W direction
#' @param SE Spatial vairability along SE-NW direction
#'
#' @return A raster with the index of anisotropy (min=0 max=1)
anisoR=function(N,NE,E,SE){
  #Standardized resultant length. This is used as anisotropy index
  #Use four rasters with directional differences: N, NE,E, SE
  #Returns a raster.
  sqrt((NE-SE)**2+(E-N)**2)/(N+NE+E+SE)
}


#' Calculate the index of anisotropy considering the spatial variability along 4 directions
#'
#' The input is represented by a list of rasters with the spatial variability index (e.g., MAD, variogram, etc.) computed
#' in four directions (N-S, NE-SW, E-W, SE-NW)
#' @param x A list of rasters with the spatial variability along 4 directions (see function anisoR())
#' @importFrom terra lapp
#'
#' @return A raster with the index of anisotropy (min=0 max=1)
anisoRL=function(x){
  #Standardized resultant length. This is used as anisotropy index
  #Use a list of four rasters with directional differences: N, NE,E, SE
  #anisoR(x[[1]],x[[2]],x[[3]],x[[4]])
  #Returns a raster.
  lapp(x, anisoR)
}

#' Calculate MAD basic indexes
#'
#' Calculate MAD basic indexes considering a specif lag and difference of order K.
#' It computes 3 indexes of roughness/image texture: isotropic/omnidirectional; direction of maximum continuity; anisotropy index.
#' The anisotropy index is based on vector dispersion approach: 0 minimum anisotropy; 1 maximum anisotropy.
#' The direction of anisotropy is in degrees according to geographical convention.
#'
#'@references
#' 1) Trevisani, S. & Rocca, M. 2015. MAD: Robust image texture analysis for applications in high resolution geomorphometry.
#' Computers and Geosciences, vol. 81, pp. 78-92.
#'
#' 2) Trevisani, S. Teza, G., Guth, P., 2023. A simplified geostatistical approach for characterizing key aspects of short-range roughness.
#' CATENA,Volume 223, ISSN 0341-8162,https://doi.org/10.1016/j.catena.2023.106927
#'
#'
#' @param inRaster The DEM/residual-dem from which to compute the indexes
#' @param kernels The kernels to be used for computing the directional differences (e.g. order 1 or 2 for various lags)
#' @param w The moving window adopted for computing the geostatistical index (i.e., MAD)
#'
#' @return A list of 3 rasters: 1)isotropic roughness; 2) direction of anisotropy;3)index of anisotropy.
#' @import terra
#' @export
#'
#'
#' @examples
#' # MAD for lag 2 with differences of order 2 using a circular search window of radius 3.
#' # Using differences of order 1, you should
#' # apply these on a detrended surface/image.
#' dem=rast(paste(system.file("extdata", package = "SurfRough"), "/trento1.tif",sep=""))
#' w=KernelCircular(3)
#' rough2c=Madscan(dem,k2ck2, w)
#' #define understandable names of layers
#' names(rough2c)=c("IsoRough","AnisoDir","AnisoR")
#' #Plot isotropic roughness
#' plot(rough2c$IsoRough)
#' #Plot anisotropy index/strenght
#' plot(rough2c$AnisoR)
#'
Madscan<-function(inRaster,kernels,w){
  #Calculate MAD basic indexes based on 4 directions in this
  #order N,NE,SE,S.
  #Returns 3 rasters: 1)isotropic roughness; 2) direction of anisotropy;
  #3)index of anisotropy.
  #If you need more directions you need to generalize
  #the functions for anisotropy.
  #With very large files and few space on disk I use a modified version
  #that works a little differently. But I do not insert it in this public library.
  #inRaster->input raster, depending from the kernels
  #it may be a detrended version (i.e., high pass filtered) or directly the DTM/image.
  #kernels->a list of kernels (e.g.,myKernels=list(N2c,NE2c,E2c,SE2c)).
  #w->search window (e.g., w=KernelCircular(3)).
  deltas=list()
  #instead of a for loop you could use lapply
  for (i in 1:length(kernels)){
    deltas[[i]]=focal(inRaster, w=data.matrix(kernels[[i]]), na.rm=FALSE,expand=F)
  }
  #directional MADs
  deltas=rast(deltas)
  dirMad=CalcMedians(deltas,w)
  rm(deltas)
  #isotropic mad
  madIso=app(dirMad,fun=mean)
  #direction of maximum continuity
  anisoDirection=anisoDirL(dirMad)
  #anisotropy computed with circular statistics
  #standardized resultant length
  anisoR=anisoRL(dirMad)
  result=c(madIso,anisoDirection,anisoR)
  names(result)=c("IsoRough","AnisoDir","AnisoR")
  result
}

###Less robust geostatistical indexes###
#' Calculate the mean of absolute values raised to an exponent found in a search window
#'
#' With this you can compute variogram and madogram (but remember that for
#' classical geostatistical indexes you need to divide the derived isotropic index by 2!)
#
#' @param deltas The values from which calculate the median of absolute values (i.e., directional differences of order K)
#' @param w The moving window used (e.g. w=KernelCircular(3))
#' @param exponent The exponent: increasing the exponent increase the sensitivity to outliers. Set 2 for Variogram and 1 for Madogram.
#' @import terra
#'
#' @return A raster with the mean of absolute values in the search window
CalcMeans=function(deltas,w,exponent){
  #Compute the means of directional
  #absolute differences elevated at an exponent.
  #Returns a raster.
  #With this you can compute variogram and madogram (but remember that for
  #classical geostatistical indexes you need to divide by 2!)
  #Deltas-> list with rasters containing directional differences (of any order...)
  #w-> the search window (e.g., w=KernelCircular(3))
  #Exponent->the exponent to consider (e.g., 2 for variogram and 1 for madogram;
  #but other exponents are possible if you need).

  means=list()
  nlayers=nlyr(deltas)
  for (i in 1:nlayers){
    means[[i]]=focal((abs(deltas[[i]]))^exponent,w,mean,na.rm=FALSE,expand=F)
  }
  return(rast(means))
  #
}

#' Calculate less robust geostatistical indexes (mean of absolute differences raised to an exponent)
#'
#' With this you can compute variogram and madogram (but remember that for
#' classical geostatistical indexes you need to divide the derived isotropic index by 2!).
#' Moreover you can calibrate the exponent in order to filter or enhance hotspots and discontinuities
#' @param inRaster The DEM/residual-dem from which to compute the indexes
#' @param kernels The kernels to be used for computing the directional differences (e.g. order 1 or 2 for various lags)
#' @param w The moving window adopted for computing the geostatistical index (i.e., MAD)
#' @param exponent The exponent: increasing the exponent increase the sensitivity to outliers. Set 2 for Variogram and 1 for Madogram.
#' @return A SpatRaster with 3 layers: 1)isotropic roughness; 2) direction of anisotropy; 3)index of anisotropy.
#' @import terra
#' @export
#'
#'
#' @examples
#' #' Variogram-like for lag 2 with differences of order 2 using a circular search window of radius 3.
#' # Using differences of order 1, you should
#' # apply these on a detrended surface/image.
#' dem=rast(paste(system.file("extdata", package = "SurfRough"), "/trento1.tif",sep=""))
#' w=KernelCircular(3)
#' rough2c=Meanscan(dem,k2ck2, w,2)
#' # give understandable names to layers
#' names(rough2c)=c("IsoRough","AnisoDir","AnisoR")
#' #Plot isotropic roughness based on variogram estimator
#' #(divide by two if you need classical estimator)
#' plot(rough2c$IsoRough)
#'
Meanscan<-function(inRaster,kernels,w,exponent){
  #calculate basic indexes based on 4 directions in this
  #order N,NE,SE,S
  #if you need more directions you need to generalize
  #the functions for anisotropy.
  #Returns 3 rasters: 1)isotropic roughness; 2) direction of anisotropy;
  #3)index of anisotropy.
  #inRaster->input raster, depending from the kernels
  #it may be a detrended version (i.e., high pass filtered) or directly the DTM.
  #kernels->a list of kernels (e.g.,myKernels=list(N2c,NE2c,E2c,SE2c))
  #w->search window (e.g., w=KernelCircular(3))
  deltas=list()
  #instead of a for loop you could use lapply
  for (i in 1:length(kernels)){
    deltas[[i]]=focal(inRaster, w=data.matrix(kernels[[i]]), na.rm=FALSE,expand=F)
  }
  #directional MADs
  deltas=rast(deltas)
  dirMad=CalcMeans(deltas,w,exponent)
  rm(deltas)
  #isotropic variability
  madIso=app(dirMad,fun=mean)
  #direction of maximum continuity
  anisoDirection=anisoDirL(dirMad)
  #anisotropy computed with circular statistics
  #standardized resultant length
  anisoR=anisoRL(dirMad)
  result=c(madIso,anisoDirection,anisoR)
  names(result)=c("IsoRough","AnisoDir","AnisoR")
  result
  #
}
###End Less robust geostatistical indexes###


###Other roughness indexes###

#Here some roughness indexes related to Vector dispersion of normals to surface

#' Compute circular variance of aspect (i.e. of the gradient vector)
#'
#' @param inraster The DEM from which compute the index
#' @param window The moving window adopted for computing the index
#' @import terra
#' @return The raster with the computed index
#' @export
#'
#' @examples
#' # Gradient vector dispersion using a circular search window of radius 3.
#' dem=rast(paste(system.file("extdata", package = "SurfRough"), "/trento1.tif",sep=""))
#' w=KernelCircular(3)
#' roughGrad=circularDispersionGV(dem,w)
#' plot(roughGrad)
circularDispersionGV=function(inraster,window){
  #Circular dispersion of gradient vectors (analogous to circular variance of aspect)
  #using the mean resultant length approach.
  #working with angles in the geographical convention
  #window-> the search window/kernel (e.g., window=KernelCircular(3))
  slope=terrain(inraster,v="slope", unit= "radians") #Use radians!
  aspect=terrain(inraster,v="aspect", unit= "radians")
  
  x=cos(slope)*cos(aspect) 
  y=cos(slope)*sin(aspect)
  z=sin(slope)
  X=focal(x, w=window,fun=sum,expand=F,na.rm=F)
  Y=focal(y, w=window,fun=sum,expand=F,na.rm=F)
  Z=focal(z, w=window,fun=sum,expand=F,na.rm=F)
  R=sqrt(X^2+Y^2+Z^2)
  return(1-(R/sum(window,na.rm=T)))
}

###
#' Compute circular variance of normal vectors to surface
#'
#'Compute circular variance of normal vectors to surface, using the resultant vector length
#' @param inraster The DEM from which compute the index
#' @param window The moving window adopted for computing the index
#' @import terra
#' @return The raster with the computed index
#' @export
#'
#' @examples
#' #
#' #Normal vector dispersion using a circular search window of radius 3.
#' dem=rast(paste(system.file("extdata", package = "SurfRough"), "/trento1.tif",sep=""))
#' w=KernelCircular(3)
#' roughVDR=circularDispersionNV(dem,w)
#' plot(roughVDR)
circularDispersionNV=function(inraster,window){
  #Circular dispersion of normal vectors using the mean resultant length approach.
  #working with angles in the geographical convention
  #window-> the search window/kernel (e.g., window=KernelCircular(3))
  slope=terrain(inraster,v="slope", unit="radians") #Use radians!
  aspect=terrain(inraster,v="aspect", unit="radians")
  
  #respect to the formulas of Davis book
  # we consider that sin(90-slope)=cos(slope)
  # with 90-slope the angle of the normal vector respect to the horizontal plane
  x=sin(slope)*cos(aspect)
  y=sin(slope)*sin(aspect)
  z=cos(slope)
  X=focal(x, w=window,fun=sum,expand=F,na.rm=F)
  Y=focal(y, w=window,fun=sum,expand=F,na.rm=F)
  Z=focal(z, w=window,fun=sum,expand=F,na.rm=F)
  R=sqrt(X^2+Y^2+Z^2)
  return(1-(R/sum(window,na.rm=T)))
}

#vector dispersion using eigenvalues
#' For computing vector dispersion using eigen values ratios
#'
#' @param x Matrix cross products
#'
#' @return The dispersion/smoothness
#' @noRd
roory<-function(x){
  #function for using the eigenvalues approach
  x=matrix(x,nrow=3,ncol=3)
  if (anyNA(x)){
    smoothness=NA
  }
  else{
    ei=eigen(x)
    l1=ei$values[1]
    l2=ei$values[2]
    smoothness=log(l1/l2)
  }
  return(smoothness)
}
#'Compute circular variance of normal vectors to surface
#'
#'Compute circular variance of normal vectors to surface, using the eigen values (only for testing, very slow)
#' @param inraster The DEM from which compute the index
#' @param window The moving window adopted for computing the index
#' @import terra
#'
#' @return The raster with the computed index
circularEigenNV=function(inraster,window){
  #normal vector dispersion using the eigenvalues approach.
  #inraster->input raster (DTM, image, etc.)
  #window-> search window/kernel
  #for example w=KernelCircular(3).
  #it is not very efficient...and with topographic data, if using 2.5D representation,
  #gives analogous results to one based on circular dispersion via resultant length (after considering the log).
  #It could be coded more elegantly and efficiently, I coded it only for testing purposes.
  #With small variations of the code you can also calculate anisotropy, considering other eigenvalues
  #s2 and s3.

  slope=terrain(inraster,v="slope", unit="radians")
  aspect=terrain(inraster,v="aspect", unit="radians") #Use radians!

  #respect to the formulas of Davis book
  # we consider that sin(90-slope)=cos(slope)
  # with 90-slope the angle of the normal vector respect to the horizontal plane
  x=sin(slope)*cos(aspect)
  y=sin(slope)*sin(aspect)
  z=cos(slope)
  #Build the matrix of cross products
  x2=x^2
  y2=y^2
  z2=z^2
  xy=x*y
  xz=x*z
  yz=y*z
  X2=focal(x2, w=window,fun=sum,expand=F,na.rm=F)
  Y2=focal(y2, w=window,fun=sum,expand=F,na.rm=F)
  Z2=focal(z2, w=window,fun=sum,expand=F,na.rm=F)
  XY=focal(xy, w=window,fun=sum,expand=F,na.rm=F)
  XZ=focal(xz, w=window,fun=sum,expand=F,na.rm=F)
  YZ=focal(yz, w=window,fun=sum,expand=F,na.rm=F)
  cp<-c(X2,XY,XZ,XY,Y2,YZ,XZ,YZ,Z2)
  cp=values(cp)
  #colnames(cp)<-(c("X2","XY","XZ","XY","Y2","YZ","XZ","YZ","Z2"))
  eigen=apply(cp,1,roory)
  #print ("done")
  rooryz=X2
  return(init(rooryz,eigen))
}


#update 3 March 2023


#' Improved TRI (with differences of order 2), removing slope dependence.
#'
#' It is essentially a roughness radial index.
#' TRIk2 modifies TRI (topographic ruggedness index) using increments of order 2, symmetrical to central pixel,
#' so as to remove the effect of local slope.
#' This version does not correct for diagonal distance.
#' It uses a 5x5 kernel, consequently 12 directional differences of order k (2)
#' are used in the estimation.
#' One could also use a 3x3 kernel using only the 4 differences centered on the central pixel
#' but the metric would be very noisy.
#' The input is the DEM (no need to detrend).
#'
#' @references
#' 1) Riley, S. J., S. D. DeGloria, and R. Elliott. 1999.
#' A terrain ruggedness index that quantifies topographic heterogeneity.
#' Intermountain Journal of Science 5:23.
#' 2) Wilson, M.F.J., O'Connell, B., Brown, C., Guinan, J.C. & Grehan, A.J. 2007.
#' Multiscale terrain analysis of multibeam bathymetry data for habitat mapping on the continental slope".
#' Marine Geodesy, vol. 30, no. 1-2, pp. 3-35.
#' 3) Trevisani S., Teza G., Guth P.L., 2023 (Preprint). Hacking the topographic ruggedness index.
#' 10.5281/zenodo.7716785
#'
#' @param x The DEM from which to compute the index
#'
#'
#' @return isotropic roughness (in the same units of input)
#' @export
#'
#'
#' @examples
#' dem=rast(paste(system.file("extdata", package = "SurfRough"), "/trento1.tif",sep=""))
#' w <- matrix(1, nrow=5, ncol=5)
#' roughTrik5x5=focal(dem, w=w, fun=Trik2)
#' plot(roughTrik5x5)
#'
Trik2=function(x){
  (
    abs(-x[1]+2*x[7]-x[13])+
      abs(-x[11]+2*x[12]-x[13])+
      abs(-x[21]+2*x[17]-x[13])+
      abs(-x[23]+2*x[18]-x[13])+
      abs(-x[25]+2*x[19]-x[13])+
      abs(-x[15]+2*x[14]-x[13])+
      abs(-x[5]+2*x[9]-x[13])+
      abs(-x[3]+2*x[8]-x[13])+
      # You could define a function using only the following 4 directional differences
      abs(-x[7]+2*x[13]-x[19])+
      abs(-x[12]+2*x[13]-x[14])+
      abs(-x[17]+2*x[13]-x[9])+
      abs(-x[18]+2*x[13]-x[8])
  )/12
}

#' RRI: Radial Roughness index
#'
#' Modified TRI, based on increments of order 2  (removing slope dependence) and correcting for diagonal distance.
#' RRI modifies TRI (topographic ruggedness index) using increments of order 2, symmetrical to the central pixel,
#' so as to remove the effect of local slope.
#' This version corrects for the diagonal distance using bilinear interpolation.
#' It uses a 5x5 kernel, consequently 12 directional differences of order k (2)
#' are used in the estimation.
#' One could also use a 3x3 kernel using only the 4 differences centered on the central pixel
#' but the metric would be very noisy.
#' The input is the DEM (no need to detrend).
#'
#' @references
#'
#' 1) Riley, S. J., S. D. DeGloria, and R. Elliott. 1999.
#' A terrain ruggedness index that quantifies topographic heterogeneity.
#'  Intermountain Journal of Science 5:23.
#' 2) Wilson, M.F.J., O'Connell, B., Brown, C., Guinan, J.C. & Grehan, A.J. 2007.
#' Multiscale terrain analysis of multibeam bathymetry data for habitat mapping on the continental slope".
#' Marine Geodesy, vol. 30, no. 1-2, pp. 3-35.
#' 3) Trevisani S., Teza G., Guth P.L., 2023. Hacking the topographic ruggedness index. Geomorphology
#' https://doi.org/10.1016/j.geomorph.2023.108838
#'
#' @param x The DEM from which to compute the index
#'
#'
#' @return isotropic roughness (in the same units of input)
#' @export
#'
#'
#' @examples
#' dem=rast(paste(system.file("extdata", package = "SurfRough"), "/trento1.tif",sep=""))
#' w <- matrix(1, nrow=5, ncol=5)
#' roughTrick5x5=focal(dem, w=w, fun=RRI)
#' plot(roughTrick5x5)
#'
RRI=function(x){
  (
    #external differences
      abs(-0.5*x[1]-0.5*x[13]-0.207106781186547*x[2]-0.207106781186547*x[6]-0.207106781186547*x[8]-0.207106781186547*x[12]+1.82842712474619*x[7])+
      abs(-x[11]+2*x[12]-x[13])+
      abs(-0.5*x[21]-0.5*x[13]-0.207106781186547*x[16]-0.207106781186547*x[22]-0.207106781186547*x[12]-0.207106781186547*x[18]+1.82842712474619*x[17])+
      abs(-x[23]+2*x[18]-x[13])+
      abs(-0.5*x[25]-0.5*x[13]-0.207106781186547*x[20]-0.207106781186547*x[24]-0.207106781186547*x[14]-0.207106781186547*x[18]+1.82842712474619*x[19])+
      abs(-x[15]+2*x[14]-x[13])+
      abs(-0.5*x[5]-0.5*x[13]-0.207106781186547*x[4]-0.207106781186547*x[10]-0.207106781186547*x[8]-0.207106781186547*x[14]+1.82842712474619*x[9])+
      abs(-x[3]+2*x[8]-x[13])+
    #You could define a function using only the following 4 directional differences
      abs(-0.5*x[7]-0.5*x[19]-0.207106781186547*x[8]-0.207106781186547*x[12]-0.207106781186547*x[14]-0.207106781186547*x[18]+1.82842712474619*x[13])+
      abs(-x[12]+2*x[13]-x[14])+
      abs(-0.5*x[17]-0.5*x[9]-0.207106781186547*x[8]-0.207106781186547*x[12]-0.207106781186547*x[14]-0.207106781186547*x[18]+1.82842712474619*x[13])+
      abs(-x[8]+2*x[13]-x[18])
  )/12
}




###End other roughness indexes###


######End Surface texture functions######
