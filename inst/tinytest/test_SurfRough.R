#14 December 2022 Test routines
library(tinytest)
#for the R package
#using outputs from original code tested manually (here I use R outputs)
#load surfRough library
#library(SurfRoughV1)
#remove all objects
rm(list=ls())
#ls()
#get_call_wd()

#set the search window
w=KernelCircular(3)
w
#load the residual DTM
rast("residual.tif")->res

#res
#ck kernel of order 1
#step 1 calculate mad based indexes for lags 1 to 8
#
#don't care about error x$.self$finalize()
mad1c=Madscan(res,k1c,w)
mad2c=Madscan(res,k2c,w)
mad4c=Madscan(res,k4c,w)
mad6c=Madscan(res,k6c,w)
mad8c=Madscan(res,k8c,w)
rast("RIso1c.tif")->MAD1cIso
rast("RIso2c.tif")->MAD2cIso
rast("RIso4c.tif")->MAD4cIso
rast("RIso6c.tif")->MAD6cIso
rast("RIso8c.tif")->MAD8cIso
rast("RAnisoDir1c.tif")->MAD1cAnD
rast("RAnisoDir2c.tif")->MAD2cAnD
rast("RAnisoDir4c.tif")->MAD4cAnD
rast("RAnisoDir6c.tif")->MAD6cAnD
rast("RAnisoDir8c.tif")->MAD8cAnD
rast("RAniso1c.tif")->MAD1cAnR
rast("RAniso2c.tif")->MAD2cAnR
rast("RAniso4c.tif")->MAD4cAnR
rast("RAniso6c.tif")->MAD6cAnR
rast("RAniso8c.tif")->MAD8cAnR
#Step 2, calcualte differences
testMAD1cIso=MAD1cIso-mad1c[[1]]
testMAD2cIso=MAD2cIso-mad2c[[1]]
testMAD4cIso=MAD4cIso-mad4c[[1]]
testMAD6cIso=MAD6cIso-mad6c[[1]]
testMAD8cIso=MAD8cIso-mad8c[[1]]
testMAD1cAnD=MAD1cAnD-mad1c[[2]]
testMAD2cAnD=MAD2cAnD-mad2c[[2]]
testMAD4cAnD=MAD4cAnD-mad4c[[2]]
testMAD6cAnD=MAD6cAnD-mad6c[[2]]
testMAD8cAnD=MAD8cAnD-mad8c[[2]]
testMAD1cAnR=MAD1cAnR-mad1c[[3]]
testMAD2cAnR=MAD2cAnR-mad2c[[3]]
testMAD4cAnR=MAD4cAnR-mad4c[[3]]
testMAD6cAnR=MAD6cAnR-mad6c[[3]]
testMAD8cAnR=MAD8cAnR-mad8c[[3]]
#step 3 check if things work putting a threshold on max differences (minor differences are possible
#due to floating point representation and rounding.
#Then for direction of anisotropy can be higher
a=0
if(minmax(abs(testMAD1cIso))[2]<1e-5) a=a+1
if(minmax(abs(testMAD2cIso))[2]<1e-5) a=a+1
if(minmax(abs(testMAD4cIso))[2]<1e-5) a=a+1
if(minmax(abs(testMAD6cIso))[2]<1e-5) a=a+1
if(minmax(abs(testMAD8cIso))[2]<1e-5) a=a+1
if (a==5) 1 else 0
expect_equal(a, 5)
if(minmax(abs(testMAD1cAnR))[2]<1e-5) a=a+1
if(minmax(abs(testMAD2cAnR))[2]<1e-5) a=a+1
if(minmax(abs(testMAD4cAnR))[2]<1e-5) a=a+1
if(minmax(abs(testMAD6cAnR))[2]<1e-5) a=a+1
if(minmax(abs(testMAD8cAnR))[2]<1e-5) a=a+1
if (a==10) 1 else 0
expect_equal(a, 10)
if(minmax(abs(testMAD1cAnD))[2]<1e-4) a=a+1
if(minmax(abs(testMAD2cAnD))[2]<1e-4) a=a+1
if(minmax(abs(testMAD4cAnD))[2]<1e-4) a=a+1
if(minmax(abs(testMAD6cAnD))[2]<1e-4) a=a+1
if(minmax(abs(testMAD8cAnD))[2]<1e-4) a=a+1
if (a==15) 1 else 0
expect_equal(a, 15)

#Test kernels of order 2
#
#rm(list=ls())
rast("trento1.tif")->trento
#set the search window
w=KernelCircular(3)
w

#Step 1 calculate index
mad05ck2=Madscan(trento,k05ck2,w)
mad1ck2=Madscan(trento,k1ck2,w)
mad2ck2=Madscan(trento,k2ck2,w)
rast("rIso05ck2.tif")->MAD05ck2Iso
rast("rIso1ck2.tif")->MAD1ck2Iso
rast("rIso2ck2.tif")->MAD2ck2Iso
rast("rAnisoDir05ck2.tif")->MAD05ck2AnD
rast("rAnisoDir1ck2.tif")->MAD1ck2AnD
rast("rAnisoDir2ck2.tif")->MAD2cK2AnD
rast("rAniso05ck2.tif")->MAD05ck2AnR
rast("rAniso1ck2.tif")->MAD1ck2AnR
rast("rAniso2ck2.tif")->MAD2ck2AnR

#Step 2 calculate differences
testMAD05ck2Iso=MAD05ck2Iso-mad05ck2[[1]]
testMAD1ck2Iso=MAD1ck2Iso-mad1ck2[[1]]
testMAD2ck2Iso=MAD2ck2Iso-mad2ck2[[1]]
testMAD05ck2AnD=MAD05ck2AnD-mad05ck2[[2]]
testMAD1ck2AnD=MAD1ck2AnD-mad1ck2[[2]]
testMAD2ck2AnD=MAD2cK2AnD-mad2ck2[[2]]
testMAD05ck2AnR=MAD05ck2AnR-mad05ck2[[3]]
testMAD1ck2AnR=MAD1ck2AnR-mad1ck2[[3]]
testMAD2ck2AnR=MAD2ck2AnR-mad2ck2[[3]]
#step 3 check if things work putting a threshold on max differences (minor differences are possible
#due to floating point representation and rounding.
#Then for direction of anisotropy can be higher
a=0
if(minmax(abs(testMAD05ck2Iso))[2]<1e-5) a=a+1
if(minmax(abs(testMAD1ck2Iso))[2]<1e-5) a=a+1
if(minmax(abs(testMAD2ck2Iso))[2]<1e-5) a=a+1
if (a==3) 1 else 0
expect_equal(a, 3)
if(minmax(abs(testMAD05ck2AnR))[2]<1e-5) a=a+1
if(minmax(abs(testMAD1ck2AnR))[2]<1e-5) a=a+1
if(minmax(abs(testMAD2ck2AnR))[2]<1e-5) a=a+1
if (a==6) 1 else 0
expect_equal(a, 6)
if(minmax(abs(testMAD05ck2AnR))[2]<1e-4) a=a+1
if(minmax(abs(testMAD1ck2AnR))[2]<1e-4) a=a+1
if(minmax(abs(testMAD2ck2AnR))[2]<1e-4) a=a+1
if (a==9) 1 else 0
expect_equal(a, 9)
# end test kernels of order 2


#Test less robust geostatistical indexes
#####test
#rm(list=ls())
#load the residual DTM
rast("residual.tif")->res
#set the search window
w=KernelCircular(3)
w
exponent=2
gamma2c=Meanscan(res,k2c,w,2)
exponent=1
mado2c=Meanscan(res,k2c,w,1)
rast("rgammaIso2c.tif")->gamma2cIso
rast("rgammaAniso2c.tif")->gamma2cAnR
rast("rgammaAnisoDir2c.tif")->gamma2cAnD
rast("rmadoIso2c.tif")->mado2cIso
rast("rmadoAniso2c.tif")->mado2cAnR
rast("rmadoAnisoDir2c.tif")->mado2cAnD
test1=gamma2c[[1]]-gamma2cIso
test2=gamma2c[[3]]-gamma2cAnR
test3=gamma2c[[2]]-gamma2cAnD
test4=mado2c[[1]]-mado2cIso
test5=mado2c[[3]]-mado2cAnR
test6=mado2c[[2]]-mado2cAnD
a=0
if(minmax(abs(test1))[2]<1e-5) a=a+1
if(minmax(abs(test2))[2]<1e-5) a=a+1
if(minmax(abs(test3))[2]<1e-4) a=a+1
if (a==3) 1 else 0
expect_equal(a, 3)
if(minmax(abs(test4))[2]<1e-5) a=a+1
if(minmax(abs(test5))[2]<1e-5) a=a+1
if(minmax(abs(test6))[2]<1e-4) a=a+1
if (a==6) 1 else 0
expect_equal(a, 6)
#rm(list=ls())
#
###
#test trik and trick roughness indexes
a=0
w<- matrix(c(1:25), nrow=5, ncol=5)
if(Trik2(w)<1e-5) a=a+1
if(RRI(w)<1e-5) a=a+1
expect_equal(a, 2)
rm(list=ls())


