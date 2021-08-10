# XCrouchBd.R

rm(list=ls())
xdata = c(-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9)
ydata = c(0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001)
Cform <- ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata)
Cstart<-list(p1=1,p2=0.2)
Cdata<-data.frame(xdata, ydata)
library(nlsj)
library(nlsr)
library(minpack.lm)


Clower<-c(0,0)
Cupper<-c(1.5, 1.5)

fjb<-nlsj(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper)
fjb
fjbm<-nlsj(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper, algorithm="marquardt")
fjbm
fnb<-nls(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper, algorithm="port")
fnb
fxb<-nlxb(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper)
fxb
cbrj<-model2rjfun(Cform, data=Cdata, pvec=Cstart)
fmb<-nlsLM(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper)
fmb
crossprod(cbrj(coef(fxb)))
deviance(fnb)
coef(fnb)
deviance(fjbm)
coef(fjbm)
crossprod(cbrj(coef(fjbm)))
deviance(fmb)
coef(fmb)
crossprod(cbrj(coef(fxb)))
coef(fxb)
# ----------------------------------------------------------------
# Now change the bounds so p2 is smaller than allowed, p1 allowed bigger
Cstart<-list(p1=1,p2=0.2)

Clower<-c(0.4, 0.2)
Cupper<-c(2.5, 0.5)

# Note that "default" does not get good bounds, marquardt better
fjb<-nlsj(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper)
fjb
fjbm<-nlsj(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper, algorithm="marquardt")
fjbm
fnb<-nls(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper, algorithm="port")
fnb
fxb<-nlxb(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper)
fxb
cbrj<-model2rjfun(Cform, data=Cdata, pvec=Cstart)
fmb<-nlsLM(Cform, data=Cdata, start=Cstart, trace=TRUE, lower=Clower, upper=Cupper)
fmb
crossprod(cbrj(coef(fxb)))
deviance(fnb)
coef(fnb)
deviance(fjbm)
coef(fjbm)
crossprod(cbrj(coef(fjbm)))
deviance(fjb)
coef(fjb)
crossprod(cbrj(coef(fjb)))
deviance(fmb)
coef(fmb)
crossprod(cbrj(coef(fxb)))
coef(fxb)
# ----------------------------------------------------------------

