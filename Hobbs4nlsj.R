# Hobbsbounded-formula.R --- change name!!
#### bounds with formula specification of problem
rm(list=ls()) # clear workspace for each major section
## Use the Hobbs Weed problem
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weedframe <- data.frame(y=weed, tt=tt)
st <- c(b1=1, b2=1, b3=1) # a default starting vector (named!)
## Unscaled model
wmodu <- y ~ b1/(1+b2*exp(-b3*tt))
## Scaled model
wmods <- y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt))
## We can provide the residual and Jacobian as functions
# Unscaled Hobbs problem
##----- Hobbsbounded unique code starts here ------
require(nlsr)
require(minpack.lm)
require(nlsj)
require(nlsralt)
traceval<-TRUE
hunlsj01md <- try(nlsj(wmodu, start=st, data=weedframe, algorithm="marquardt", trace=traceval, control=nlsj.control(lamda=1e-4)))
print(hunlsj01md)

hunlsj01dn <- try(nlsj(wmodu, start=st, data=weedframe, algorithm="default", trace=traceval, control=nlsj.control(derivmeth="numericDeriv")))
print(hunlsj01dn)

hunlsj01dd <- try(nlsj(wmodu, start=st, data=weedframe, algorithm="default", trace=traceval, control=nlsj.control(derivmeth="default")))
print(hunlsj01dd)

hunlxbx01<-nlxbx(wmodu, start=st, data=weedframe, trace=traceval, control=list(watch=FALSE))
print(hunlxbx01)

hunlxb01<-nlxb(wmodu, start=st, data=weedframe, trace=traceval, control=list(watch=FALSE))
print(hunlxb01)

hump01<-nlsLM(wmodu, start=st, data=weedframe, trace=traceval, control=list(watch=FALSE))
summary(hump01)

hunls01 <- try(nls(wmodu, start=st, data=weedframe, trace=traceval))
summary(hunls01)

### ==============================================================
stgood<-c(b1=200, b2=50, b3=0.3)
hunlsj02md <- try(nlsj(wmodu, start=stgood, data=weedframe, algorithm="marquardt", trace=traceval, control=nlsj.control(lamda=1e-4)))
print(hunlsj02md)

hunlsj02dn <- try(nlsj(wmodu, start=stgood, data=weedframe, algorithm="default", trace=traceval, control=nlsj.control(derivmeth="numericDeriv")))
print(hunlsj02dn)

hunlsj02dd <- try(nlsj(wmodu, start=stgood, data=weedframe, algorithm="default", trace=traceval, control=nlsj.control(derivmeth="default")))
print(hunlsj02dd)

hunlxbx02<-nlxbx(wmodu, start=stgood, data=weedframe, trace=traceval, control=list(watch=FALSE))
print(hunlxbx02)

hunlxb02<-nlxb(wmodu, start=stgood, data=weedframe, trace=traceval, control=list(watch=FALSE))
print(hunlxb02)

hump02<-nlsLM(wmodu, start=stgood, data=weedframe, trace=traceval, control=list(watch=FALSE))
summary(hump02)

hunls02 <- try(nls(wmodu, start=stgood, data=weedframe, trace=traceval))
summary(hunls02)

### ==============================================================
cat("Infeasible start test\n")
st1inf <- c(b1=4, b2=4, b3=4) # OUT OF BOUNDS 
lb<-c(0, 0, 0)
ub<-c(2, 6, 3)

hunlsj03md <- try(nlsj(wmodu, start=st1inf, data=weedframe, algorithm="marquardt", lower=lb, upper=ub, trace=traceval, control=nlsj.control(lamda=1e-4)))
print(hunlsj03md)

hunlsj03dn <- try(nlsj(wmodu, start=st1inf, data=weedframe, algorithm="default", lower=lb, upper=ub, trace=traceval, control=nlsj.control(derivmeth="numericDeriv")))
print(hunlsj03dn)

hunlsj03dd <- try(nlsj(wmodu, start=st1inf, data=weedframe, algorithm="default", lower=lb, upper=ub, trace=traceval, control=nlsj.control(derivmeth="default")))
print(hunlsj03dd)

hunlxbx03<-nlxbx(wmodu, start=st1inf, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
print(hunlxbx03)

hunlxb03<-nlxb(wmodu, start=st1inf, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
print(hunlxb03)

hump03<-nlsLM(wmodu, start=st1inf, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
summary(hump03)

hunls03 <- try(nls(wmodu, start=st1inf, data=weedframe, lower=lb, upper=ub, algorithm="port", trace=traceval))
summary(hunls03)

### ==============================================================
# feasible start i.e. on or within bounds
## start = st = c(1, 1, 1)
lb<-c(0, 0, 0)
ub<-c(2, 6, 3)

hunlsj04md <- try(nlsj(wmodu, start=st1inf, data=weedframe, algorithm="marquardt", lower=lb, upper=ub, trace=traceval, control=nlsj.control(lamda=1e-4)))
print(hunlsj04md)

hunlsj04dn <- try(nlsj(wmodu, start=st1inf, data=weedframe, algorithm="default", lower=lb, upper=ub, trace=traceval, control=nlsj.control(derivmeth="numericDeriv")))
print(hunlsj04dn)

hunlsj04dd <- try(nlsj(wmodu, start=st1inf, data=weedframe, algorithm="default", lower=lb, upper=ub, trace=traceval, control=nlsj.control(derivmeth="default")))
print(hunlsj04dd)

hunlxbx04<-nlxbx(wmodu, start=st1inf, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
print(hunlxbx04)

hunlxb04<-nlxb(wmodu, start=st1inf, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
print(hunlxb04)

hump04<-nlsLM(wmodu, start=st1inf, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
summary(hump04)

hunls04 <- try(nls(wmodu, start=st1inf, data=weedframe, lower=lb, upper=ub,
algorithm="port", trace=traceval))
summary(hunls04)


### ==============================================================
# Single number bounds
lb <- 0
ub <- 3

hunlsj05md <- try(nlsj(wmodu, start=st1inf, data=weedframe, algorithm="marquardt", lower=lb, upper=ub, trace=traceval, control=nlsj.control(lamda=1e-4)))
print(hunlsj05md)

hunlsj05dn <- try(nlsj(wmodu, start=st1inf, data=weedframe, algorithm="default", lower=lb, upper=ub, trace=traceval, control=nlsj.control(derivmeth="numericDeriv")))
print(hunlsj05dn)

hunlsj05dd <- try(nlsj(wmodu, start=st1inf, data=weedframe, algorithm="default", lower=lb, upper=ub, trace=traceval, control=nlsj.control(derivmeth="default")))
print(hunlsj05dd)

hunlxbx05<-nlxbx(wmodu, start=st1inf, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
print(hunlxbx05)

hunlxb05<-nlxb(wmodu, start=st1inf, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
print(hunlxb05)

hump05<-nlsLM(wmodu, start=st1inf, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
summary(hump05)

hunls05 <- try(nls(wmodu, start=st1inf, data=weedframe, lower=lb, upper=ub,
algorithm="port", trace=traceval))
summary(hunls05)


### ==============================================================
# Single number bounds
lb <- c(0,0,0)
ub <- c(3,3,3)

hunlsj06md <- try(nlsj(wmods, start=st, data=weedframe, algorithm="marquardt", lower=lb, upper=ub, trace=traceval, control=nlsj.control(lamda=1e-4)))
print(hunlsj06md)

hunlsj06dn <- try(nlsj(wmods, start=st, data=weedframe, algorithm="default", lower=lb, upper=ub, trace=traceval, control=nlsj.control(derivmeth="numericDeriv")))
print(hunlsj06dn)

hunlsj06dd <- try(nlsj(wmods, start=st, data=weedframe, algorithm="default", lower=lb, upper=ub, trace=traceval, control=nlsj.control(derivmeth="default")))
print(hunlsj06dd)

hunlxbx06<-nlxbx(wmods, start=st, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
print(hunlxbx06)

hunlxb06<-nlxb(wmods, start=st, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
print(hunlxb06)

hump06<-nlsLM(wmods, start=st, data=weedframe, trace=traceval, lower=lb, upper=ub, control=list(watch=FALSE))
summary(hump06)

hunls06 <- try(nls(wmods, start=st, data=weedframe, lower=lb, upper=ub,
algorithm="port", trace=traceval))
summary(hunls06)

