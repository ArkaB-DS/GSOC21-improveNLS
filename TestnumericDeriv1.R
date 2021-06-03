## TestnumericDeriv1.R 
# Data for Hobbs problem
rm(list=ls())
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
            38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

## library(nlspkg) # need these "installed" (or built in Rstudio)
## library(nlsalt)

# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)
eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))
weeddata1  <-  data.frame(y=ydat, tt=tdat)
weedenv <- list2env(weeddata1)
weedenv$b1 <- start1[[1]]
weedenv$b2 <- start1[[2]]
weedenv$b3 <- start1[[3]]
## ?? We seem to need the rexpr form.
rexpr<-call("-",eunsc[[3]], eunsc[[2]])
theta <- c("b1", "b2", "b3")
ndeunsc0<-nlspkg::numericDeriv(rexpr, theta, rho=weedenv)
print(ndeunsc0) # includes gradient
print(sum(ndeunsc0^2))

ndeunsc0c<-nlspkg::numericDeriv(rexpr, theta, rho=weedenv, central=TRUE)
print(ndeunsc0c) # includes gradient
print(sum(ndeunsc0c^2))

ndeunsc1<-nlsalt::numericDeriv(rexpr, theta, rho=weedenv)
print(ndeunsc1) # includes gradient
print(sum(ndeunsc1^2))

ndeunsc1c<-nlsalt::numericDeriv(rexpr, theta, rho=weedenv, central=TRUE)
print(ndeunsc1c) # includes gradient
print(sum(ndeunsc1c^2))

