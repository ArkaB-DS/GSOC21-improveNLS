## HobbsTiming.R -- computes jacobian of Hobbs at 1,1,1 using numericDeriv()
# rm(list=ls()) # clear the workspace (does not completely remove things, need to restart R sometimes)

traceval  <-  FALSE  # traceval set TRUE to debug or give full history

# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
            38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)
eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))
weeddata1  <-  data.frame(y=ydat, tt=tdat)

library(microbenchmark)
library(nlsralt)
library(nlsr)
tnlsrh1<-microbenchmark(nlsrh1<-nlxb(eunsc, data=weeddata1, start=start1))
tnlsrxh1<-microbenchmark(nlsrxh1<-nlxbx(eunsc, data=weeddata1, start=start1))
nlsrh1
nlsrxh1
tnlsrh1
tnlsrxh1

