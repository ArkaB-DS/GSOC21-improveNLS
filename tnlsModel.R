## tnlsModel.R -- run nlsModel and examine results
# rm(list=ls()) # clear the workspace (does not completely remove things, need to restart R sometimes)
# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
            38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)
eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))
weeddata1  <-  data.frame(y=ydat, tt=tdat)
weedenv <- list2env(weeddata1)
weedenv$b1 <- start1[[1]]
weedenv$b2 <- start1[[2]]
weedenv$b3 <- start1[[3]]
rexpr<-call("-",eunsc[[3]], eunsc[[2]])
theta <- c("b1", "b2", "b3")
library(nlsalt) # ?? needed because base R does not export nlsModel()
nmod1<-nlsModel(form=eunsc, data=weeddata1, start=start1, wts=NULL, upper=NULL, scaleOffset = 0, nDcentral = FALSE)
str(nmod1)
ls.str(nmod1)
print(nmod1)
