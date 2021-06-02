## ExDerivs.R -- a file to explore different aspects of R in calculating Jacobians
## This is an experimental file for learning. 2021-5-25
  
rm(list=ls()) # clear the workspace (does not completely remove things, need to restart R sometimes)

# library(nlsr)
# library(numDeriv)
# library(microbenchmark)
# 
# printsum <- function(xx){ print(summary(xx))}
# traceval  <-  TRUE  # traceval set TRUE to debug or give full history

# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
            38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)
eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))
# str(eunsc)
# ceunsc <- " y ~ b1/(1+b2*exp(-b3*tt))"
# str(ceunsc)
# print(as.formula(ceunsc)==eunsc)

# cat("LOCAL DATA IN DATA FRAMES\n")
weeddata1  <-  data.frame(y=ydat, tt=tdat)
weedenv <- list2env(weeddata1)
## weedenv$b1 <- start1$b1
## Error in start1$b1 : $ operator is invalid for atomic vectors
weedenv$b1 <- start1[[1]]
weedenv$b2 <- start1[[2]]
weedenv$b3 <- start1[[3]]
# ls.str(weedenv)
rexpr<-call("-",eunsc[[3]], eunsc[[2]])
r0<-eval(rexpr, weedenv)
cat("Sumsquares at 1,1,1 is ",sum(r0^2),"\n")## Another way
# ldata<-list2env(as.list(start1),envir=weedenv)
# ldata
# ls.str(ldata)
# eval(rexpr,envir=ldata)
# 
# ## How to get a model frame? Do we need it?
# 
# ## Error in model.frame.default(eunsc, ldata) : variable lengths differ (found for 'b1')
# mfeunsc<-model.frame(eunsc, ldata)
# 
# y<-ydat
# tt<-tdat
# b1<-1
# b2<-1
# b3<-1
# ## Error in model.frame.default(eunsc) : variable lengths differ (found for 'b1')
# mfeunsc<-model.frame(eunsc)
# get_all_vars(eunsc)
# eval(eunsc, ldata)
# str(mfeunsc)
# 
# ## 2021-5-25 try derivatives
# print(eunsc)
# ## nlsr::model2rjfun forms a function with gradient (jacobian) attribute
# funsc <- model2rjfun(eunsc, start1, data=weeddata1) # from nlsr: creates a function
# tmodel2rjfun <- microbenchmark(model2rjfun(eunsc, start1, data=weeddata1))
# print(tmodel2rjfun)
# print(funsc)
# print(funsc(start1))
# print(environment(funsc))
# print(ls.str(environment(funsc)))
# print(ls(environment(funsc)$data))
# eval(eunsc, environment(funsc))
# vfunsc<-funsc(start1)
# print(vfunsc)
# tfunsc<-microbenchmark(funsc(start1))
# print(tfunsc)
# 
# # following gives an error "theta should be of type character"
# nDfunsc<-try(numericDeriv(eunsc, start1, weeddata1))
# print(nDfunsc)
# print(start1)
# 
theta <- c("b1", "b2", "b3")
# Following gives error
## Error in numericDeriv(eunsc, theta, weedenv) : 
## 'language' object cannot be coerced to type 'double'

expr <- eunsc
rho <- weedenv
rexpr<-call("-",eunsc[[3]], eunsc[[2]])
res0<-eval(rexpr, rho) # the base residuals
res0
nt <- length(theta) # number of parameters
mr <- length(res0) # number of residuals
prm<-as.vector(mget(theta,envir=rho))
JJ <- matrix(NA, nrow=mr, ncol=nt)
colnames(JJ)<-theta # May not be necessary
str(prm)
eps<-.Machine$double.eps
dir<-rep(1,nt)
for (j in 1:nt){
    origPar <- prm[[j]]
    cat("prm[[",j,"]]=",prm[[j]],"\n")
    xx <- abs(origPar)
    delta <- if (xx == 0.0) {eps} else { xx*eps }
    prm[j]<- origPar * delta * dir[j]
    assign(theta[j], prm[[j]], envir=rho) # may be able to make more efficient later??
    res1 <- eval(rexpr, rho) # new residuals (forward step)
    for (i in 1:mr){
        JJ[i,j] <- dir[j]*(res1[i]-res0[i])/delta
    }
}  
JJ
