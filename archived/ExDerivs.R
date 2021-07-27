## ExDerivs.R -- a file to explore different aspects of R in calculating Jacobians
## This is an experimental file for learning. 2021-5-25
rm(list=ls()) # clear the workspace (does not completely remove things, need to restart R sometimes)

# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
            38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)
eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))
weeddata1  <-  data.frame(y=ydat, tt=tdat)
weedenv <- list2env(weeddata1)
## weedenv$b1 <- start1$b1 ## Error in start1$b1 : $ operator is invalid for atomic vectors
weedenv$b1 <- start1[[1]]
weedenv$b2 <- start1[[2]]
weedenv$b3 <- start1[[3]]
rexpr<-call("-",eunsc[[3]], eunsc[[2]])
r0<-eval(rexpr, weedenv)
cat("Sumsquares at 1,1,1 is ",sum(r0^2),"\n")## Another way
expr <- eunsc
rho <- weedenv
rexpr<-call("-",eunsc[[3]], eunsc[[2]])
res0<-eval(rexpr, rho) # the base residuals
res0
## Try the numericDeriv option
theta<-names(start1)
nDnls<-numericDeriv(rexpr, theta, weedenv)
nDnls
nt <- length(theta) # number of parameters
mr <- length(res0) # number of residuals
prm<-as.vector(mget(theta,envir=rho))
JJ <- matrix(NA, nrow=mr, ncol=nt)
colnames(JJ)<-theta # May not be necessary
str(prm)
eps<-.Machine$double.eps
eps<-sqrt(eps) # see call in nls.R
dir<-rep(1,nt) # dir could be a vector in general??
for (j in 1:nt){
    origPar<-get(theta[j],rho)
    cat(theta[j]," = prm[[",j,"]]=",prm[[j]],"\n")
    xx <- abs(origPar)
    delta <- if (xx == 0.0) {eps} else { xx*eps }
    prmx<-origPar+delta*dir[j]
    assign(theta[j],prmx,rho)
    res1 <- eval(rexpr, rho) # new residuals (forward step)
    JJ[,j] <- dir[j]*(res1-res0)/delta
    assign(theta[j], origPar, envir=rho) # reset the parameter
}  
JJ
JJnD<-attr(nDnls,"gradient")
# Show that we have the same value
max(abs(JJ-JJnD))
rjfn<-nlsr::model2rjfun(eunsc, data=weeddata1, pvec=start1)
Ja <- attr(rjfn(start1),"gradient") # analytic derivatives
Ja
max(abs(Ja-JJ))
