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
jumq1m <- try(nlsj(wmodu, start=st, data=weedframe, algorithm="marquardt", trace=traceval, control=nlsj.control(lamda=1e-4)))
print(jumq1m)
tmp <- readline("next")
jr1m<-nlxbx(wmodu, start=st, data=weedframe, trace=traceval, control=list(watch=FALSE))
print(jr1m)
cat("Default usually fails from bad start\n")
tmp <- readline("next")
jumq1d <- try(nlsj(wmodu, start=st, data=weedframe, algorithm="default", trace=traceval))
print(jumq1d)
n01d <- try(nls(wmodu, start=st, data=weedframe, trace=traceval))
print(n01d)

stgood<-c(b1=200, b2=50, b3=0.3)
tmp <- readline("next")
jumq2m <- try(nlsj(wmodu, start=stgood, data=weedframe, algorithm="marquardt", trace=traceval, control=nlsj.control(lamda=1e-4)))
print(jumq2m)
tmp <- readline("next")
jumq2d <- try(nlsj(wmodu, start=stgood, data=weedframe, algorithm="default", trace=traceval))
print(jumq2d)
jumq2d <- try(nlsj(wmodu, start=stgood, data=weedframe, algorithm="default", trace=traceval))
print(jumq2d)
n02d <- try(nls(wmodu, start=stgood, data=weedframe, trace=traceval))
print(n02d)

cat("Infeasible start test\n")
start1inf <- c(b1=4, b2=4, b3=4) # b3 OUT OF BOUNDS for next few tries
anxb2i <- try(nlxb(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                   trace=traceval))
print(anxb2i)
## Infeasible start! No warning message!
anlM2i <- try(nlsLM(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                    trace=traceval))
print(anlM2i)
# nls gives warnings
anls2ia <- try(nls(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                   trace=traceval))
anls2i <- try(nls(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                   algorithm="port", trace=traceval))
print(anls2i)


# feasible start i.e. on or within bounds
start1<-st
anxb1 <- try(nlxb(wmodu, start=start1, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                   trace=traceval))
print(anxb1)
anlM1 <- try(nlsLM(wmodu, start=start1, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                    trace=traceval))
print(anlM1)
# nls gives warnings
anls1 <- try(nls(wmodu, start=start1, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                  algorithm="port", trace=traceval))
print(anls1)

# Hobbs scaled problem with bounds, formula specification
## cat("BUT ... nls() seems to do better from the TRACE information\n")
anlshob1b <- nls(wmods, start=start1, trace=traceval, data=weedframe, lower=c(0,0,0),
             upper=c(2,6,3), algorithm='port')
#  check the answer
print(anlshob1b)
cat("More precisely...crossprod(resid(anlsb))=",crossprod(resid(anlshob1b)),"\n")
## nlsLM seems NOT to work with bounds ?? try more examples?
anlsLM1b <- nlsLM(wmods, start=start1, data=weedframe, lower=c(0,0,0),
                 upper=c(2,6,3))
print(anlsLM1b)
newst<-coef(anlsLM1b)
print(newst)
anlsLM1bx <- nlsLM(wmods, start=newst, data=weedframe, lower=c(0,0,0),
                   upper=c(2,6,3))
print(anlsLM1bx)
anlx1bx <- nlxb(wmods, start=newst, data=weedframe, lower=c(0,0,0),
                   upper=c(2,6,3))
print(anlx1bx)

anlx1b <- nlxb(wmods, start=start1, data=weedframe, lower=c(0,0,0),
                  upper=c(2,6,3))
print(anlx1b)

cat("Single number bounds:\n")

anlx1b1 <- try(nlxb(wmods, start=start1, data=weedframe, lower=0,
                upper=3))
print(anlx1b1)

anls1b1 <- try(nls(wmods, start=start1, data=weedframe, lower=0,
                upper=3, algorithm="port"))
print(anls1b1)

anlsLM1b1 <- try(nlsLM(wmods, start=start1, data=weedframe, lower=0,
                   upper=3))
print(anlsLM1b1)
anlsLM1b1a <- try(nlsLM(wmods, start=start1, data=weedframe, lower=c(0,0,0),
                   upper=c(3,3,3)))
print(anlsLM1b1a)

anlm1b <- nls.lm(par=start1, fn=shobbs.res, jac=shobbs.jac, lower=c(0,0,0),
                  upper=c(2,6,3))
print(anlm1b)

anlf1b <- nlfb(start=start1, resfn=shobbs.res, jacfn=shobbs.jac, lower=c(0,0,0),
                upper=c(2,6,3))
print(anlf1b)

anlm1bx <- nls.lm(par=newst, fn=shobbs.res, jac=shobbs.jac, lower=c(0,0,0),
                   upper=c(2,6,3))
print(anlm1bx)

anlf1bx <- nlfb(start=newst, resfn=shobbs.res, jacfn=shobbs.jac, lower=c(0,0,0),
                upper=c(2,6,3))
print(anlf1bx)

shobbs.fn <- function(x) {
   resid <- shobbs.res(x)
   val<-as.numeric(crossprod(resid))
}

require(dfoptim)
anmkb1b <- nmkb(par=start1, fn=shobbs.fn, lower=c(0,0,0), upper=c(2,6,3))
anmkb1b
