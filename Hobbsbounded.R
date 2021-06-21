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
hobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
  # This variant uses looping
  if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
  y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
           38.558, 50.156, 62.948, 75.995, 91.972)
  tt  <-  1:12
  res  <-  x[1]/(1+x[2]*exp(-x[3]*tt)) - y
}

hobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
  jj  <-  matrix(0.0, 12, 3)
  tt  <-  1:12
  yy  <-  exp(-x[3]*tt)
  zz  <-  1.0/(1+x[2]*yy)
  jj[tt,1]   <-   zz
  jj[tt,2]   <-   -x[1]*zz*zz*yy
  jj[tt,3]   <-   x[1]*zz*zz*yy*x[2]*tt
  attr(jj, "gradient") <- jj
  jj
}
# Scaled Hobbs problem
shobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
  # This variant uses looping
  if(length(x) != 3) stop("shobbs.res -- parameter vector n!=3")
  y  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
           38.558, 50.156, 62.948, 75.995, 91.972)
  tt  <-  1:12
  res  <-  100.0*x[1]/(1+x[2]*10.*exp(-0.1*x[3]*tt)) - y
}

shobbs.jac  <-  function(x) { # scaled Hobbs weeds problem -- Jacobian
  jj  <-  matrix(0.0, 12, 3)
  tt  <-  1:12
  yy  <-  exp(-0.1*x[3]*tt)
  zz  <-  100.0/(1+10.*x[2]*yy)
  jj[tt,1]   <-   zz
  jj[tt,2]   <-   -0.1*x[1]*zz*zz*yy
  jj[tt,3]   <-   0.01*x[1]*zz*zz*yy*x[2]*tt
  attr(jj, "gradient") <- jj
  jj
}

##----- Hobbsbounded unique code starts here ------
require(nlsr)
require(minpack.lm)
traceval<-FALSE
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
