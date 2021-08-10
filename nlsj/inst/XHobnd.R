## XHobnd -- bounds test example for nlsj using Hobbs
rm(list=ls()) # clear workspace for each major section
## Use the Hobbs Weed problem
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weedframe <- data.frame(y=weed, tt=tt)
st <- c(b1=1, b2=1, b3=1) # a default starting vector (named!)
wmodu <- y ~ b1/(1+b2*exp(-b3*tt))
require(nlsr)
require(minpack.lm)
require(nlsj)
traceval<-TRUE
cat("Infeasible start test\n")
start1inf <- c(b1=4, b2=4, b3=4) # b3 OUT OF BOUNDS for next few tries
aji <- try(nlsj(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                   trace=traceval))
print(aji)
axi <- try(nlxb(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                   trace=traceval))
print(axi)
## Infeasible start! No warning message!
#anlM2i <- try(nlsLM(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
#                    trace=traceval))
#print(anlM2i)
# nls gives warnings
#anls2ia <- try(nls(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
#                   trace=traceval))
#anls2i <- try(nls(wmodu, start=start1inf, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
#                   algorithm="port", trace=traceval))
#print(anls2i)


# feasible start i.e. on or within bounds
start1<-st
anxb1 <- try(nlxb(wmodu, start=start1, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                   trace=traceval))
print(anxb1)
# Note that nlsj default (Gauss-Newton) fails to get "best" bounded answer
ajb1 <- try(nlsj(wmodu, start=start1, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                 trace=traceval))
print(ajb1)
# But marquardt variant works
ajb1m <- try(nlsj(wmodu, start=start1, data=weedframe, algorithm="marquardt", lower=c(0,0,0), upper=c(2, 6, 3), 
                 trace=traceval, control=nlsj.control(derivmeth="default")))
print(ajb1m)
# minpack.lm::nlsLM seems to get close in iterations, but then fails "singular gradient"
anlM1 <- try(nlsLM(wmodu, start=start1, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                    trace=traceval))
print(anlM1)
# nls gives warnings -- What is "singular convergence"? Does not return answer.
# Note "port" displays 1/2 sumsquares in trace info rather than ss
anls1 <- try(nls(wmodu, start=start1, data=weedframe, lower=c(0,0,0), upper=c(2, 6, 3), 
                  algorithm="port", trace=traceval))
print(anls1)
summary(anls1)
# Prepare some experiments with the residual/jacobian function
hbrj<-model2rjfun(wmodu, pvec=st, data=weedframe)
crossprod(hbrj(c(b1=2, b2=0, b3=2.57246)))
crossprod(hbrj(start1))
crossprod(hbrj(c(b1=2, b2=0, b3=3)))
11349.121*2
# ----------------------------------------------------------------- #
