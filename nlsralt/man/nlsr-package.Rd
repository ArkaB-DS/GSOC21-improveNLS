\name{nlsr-package}
\alias{nlsr-package}
\alias{nlsr}
\docType{package}
\title{
Tools for solving nonlinear least squares problems. UNDER DEVELOPMENT.
}
\description{
The package provides some tools related to using the Nash variant
of Marquardt's algorithm for nonlinear least squares. 
}
\details{
\tabular{ll}{
Package: \tab nlsr\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-03-05\cr
License: \tab GPL-2\cr
}
This package includes methods for solving nonlinear least squares problems
specified by a modeling expression and given a starting vector of named 
paramters. Note: You must provide an expression of the form
    lhs ~ rhsexpression
so that the residual expression 
    rhsexpression - lhs
can be computed. The expression can be enclosed in quotes, and this seems to give
fewer difficulties with R. Data variables must already be defined, either within 
the parent environment or else in the dot-arguments. Other symbolic elements in
the modeling expression must be standard functions or else parameters that are 
named in the start vector.

The main functions in \code{nlsr} are:

\describe{
\item{nlfb}{Nash variant of the Marquardt procedure for nonlinear least squares,
	with bounds constraints, using a residual and optionally Jacobian
	described as \code{R} functions.}

\item{nlxb}{Nash variant of the Marquardt procedure for nonlinear least squares,
	with bounds constraints, using an expression to describe the residual via
        an \code{R} modeling expression. The Jacobian is computed via symbolic
	differentiation.}
            
\item{wrapnlsr}{Uses \code{nlxb} to solve nonlinear least squares then calls \code{nls()} to create an object of type nls.}

\item{model2rjfun}{returns a function with header \code{function(prm)}, which
        evaluates the residuals (and if jacobian is TRUE the
        Jacobian matrix) of the model at \code{prm}.  The residuals are defined 
        to be the right hand side of \code{modelformula} minus the left hand side.}

\item{model2ssgrfun}{returns a function with header \code{function(prm)}, which
        evaluates the sum of squared residuals (and if gradient is \code{TRUE} the
        gradient vector) of the model at \code{prm}.} 

\item{modelexpr}{returns the expression used to calculate the vector of 
        residuals (and possibly the Jacobian) used in the previous functions.}
}


}
\author{ John C Nash and Duncan Murdoch

Maintainer:  <nashjc@uottawa.ca>
}

\references{
     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

     Nash, J. C. (2014) _Nonlinear Parameter Optimization Using R Tools._
     Wiley 

}
\keyword{ nls }
\keyword{ nonlinear least squares }
\seealso{\code{nls} }
\examples{

rm(list=ls())
# require(nlsr)

traceval  <-  TRUE  # traceval set TRUE to debug or give full history

# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)
startf1  <-  c(b1=1, b2=1, b3=.1)

eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))

cat("LOCAL DATA IN DATA FRAMES\n")
weeddata1  <-  data.frame(y=ydat, tt=tdat)
weeddata2  <-  data.frame(y=1.5*ydat, tt=tdat)

anlxb1  <-  try(nlxb(eunsc, start=start1, trace=traceval, data=weeddata1))
print(anlxb1)

# illustrate predict
newdta <- colMeans(weeddata1)
newdta["tt"]<-25 # This only works for 1D example -- CAUTION
predn <- predict(anlxb1, as.list(newdta))
print(predn)

anlxb2  <-  try(nlxb(eunsc, start=start1, trace=traceval, data=weeddata2))
print(anlxb2)


escal  <-   y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt))
suneasy  <-  c(b1=200, b2=50, b3=0.3)
ssceasy  <-  c(b1=2, b2=5, b3=3)
st1scal  <-  c(b1=100, b2=10, b3=0.1)



shobbs.res  <-  function(x){ # scaled Hobbs weeds problem -- residual
# This variant uses looping
    if(length(x) != 3) stop("hobbs.res -- parameter vector n!=3")
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
    attr(jj,"gradient")<-jj # needed for nlfb to function seamlessly with nlxb
    return(jj)
}

cat("try nlfb\n")
st  <-  c(b1=1, b2=1, b3=1)
low  <-  -Inf
up <- Inf


\dontrun{

ans1 <- nlfb(st, shobbs.res, shobbs.jac, trace=traceval)
ans1
cat("No jacobian function -- use internal approximation\n")
ans1n <- nlfb(st, shobbs.res, trace=TRUE, control=list(watch=TRUE)) # NO jacfn
ans1n

# tmp <- readline("Try with bounds at 2")
time2 <- system.time(ans2 <- nlfb(st, shobbs.res, shobbs.jac, upper=c(2,2,2), 
                                  trace=traceval))
ans2
time2


} # end dontrun



cat("BOUNDS")
st2s <- c(b1=1, b2=1, b3=1)

\dontrun{

an1qb1 <- try(nlxb(escal, start=st2s, trace=traceval, data=weeddata1, 
  lower=c(0,0,0), upper=c(2, 6, 3), control=list(watch=FALSE)))
print(an1qb1)
tmp <- readline("next")
ans2 <- nlfb(st2s,shobbs.res, shobbs.jac, lower=c(0,0,0), upper=c(2, 6, 3), 
   trace=traceval, control=list(watch=FALSE))
print(ans2)

cat("BUT ... nls() seems to do better from the TRACE information\n")
anlsb <- nls(escal, start=st2s, trace=traceval, data=weeddata1, lower=c(0,0,0),
     upper=c(2,6,3), algorithm='port')
cat("However, let us check the answer\n")
print(anlsb)
cat("BUT...crossprod(resid(anlsb))=",crossprod(resid(anlsb)),"\n")

} # end dontrun


tmp <- readline("next")
cat("Try wrapnlsr\n")
traceval <- TRUE
# Data for Hobbs problem
ydat <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat <- seq_along(ydat) # for testing
start1 <- c(b1=1, b2=1, b3=1)
escal <-  y ~ 100*b1/(1+10*b2*exp(-0.1*b3*tt))
up1 <- c(2,6,3)
up2 <- c(1, 5, 9)

weeddata1 <- data.frame(y=ydat, tt=tdat)


\dontrun{

an1w <- try(wrapnlsr(escal, start=start1, trace=traceval, data=weeddata1))
print(an1w)



cat("BOUNDED wrapnlsr\n")

an1wb <- try(wrapnlsr(escal, start=start1, trace=traceval, data=weeddata1, upper=up1))
print(an1wb)


cat("BOUNDED wrapnlsr\n")

an2wb <- try(wrapnlsr(escal, start=start1, trace=traceval, data=weeddata1, upper=up2))
print(an2wb)

cat("TRY MASKS ONLY\n")

an1xm3 <- try(nlxb(escal, start1, trace=traceval, data=weeddata1, 
                   masked=c("b3")))
printsum(an1xm3)
#an1fm3 <- try(nlfb(start1, shobbs.res, shobbs.jac, trace=traceval, 
#                   data=weeddata1, maskidx=c(3)))
an1fm3 <- try(nlfb(start1, shobbs.res, shobbs.jac, trace=traceval, 
                   data=weeddata1, maskidx=c(3)))printsum(an1fm3)

an1xm1 <- try(nlxb
(escal, start1, trace=traceval, data=weeddata1, 
                   masked=c("b1")))
printsum(an1xm1)
#an1fm1 <- try(nlfb(start1, shobbs.res, shobbs.jac, trace=traceval, 
 an1fm1 <- try(nlfb(start1, shobbs.res, shobbs.jac, trace=traceval, 
                   data=weeddata1, maskidx=c(1)))
printsum(an1fm1)

} # end dontrun


# Need to check when all parameters masked.

}
