\name{resgr}
\alias{resgr}
\title{Compute gradient from residuals and Jacobian.
}
\description{
   For a nonlinear model originally expressed as an expression of the form
         lhs ~ formula_for_rhs
   assume we have a resfn and jacfn that compute the residuals and the 
   Jacobian at a set of parameters. This routine computes the gradient, 
   that is, t(Jacobian) . residuals. 
}
\usage{
   resgr(prm, resfn, jacfn, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prm}{
     A parameter vector. For our example, we could use
        start=c(b1=1, b2=2.345, b3=0.123)
     However, the names are NOT used, only positions in the vector.
}
  \item{resfn}{
     A function to compute the residuals of our model at a parameter vector.
}
  \item{jacfn}{
     A function to compute the Jacobian of the residuals at a paramter vector.
}
  \item{\dots}{
     Any data needed for computation of the residual vector from the expression
       rhsexpression - lhsvar. Note that this is the negative of the usual residual,
     but the sum of squares is the same.
}
}
\details{
    \code{resgr} calls resfn to compute residuals and jacfn to compute the
    Jacobian at the parameters \code{prm} using external data in the dot arguments.
    It then computes the gradient using t(Jacobian) . residuals.

    Note that it appears awkward to use this function in calls to optimization
    routines. The author would like to learn why. 
}
\value{
    The numeric vector with the gradient of the sum of squares at the paramters.
}
\references{

     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

}
\author{
John C Nash <nashjc@uottawa.ca>
}

\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}

\examples{
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
  attr(jj, "gradient") <- jj
  jj
}

st  <-  c(b1=1, b2=1, b3=1)
RG <- resgr(st, shobbs.res, shobbs.jac)
RG
}
\keyword{ nonlinear least squares }

