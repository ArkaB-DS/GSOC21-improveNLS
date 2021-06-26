\name{resss}
\alias{resss}
\title{Compute sum of squares from residuals via the residual function.
}
\description{
   For a nonlinear model originally expressed as an expression of the form
         \code{lhs ~ formula_for_rhs}
   assume we have a \code{resfn} and \code{jacfn} that compute the residuals and the 
   Jacobian at a set of parameters. This routine computes the sum of squares of
   the residuals.
}
\usage{
   resss(prm, resfn, ...)
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
  \item{\dots}{
     Any data needed for computation of the residual vector from the expression
       rhsexpression - lhsvar. Note that this is the negative of the usual residual,
     but the sum of squares is the same.
}
}
\details{
    \code{resss} calls \code{resfn} to compute residuals and then uses \code{crossprod}
    to compute the sum of squares. 

    At 2012-4-26 there is no checking for errors. The evaluations of residuals and
    the cross product could be wrapped in \code{try()} if the evaluation could be
    inadmissible.
}
\value{
    The scalar numeric value of the sum of squares at the paramters.
}
\references{
%% ~put references to the literature/web site here ~
     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

     others!!

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
 
st  <-  c(b1=1, b2=1, b3=1)

firstss<-resss(st, shobbs.res)
# The sum of squares of the scaled Hobbs function at parameters
st
firstss
# now illustrate how to get solution via optimization

tf <- function(prm){
  val <- resss(prm, shobbs.res)
}
testop <- optim(st, tf, control=list(trace=1))
testop
}
\keyword{ nonlinear least squares }

