\name{model2rjfunx}
\alias{model2rjfunx}
\alias{model2ssgrfunx}
\title{
Create functions to calculate the residual vector or the sum of squares, possibly
with derivatives.
}
\description{
Create functions to calculate the residual vector or the sum of squares, possibly
with derivatives.
}
\usage{
model2rjfunx(modelformula, pvec, data = NULL, jacobian = TRUE, testresult = TRUE, ...)
model2ssgrfunx(modelformula, pvec, data = NULL, gradient = TRUE, 
	      testresult = TRUE, ...)
}
\arguments{
  \item{modelformula}{
A formula describing a nonlinear regression model.
}
  \item{pvec}{
A vector of parameters.
}
  \item{data}{
A dataframe, list or environment holding data used in the calculation.
}
  \item{jacobian}{
Whether to compute the Jacobian matrix.
}
  \item{gradient}{
Whether to compute the Jacobian matrix.
}
  \item{testresult}{
Whether to test the function by evaluating it at \code{pvec}.
}
  \item{...}{
Additional data needed to compute the residuals and Jacobian.
}
}
\details{
If \code{pvec} does not have names, the parameters will have names
generated in the form \samp{p<n>}, e.g. \code{p1, p2}.  Names that appear in
\code{pvec} will be taken to be parameters of the model.  

The \code{data} argument may be a dataframe, list or environment, or \code{NULL}.
If it is not an environment, one will be constructed using the components
of \code{data} with parent environment set to be
the environment of \code{modelformula}.  

The search order for variables in the formula will be as follows:
Variables specified in the \code{...} argument will always take highest precedence.
If \code{data} is given and is an environment, the search continues only in that environment and its
parent, grandparent, etc.   Otherwise \code{data} is searched
next, followed by the environment of \code{modelformula} (and its parents, grandparents, etc.)

}
\value{

\code{model2rjfunx} returns a function with header \code{function(prm)}, which
evaluates the residuals (and if \code{jacobian} is \code{TRUE} the
Jacobian matrix) of the model at \code{prm}.  The residuals are defined to be
the right hand side of \code{modelformula} minus the left hand side.

\code{model2ssgrfunx} returns a function with header \code{function(prm)}, which
evaluates the sum of squared residuals (and if \code{gradient} is \code{TRUE} the
gradient vector) of the model at \code{prm}. 

}
\author{
John Nash and Duncan Murdoch
}

\seealso{
\code{\link{nls}}
}
\examples{

  y <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 38.558, 
      50.156, 62.948, 75.995, 91.972) 
      
  tt <- seq_along(y)  # for testing
  mydata <- data.frame(y = y, tt = tt)
  f <- y ~ b1/(1 + b2 * exp(-1 * b3 * tt)) + shift
  p <- c(b1 = 1, b2 = 1, b3 = 1)
# This fails
  first<-try(rjfn <- model2rjfunx(f, p, data = mydata, shift=100))
  if (class(first) != "try-error") {
    rjfn(p)
    myexp <- modelexpr(rjfn)
    cat("myexp:")
    print(myexp)
  } else { cat("rjfn try-error\n") }

#  ssgrfn <- model2ssgrfun(f, p, data = mydata, shift=100)
#  ssgrfn(p)
#  ssgrfn2 <- model2ssgrfun(f, p, data = mydata, shift=50)
#  ssgrfn2(p)
}
\keyword{ nonlinear }
