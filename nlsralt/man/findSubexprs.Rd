\name{findSubexprs}
\alias{findSubexprs}
\title{
Find common subexpressions
}
\description{
This function finds common subexpressions in an expression vector
so that duplicate computation can be avoided.
}
\usage{
findSubexprs(expr, simplify = FALSE, tag = ".expr", verbose = FALSE, ...)
}
\arguments{
  \item{expr}{
An expression vector or language object.
}
  \item{simplify}{
Whether to call \code{\link{nlsSimplify}} on each subexpression before 
looking for common subexpressions.
}
  \item{tag}{
The prefix to use for locally created variables.
}
  \item{verbose}{
If \code{TRUE}, diagnostics are printed as simplifications are
recognized.
}
  \item{...}{
Additional parameters to pass to \code{\link{nlsSimplify}}.
Used only if \code{simplify} is \code{TRUE}.
}
}
\details{
This function identifies all repeated subexpressions in an expression vector, 
and stores them in locally created variables.  It is used by \code{\link{fnDeriv}}
to share common subexpressions between expression evaluations and gradient
evaluations, for example.

If \code{simplify} is \code{TRUE}, the assumptions behind the simplifications done
by \code{\link{nlsSimplify}} must be valid for the result to match the input.  With
the default simplifications, this means that all variables should take finite
real values.
}
\value{
A language object which evaluates to an expression vector which would
evaluate to the same result as the original vector with less 
duplicated code but more storage of intermediate results.

}

\author{
Duncan Murdoch
}
\seealso{
\code{\link{deriv}} in the \code{stats} package, \code{\link{nlsSimplify}}
}
\examples{
findSubexprs(expression(x^2, x-y, y^2-x^2))
}
\keyword{ math }
