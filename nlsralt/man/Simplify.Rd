\name{nlsSimplify}
\alias{newSimplification}
\alias{nlsSimplify}
\alias{sysSimplifications}
\alias{isFALSE}
\alias{isZERO}
\alias{isONE}
\alias{isMINUSONE}
\alias{isCALL}
\title{
Functions to simplify expressions.
}
\description{
\code{nlsSimplify} simplifies expressions according to rules specified
by \code{newSimplification}.
}
\usage{
nlsSimplify(expr, simpEnv = sysSimplifications, verbose = FALSE)
newSimplification(expr, test, simplification, do_eval = FALSE, 
    simpEnv = sysSimplifications)
sysSimplifications
isFALSE(x)
isZERO(x)
isONE(x)
isMINUSONE(x)
isCALL(x, name)
}
\arguments{
  \item{expr}{
An expression to simplify; this should be a language object.
}
  \item{simpEnv}{
An environment holding the registered simplifications.
}  
  \item{verbose}{
If \code{TRUE}, diagnostics are printed as simplifications are recognized.
}
  \item{test}{
An expression giving a test to apply to decide whether this simplification
applies.
}
  \item{simplification}{
The new expression to use to replace the original one.
}
  \item{do_eval}{
Whether to evaluate the new expression (to give an even newer
expression) to use as the simplification.
}

  \item{x}{
An expression to test.
}
  \item{name}{
The name of a function as a character string.
}
}
\details{
The \code{nlsSimplify} function uses simple rules to simplify
expressions.  The simplification is aimed at the needs of this
package, so the built-in simplificatinos assume that variables and expressions 
have finite real values. For example, \code{0*expr} will simplify to \code{0}
regardless of the value of \code{expr}.  (The name is
\code{nlsSimplify} to avoid a clash with the \code{Simplify}
function in the \code{Deriv} package.)

\code{newSimplification} adds a new simplification pattern to the
registered collection.  The tests are applied to function calls with
the same function and number of parameters, in order as specified.
Users may specify their own environment (perhaps parented by
\code{sysSimplifications}) to hold rules if they wish to override the
standard rules.

The \code{isFALSE}, \code{isZERO}, \code{isONE}, and \code{isMINUSONE}
functions are simple functions to test whether expressions are simple constants,
similar to \code{\link{isTRUE}}.

The \code{isCALL} function tests whether an expression is a call to a particular
function.

To handle functions which act differently depending on which arguments
are present, \code{nlsSimplify} will simplify the expression \code{missing(.MissingVal)}
to \code{TRUE}.  This is used in the definition of the derivative
for \code{x - y}, where the unary minus is seen as a missing \code{y}
value.
}
\note{
The \code{isFALSE} function was added to base \R in version 3.5.0;
starting with that version, \code{nlsr::isFALSE} is simply a copy
of it.
}
\value{
\code{nlsSimplify} returns a simplification of the expression, assuming
that variables and functions take real values.

The \code{newSimplification} function is called for the side
effect of recording a new simplification rule.

If \code{expr} or \code{deriv} is missing, \code{newSimplification()} will
report on the currently saved simplifications in \code{simpEnv}.
}
\author{
Duncan Murdoch
}
\seealso{
\code{\link{nlsDeriv}}, which makes use of \code{nlsSimplify}.
}
\examples{
nlsSimplify(quote(a + 0))
nlsSimplify(quote(exp(1)), verbose = TRUE)

nlsSimplify(quote(sqrt(a + b)))  # standard rule
myrules <- new.env(parent = sysSimplifications)
newSimplification(sqrt(a), TRUE, a^0.5, simpEnv = myrules)
nlsSimplify(quote(sqrt(a + b)), simpEnv = myrules)
}
\keyword{ math }
