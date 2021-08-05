% File src/library/stats/man/nls.control.Rd
% Part of the R package, https://www.R-project.org
% Copyright 1995-2020 R Core Team
% Distributed under GPL 2 or later

\name{nlsj.control}
\alias{nlsj.control}
\title{Control the Iterations in nlsj}
\description{
  Allow the user to set some characteristics of the \code{\link{nls}}
  nonlinear least squares algorithm.
}
\usage{
nlsj.control(maxiter = 500, tol = 0.00001, minFactor = 1/1024,
	printEval = FALSE, warnOnly = FALSE, scaleOffset = 0,
    nDcentral = FALSE, watch = FALSE, phi = 1, lamda = 0, 
	offset = 100, laminc = 10, lamdec = 0.4, resmax = 10000, 
	rofftest = TRUE, smallsstest = TRUE,
	derivmeth="numericDeriv", altderivmeth="numericDeriv",
	trace=FALSE)
}
\arguments{
  \item{maxiter}{A positive integer specifying the maximum number of
    iterations allowed.}
  \item{tol}{A positive numeric value specifying the tolerance level for
    the relative offset convergence criterion.}
  \item{minFactor}{A positive numeric value specifying the minimum
    step-size factor allowed on any step in the iteration.  The
    increment is calculated with a Gauss-Newton algorithm and
    successively halved until the residual sum of squares has been
    decreased or until the step-size factor has been reduced below this
    limit.}
  \item{printEval}{a logical specifying whether the number of evaluations
    (steps in the gradient direction taken each iteration) is printed.}
  \item{warnOnly}{a logical specifying whether \code{\link{nls}()} should
    return instead of signalling an error in the case of termination
    before convergence.
    Termination before convergence happens upon completion of \code{maxiter}
    iterations, in the case of a singular gradient, and in the case that the
    step-size factor is reduced below \code{minFactor}.}
  \item{scaleOffset}{a constant to be added to the denominator of the relative
    offset convergence criterion calculation to avoid a zero divide in the case
    where the fit of a model to data is very close.  The default value of
    \code{0} keeps the legacy behaviour of \code{nls()}.  A value such as
    \code{1} seems to work for problems of reasonable scale with very small
    residuals.}
  \item{nDcentral}{only when \emph{numerical} derivatives are used:
    \code{\link{logical}} indicating if \emph{central} differences
    should be employed, i.e., \code{\link{numericDeriv}(*, central=TRUE)}
    be used.}
  \item{watch}{If TRUE, allows iterations to be interactive. Default = FALSE.
     NOT currently active.}
  \item{phi}{Number to use in a diagonal added to the augmented J' J Marquardt matrix
    in the Nash-Marquardt stabilization of the Gauss-Newton algorithm. Default = 1}
  \item{lamda}{Initial value for the lamda parameter used in Levenberg-Marquardt
    stabilization (note false spelling of lambda as a historical artifact that 
    often detects copying). When zero, the Gauss-Newton method is used. 
    ?? May need to coordinate with algorithm parameter of nlsj()}
  \item{offset}{Number to use in a tolerance-free comparison of two numbers. Default = 100.
      Items a and b are taken as equal if (a + offset) == (b + offset).}
  \item{laminc}{Factor to use to increase the Marquardt lamda parameter. Default = 10}
  \item{lamdec}{Factor to use to decrease the Marquardt lamda parameter. Default = 0.4} 
  \item{resmax}{Maximum number of residual evaluations allowed. Not used in nls(). Default 10000
       suggested is likely too big for general use.}
  \item{rofftest}{When TRUE (default), use the relative offset convergence criterion, 
    modified if appropriate with scaleOffset to deal with small residuals.}
  \item{smallsstest}{TRUE if we check for a very small sum of squares as an indicator
     for terminating the iterations. nls() only uses the relative offset test, and then
     only with scaleOffset=0.0. Nevertheless, suggest default of TRUE.}
  \item{derivmeth}{The method to computer derivatives for the Jacobian. To mimic nls(),
    this defaults to "numericDeriv".}
    \item{altderivmeth}{An alternate derivative method if "derivmeth" is infeasible,
      e.g., if derivmeth is "default" (analytic), then the alternative "numericDeriv"
      would call that routine. In future, there may be other methods.??}
  \item{trace}{When TRUE, allows progress information to be printed. ?? NEED
     to coordinate with argument "trace". Default=FALSE)}
}
\value{
  A \code{\link{list}} with components
  
  \item{maxiter}{}
  \item{tol}{}
  \item{minFactor}{}
  \item{printEval}{}
  \item{warnOnly}{}
  \item{scaleOffset}{}
  \item{nDcentreal}{}
  \item{watch}{}
  \item{phi}{}
  \item{lamda}{}
  \item{offset}{}
  \item{laminc}{}
  \item{lamdec}{}
  \item{resmax}{}
  \item{rofftest}{}
  \item{smallsstest}{}
  \item{derivmeth}{}
  \item{altderivmeth}{}
  \item{trace}{}
  
  with meanings as explained under \sQuote{Arguments}.
}
\references{
  Bates, D. M. and Watts, D. G. (1988),
  \emph{Nonlinear Regression Analysis and Its Applications}, Wiley.
}
\author{Douglas Bates and Saikat DebRoy; John C. Nash for extensions
  made in the 2020 Improvements to nls() Google Summer of Code 
  project with Arkajyoti Bhattacharjee.}
\seealso{
  \code{\link{nls}}
}
\examples{
nlsj.control(minFactor = 1/2048)
}
\keyword{nonlinear}
\keyword{regression}
\keyword{models}
