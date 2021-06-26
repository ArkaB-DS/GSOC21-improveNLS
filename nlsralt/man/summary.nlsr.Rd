\name{summary.nlsr}
\alias{summary.nlsr}
\title{Summary output for nlsr object.
}
\description{
    Provide a summary output (but involving some serious computations!) of
    an object of class nlsr from \code{nlxb} or \code{nlfb} from package
    \code{nlsr}.
}
\usage{
    \method{summary}{nlsr}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
     An object of class 'nlsr'
}
  \item{\dots}{
     Currently ignored.
}
}
\details{
    \code{summary.nlsr} performs a summary method for an object of class 'nlsr' that
    has been created by a routine such as \code{nlfb} or \code{nlxb} for nonlinear
    least squares problems.

    Issue: When there are bounded parameters, \code{nls} returns a Standard Error for each of
    the parameters. However, this summary does NOT have a Jacobian value (it is set to 0)
    for columns where a parameter is masked or at (or very close to) a bound. See the 
    \code{R} code for the determination of whether we are at a bound. In this case, 
    users may wish to look in the \file{inst/dev-codes} directory of this package, 
    where there is a script \file{seboundsnlsrx.R} that computes the \code{nls()} 
    standard errors for comparison on a simple problem. 

    Issue: The \code{printsum()} of this object includes the singular values of the Jacobian. 
    These are displayed, one per coefficient row, with the coefficients. However, the 
    Jacobian singular values do NOT have a direct correspondence to the coefficients 
    on whose display row they appear. It simply happens that there are as many Jacobian 
    singular values as coefficients, and this is a convenient place to display them.
    The same issue applies to the gradient components.

}
\value{
    returns an invisible copy of the \code{nlsr} object.
}
\author{
John C Nash <nashjc@uottawa.ca>
}
\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}
\keyword{ nonlinear least squares }

