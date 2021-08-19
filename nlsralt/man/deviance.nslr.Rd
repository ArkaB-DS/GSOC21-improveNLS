\name{deviance.nlsr}
\alias{deviance.nlsr}
\title{Output sum of squares (deviance) for nlsr object.
}
\description{
    \code{deviance.nlsr} extracts the ssquares value for a model 
    estimated by \code{nlxb} or \code{nlfb} in the \code{nlsr} structured
    \code{object}. 
}
\usage{
    \method{deviance}{nlsr}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
     An object of class 'nlsr'
}
  \item{\dots}{
     Any data needed for the function. We do not know of any!
}
}
\details{
    \code{coef.nlsr} extracts and displays the coefficients for a model 
    estimated by \code{nlxb} or \code{nlfb}. 
}
\value{
    returns the coefficients from the nlsr object.
}
\author{
John C Nash <nashjc@uottawa.ca>
}
\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}
\keyword{ nonlinear least squares }
