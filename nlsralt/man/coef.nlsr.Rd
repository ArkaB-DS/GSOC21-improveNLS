\name{coef.nlsr}
\alias{coef.nlsr}
\title{Output model coefficients for nlsr object.
}
\description{
    \code{coef.nlsr} extracts and displays the coefficients for a model 
    estimated by \code{nlxb} or \code{nlfb} in the \code{nlsr} structured
    \code{object}. 
}
\usage{
    \method{coef}{nlsr}(object, ...)
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
