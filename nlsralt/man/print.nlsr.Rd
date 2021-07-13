\name{print.nlsr}
\alias{print.nlsr}
\title{Print method for an object of class nlsr.
}
\description{
    Print summary output (but involving some serious computations!) of
    an object of class nlsr from \code{nlxb} or \code{nlfb} from package
    \code{nlsr}.
}
\usage{
    \method{print}{nlsr}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
     An object of class 'nlsr'
}
  \item{\dots}{
     Any data needed for the function. We do not know of any!
}
}
\details{
    \code{printsum.nlsr} performs a print method for an object of class 'nlsr' that
    has been created by a routine such as \code{nlfb} or \code{nlxb} for nonlinear
    least squares problems.
}
\value{
    Invisibly returns the input object.
}
\author{
John C Nash <nashjc@uottawa.ca>
}

\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}
\keyword{ nonlinear least squares }

