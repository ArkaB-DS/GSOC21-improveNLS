\name{predict.nlsr}
\alias{predict.nlsr}
\title{Predictions for models specified as a formula of style y ~ (something)
}
\description{
    Function to allow predictions from nonlinear models estimated with
    \code{nlxb} from package \code{nlsr} if the model is specified by
    a formula with the structure y ~ (something).
}
\usage{
## S3 method for class 'nlsr'
 \method{predict}{nlsr}(object, newdata,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
     The output object of executing \code{nlxb} to estimate the model.
     An object of class 'nlsr'.}
  \item{newdata}{
     A named list containing the new data. This should be derived from
     a data frame of the same structure as the data used to estimate the
     model.}
  \item{\dots}{
     Any data needed for the function. We do not know of any! This is NOT
     currently used.
}
}
\details{

    To be added.
}
\value{
    Returns the predictions.
}
\author{
John C Nash <nashjc@uottawa.ca>
}
\examples{

cat("See examples in nlsr-package.Rd and in vignettes.\n")

}
\seealso{
    Function \code{nls()}.
}
\keyword{ nonlinear least squares }

