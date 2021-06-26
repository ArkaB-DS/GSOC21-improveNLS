\name{wrapnlsr}
\alias{wrapnlsr}
\title{Provides class nls solution to a nonlinear least squares solution 
     using the Nash Marquardt tools.
}
\description{
   Given a nonlinear model expressed as an expression of the form
         \code{lhs ~ formula_for_rhs}
   and a start vector where parameters used in the model formula are named,
   attempts to find the minimum of the residual sum of squares using the
   Nash variant (Nash, 1979) of the Marquardt algorithm, where the linear 
   sub-problem is solved by a qr method. The resulting solution is fed
   into the \code{nls()} function in an attempt to get the nls class solution.
}
\usage{
   wrapnlsr(formula, start, trace=FALSE, data, lower=-Inf, upper=Inf, 
           control=list(), ...)
}
\arguments{
  \item{formula}{
     This is a modeling formula of the form (as in \code{nls})
           lhsvar ~ rhsexpression
     for example,
           y ~ b1/(1+b2*exp(-b3*tt))
     You may also give this as a string.
}
  \item{start}{
     A named parameter vector. For our example, we could use
        start=c(b1=1, b2=2.345, b3=0.123)
}
  \item{trace}{
     Logical \code{TRUE} if we want intermediate progress to be reported. Default is \code{FALSE}.
} 
  \item{data}{
     A data frame containing the data of the variables in the formula. This data may,
     however, be supplied directly in the parent frame. 
} 
  \item{lower}{
     Lower bounds on the parameters. If a single number, this will be applied to all
     parameters. Default \code{-Inf}.
} 
  \item{upper}{
     Upper bounds on the parameters. If a single number, this will be applied to all
     parameters. Default \code{Inf}.
} 
  \item{control}{ 
     A list of controls for the algorithm. These are as for \code{nlxb()}.

}
  \item{\dots}{
     Any data needed for computation of the residual vector from the expression
       rhsexpression - lhsvar. Note that this is the negative of the usual residual,
     but the sum of squares is the same.
}
}
\details{
    \code{wrapnlsr} first attempts to solve the nonlinear sum of squares problem by using
    \code{nlsmnq}, then takes the parameters from that method to call \code{nls}.
}
\value{
  An object of type nls.}
\author{
John C Nash <nashjc@uottawa.ca>
}
\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}
\examples{
cat("See more examples in nlmrt-package.Rd\n")
cat("kvanderpoel.R test of wrapnlsr\n")
# require(nlsr)
x<-c(1,3,5,7)
y<-c(37.98,11.68,3.65,3.93)
pks28<-data.frame(x=x,y=y)
fit0<-try(nls(y~(a+b*exp(1)^(-c*x)), data=pks28, start=c(a=0,b=1,c=1), 
          trace=TRUE))
print(fit0)
cat("\n\n")
fit1<-nlxb(y~(a+b*exp(-c*x)), data=pks28, start=c(a=0,b=1,c=1), trace = TRUE)
print(fit1) 
cat("\n\nor better\n")
fit2<-wrapnlsr(y~(a+b*exp(-c*x)), data=pks28, start=c(a=0,b=1,c=1), 
              lower=-Inf, upper=Inf, trace = TRUE)

}
\keyword{ nonlinear least squares }

