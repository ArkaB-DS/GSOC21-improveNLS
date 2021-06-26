\name{nlxb}
\alias{nlxb}
\title{Nash variant of Marquardt nonlinear least squares solution via
     qr linear solver.
}
\description{
   Given a nonlinear model expressed as an expression of the form
         \code{lhs ~ formula_for_rhs}
   and a start vector where parameters used in the model formula are named,
   attempts to find the minimum of the residual sum of squares using the
   Nash variant (Nash, 1979) of the Marquardt algorithm, where the linear 
   sub-problem is solved by a qr method. This is a restructured version
   of a function by the same name from package \code{nlmrt} which is now
   deprecated.
}
\usage{
   nlxb(formula, start, trace=FALSE, data, lower=-Inf, upper=Inf,
        masked=NULL, weights=NULL, control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
     This is a modeling formula of the form (as in \code{nls})
           \code{lhsvar ~ rhsexpression}
     for example,
           \code{y ~ b1/(1+b2*exp(-b3*tt))}
     You may also give this as a string.
     Note that the residuals are computed within this code using
           \code{residual <- rhsexpression - lhsvar}
     which is the negative of the usual choice, but the sum of squares is the same.
}
  \item{start}{
     A named parameter vector. For our example, we could use
        \code{start=c(b1=1, b2=2.345, b3=0.123)}
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
  \item{masked}{
     Character vector of quoted parameter names. These parameters will NOT be altered
     by the algorithm. Masks may also be defined by setting lower and upper bounds
     equal for the parameters to be fixed. Note that the starting parameter value
     must also be the same as the lower and upper bound value.
} 
  \item{weights}{
     A vector of fixed weights. The objective function that will be minimized is the
     sum of squares where each residual is multiplied by the square root of the 
     corresponding weight. Default \code{NULL} implies unit weights.
} 

   \item{control}{ 
     A list of controls for the algorithm. These are:
     \describe{
      \item{\code{watch}}{Monitor progress if \code{TRUE}. Default is \code{FALSE}.}
      \item{\code{phi}}{Default is \code{phi=1}, which adds \code{phi*Identity} to 
             Jacobian inner product.}
      \item{\code{lamda}}{Initial Marquardt adjustment (Default \code{0.0001}). 
             Odd spelling is deliberate.}
      \item{\code{offset}}{Shift to test for floating-point equality. 
             Default is \code{100}.}
      \item{\code{laminc}}{Factor to use to increase \code{lamda}. Default is \code{10}.}
      \item{\code{lamdec}}{Factor to use to decrease \code{lamda} is \code{lamdec/laminc}. 
             Default \code{lamdec=4}.}
      \item{\code{femax}}{Maximum function (sum of squares) evaluations. 
             Default is \code{10000}, which is extremely aggressive.}
      \item{\code{jemax}}{Maximum number of Jacobian evaluations. 
             Default is \code{5000}.}
      \item{\code{rofftest}}{Default is \code{TRUE}. Use a termination test of the relative offset 
            orthogonality type. Useful for nonlinear regression problems.}
      \item{\code{smallsstest}}{Default is \code{TRUE}. Exit the function if the sum of squares
            falls below \code{(100 * .Machine$double.eps)^4} times the initial sumsquares. This 
            is a test for a ``small'' sum of squares, but there are problems which are 
            very extreme for which this control needs to be set \code{FALSE}.}
     }
}

}
\details{
    \code{nlxb} attempts to solve the nonlinear sum of squares problem by using
    a variant of Marquardt's approach to stabilizing the Gauss-Newton method using
    the Levenberg-Marquardt adjustment. This is explained in Nash (1979 or 1990) in
    the sections that discuss Algorithm 23. 

    In this code, we solve the (adjusted) Marquardt equations by use of the 
    \code{qr.solve()}. Rather than forming the \verb{J'J + lambda*D} matrix, we augment
    the \code{J} matrix with extra rows and the \code{y} vector with null elements. 
}
\value{
  A list of the following items
  \item{coefficients}{A named vector giving the parameter values at the supposed solution.}
  \item{ssquares}{The sum of squared residuals at this set of parameters.}
  \item{resid}{The residual vector at the returned parameters.}
  \item{jacobian}{The jacobian matrix (partial derivatives of residuals 
           w.r.t. the parameters) at the returned parameters.}
  \item{feval}{The number of residual evaluations (sum of squares computations) used.}
  \item{jeval}{The number of Jacobian evaluations used.}
}
\references{
%% ~put references to the literature/web site here ~
     Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
     Linear Algebra and Function Minimisation._ Adam Hilger./Institute
     of Physics Publications

     others!!

}
\author{
John C Nash <nashjc@uottawa.ca>
}

\seealso{
    Function \code{nls()}, packages \code{\link{optim}} and \code{optimx}.
}

\examples{

cat("See examples in nlsr-package.Rd\n")


}
\keyword{ nonlinear least squares }

