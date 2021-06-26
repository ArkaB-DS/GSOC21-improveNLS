\name{nlfb}
\alias{nlfb}
\title{Nash variant of Marquardt nonlinear least squares solution via
     qr linear solver.
}
\description{
   Given a nonlinear model expressed as a vector valued residual function
   \code{resfn} and a start vector of parameter values for that function, 
   attempts to find the minimum of the residual sum of squares using the
   Nash variant (Nash, 1979) of the Marquardt algorithm, where the linear 
   sub-problem is solved by a qr method. This is a restructured version
   of a function by the same name from package \code{nlmrt} which is now
   deprecated.
}
\usage{
   nlfb(start, resfn, jacfn=NULL, trace=FALSE, lower=-Inf, upper=Inf, 
         maskidx=NULL, weights=NULL, data=NULL, control, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resfn}{
     A function that evaluates the residual vector for computing the elements of
     the sum of squares function at the set of parameters \code{start}. Where this
     function is created by actions on a formula or expression in \code{nlxb}, this
     residual vector will be created by evaluation of the 'model - data', rather than
     the conventional 'data - model' approach. The sum of squares is the same.
}
  \item{jacfn}{
     A function that evaluates the Jacobian of the sum of squares function, that
     is, the matrix of partial derivatives of the residuals with respect to each
     of the parameters. If NULL (default), uses an approximation.

     The Jacobian MUST be returned as the attribute "gradient" of this function,
     allowing \code{jacfn} to have the same name and be the same code block
     as \code{resfn}, which may permit some efficiencies of computation.
}
  \item{start}{
     A named parameter vector. For our example, we could use
        start=c(b1=1, b2=2.345, b3=0.123)
     \code{nls()} takes a list, and that is permitted here also.
}
  \item{trace}{
     Logical TRUE if we want intermediate progress to be reported. Default is FALSE.
} 
  \item{lower}{
     Lower bounds on the parameters. If a single number, this will be applied to all
     parameters. Default -Inf.
} 
  \item{upper}{
     Upper bounds on the parameters. If a single number, this will be applied to all
     parameters. Default Inf.
} 
  \item{maskidx}{
     Vector if indices of the parameters to be masked. These parameters will NOT be altered
     by the algorithm. Note that the mechanism here is different from that in \code{nlxb}
     which uses the names of the parameters.

} 

  \item{weights}{
     A vector of fixed weights. The objective function that will be minimized is the
     sum of squares where each residual is multiplied by the square root of the 
     corresponding weight.  Default (NULL) implies unit weights.
} 

  \item{data}{
     Data frame of variables used by resfn and jacfn to compute the
     required residuals and Jacobian.
}

  \item{control}{ 
     A list of controls for the algorithm. These are:
     \describe{
      \item{\code{watch}}{Monitor progress if TRUE. Default is FALSE.}
      \item{\code{phi}}{Default is phi=1, which adds phi*Identity to Jacobian inner product.}
      \item{\code{lamda}}{Initial Marquardt adjustment (Default 0.0001). Odd spelling is deliberate.}
      \item{\code{offset}}{Shift to test for floating-point equality. Default is 100.}
      \item{\code{laminc}}{Factor to use to increase lamda. Default is 10.}
      \item{\code{lamdec}}{Factor to use to decrease lamda is lamdec/laminc. Default lamdec=4.}
      \item{\code{femax}}{Maximum function (sum of squares) evaluations. Default is 10000,
            which is extremely aggressive.}
      \item{\code{jemax}}{Maximum number of Jacobian evaluations. Default is 5000.}
      \item{\code{ndstep}}{Stepsize to use to computer numerical Jacobian approximatin. Default
            is 1e-7.}
      \item{\code{rofftest}}{Default is TRUE. Use a termination test of the relative offset 
            orthogonality type. Useful for nonlinear regression problems.}
      \item{\code{smallsstest}}{Default is TRUE. Exit the function if the sum of squares falls
            below (100 * .Machine$double.eps)^4 times the initial sumsquares. This is a test
            for a ``small'' sum of squares, but there are problems which are very extreme 
            for which this control needs to be set FALSE.}
     }
}

  \item{\dots}{
     Any data needed for computation of the residual vector from the expression
       rhsexpression - lhsvar. Note that this is the negative of the usual residual,
     but the sum of squares is the same. It is not clear how the dot variables should
     be used, since data should be in 'data'.
}
}
\details{
    \code{nlfb} attempts to solve the nonlinear sum of squares problem by using
    a variant of Marquardt's approach to stabilizing the Gauss-Newton method using
    the Levenberg-Marquardt adjustment. This is explained in Nash (1979 or 1990) in
    the sections that discuss Algorithm 23.

    In this code, we solve the (adjusted) Marquardt equations by use of the 
    \code{qr.solve()}. Rather than forming the J'J + lambda*D matrix, we augment
    the J matrix with extra rows and the y vector with null elements. 
}
\value{
  A list of the following items
  \item{coefficients}{A named vector giving the parameter values at the supposed solution.}
  \item{ssquares}{The sum of squared residuals at this set of parameters.}
  \item{resid}{The residual vector at the returned parameters.}
  \item{jacobian}{The jacobian matrix (partial derivatives of residuals w.r.t. the parameters)
           at the returned parameters.}
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
cat("See examples in nls-package.Rd\n")

}
\keyword{ nonlinear least squares }

