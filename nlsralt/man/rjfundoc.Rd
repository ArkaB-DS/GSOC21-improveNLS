\name{rjfundoc}
\alias{rjfundoc}
\alias{print.rjfundoc}
\title{Document functions that evaluate residuals at 
particular parameter locations.
}
\description{
Output a description of the model and the data that was used at the time
of its creation to the console and optionally to a file. The purpose of 
this function is to provide a record of the details underlying the function
\code{fun} as well as to aid users wishing to create objective functions for
optimization tools.
}
\usage{
rjfundoc(fun, savefile=NULL)
\method{print}{rjfundoc}(x, ...)
}
\arguments{
  \item{fun}{
A function produced by \code{model2rjfun}.
}

  \item{savefile}{A character string or connection giving a location in which
to record the \code{rjfundoc} output. 
}
  \item{x, ...}{Object to print, and other arguments that will be ignored.
}

}
\details{
\code{rjfundoc} displays the contents of the environment associated with
\code{fun}. 
}

\value{
\code{rjfundoc} returns a list of class \code{"rjfundoc"} containing values
extracted from the environment of \code{fun}.
}

\author{
John Nash and Duncan Murdoch
}

\seealso{
\code{\link{nls}}
}
\examples{

## require(nlsr)
traceval  <-  TRUE  # traceval set TRUE to debug or give full history
# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
          38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing
# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)
weeddata2  <-  data.frame(y=1.5*ydat, tt=tdat)
escal <-  y ~ ms*b1/(1+b2*exp(-b3*tt))
ms <- 1
weeddata3<-weeddata2
weeddata3$y<-0
grs3<-model2rjfun(escal, start1, data=weeddata3, ms=ms)
res3<-grs3(start1)
res3
rjfundoc(grs3)

# Now a different value of ms
ms<-2
grs3b<-model2rjfun(escal, start1, data=weeddata3, ms=ms)
res3b<-grs3b(start1)
res3b
rjfundoc(grs3b)
# rjfundoc(grs3b, savefile="grs3save.txt") ## to save the output

}
\keyword{ nonlinear }
