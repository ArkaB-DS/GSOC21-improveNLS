# NLSProbName: hobbsx002
# NLSProbDescription: {The Hobbs weed infestation problem to estimate a 
#    3-parameter logistic S curve in its unscaled form from a reasonably
#    easy starting point of (200, 50, 0.3)
#    Using data in workspace.
#  }
## Arkajyoti: Do we want to use names like these so we can find things
## with programs? There may be useful ideas like ROxygen, but I haven't
## used them. 

# Use the Hobbs Weed data
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
y <- weed # need to have data correspond to the formula

NLStestdata <- NULL
NLSstart <- c(b1=200, b2=50, b3=0.3) # a default starting vector (named!)
## Unscaled model
NLSformula <- y ~ b1/(1+b2*exp(-b3*tt))
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- "(formula=NLSformula, start=NLSstart)"
NLSexpout <- "Nonlinear regression model
  model: y ~ b1/(1 + b2 * exp(-b3 * tt))
   data: NLStestdata
     b1      b2      b3 
196.186  49.092   0.314 
 residual sum-of-squares: 2.59

Number of iterations to convergence: 4 
Achieved convergence tolerance: 1.92e-07"

NLSanswer <- "nlsr object: x 
residual sumsquares =  2.5873  on  12 observations
    after  6    Jacobian and  7 function evaluations
  name            coeff          SE       tstat      pval      gradient    JSingval   
b1               196.186         11.31      17.35  3.167e-08  -1.334e-10        1011  
b2               49.0916         1.688      29.08  3.284e-10  -3.589e-09      0.4605  
b3               0.31357      0.006863      45.69  5.768e-12    4.09e-07     0.04714  
"

## Include a line to indicate the end of a problem. This lets us put
## many such files together if we wish.
## --------------------------------------------------------------- ##

