# NLSProbName: Hobbs_4.R
# NLSProbDescription: {The Hobbs weed infestation problem to estimate a 
#    3-parameter logistic S curve in its unscaled form from a  
# starting point of (1,1,.1) and lower=(0,0,0) and upper=(500,100,.25)
# }
# bound on b3 is lower than unconstrained solution value.

NLSstart <- list(b1=1,b2=1,b3=.1) # feasible start
NLSlower <- c(0,0,0)
NLSupper <- c(b1=500,b2=100,b3=.25)
# residual sumsquares =  29.993  on  12 observations
# coef(bestsol) = c(b1=500.000, b2=87.942, b3=0.250)
# Note: nls(algorithm="port") gives residual sumsquares = 30, at 500.00  87.94   0.25 
