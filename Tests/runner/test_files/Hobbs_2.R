# NLSProbName: Hobbs_2.R
# NLSProbDescription: {The Hobbs weed infestation problem to estimate a 
#    3-parameter logistic S curve in its unscaled form from a  
# starting point of (1,1,.1) and lower=(0,0,0) and upper=(500,100,5)
# }
# feasible bounds containing solution

NLSstart <- list(b1=1,b2=1,b3=.1)
NLSlower <- c(0,0,0)
NLSupper <- c(b1=500,b2=100,b3=5)
# residual sumsquares =  2.5873  on  12 observations
# #  coef(bestsol) = c(b1=196.18626,  b2=49.09164, b3= 0.31357)
