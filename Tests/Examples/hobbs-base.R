# NLSProbName: hobbs-base
# NLSProbDescription: {The Hobbs weed infestation problem to estimate a 
#    3-parameter logistic S curve in its unscaled form from a reasonably
#    easy starting point of (200, 50, 0.3)
#  }
## Arkajyoti: Do we want to use names like these so we can find things
## with programs? There may be useful ideas like ROxygen, but I haven't
## used them. 

# Use the Hobbs Weed data
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12

NLStestdata <- data.frame(y=weed, tt=tt) # should we use standard name?
NLSstart <- c(b1=200, b2=50, b3=0.3) # an easy starting vector (named!)
## Unscaled model
NLSformula <- y ~ b1/(1+b2*exp(-b3*tt))
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- "(formula=NLSformula, data=NLStestdata, start=NLSstart)"
