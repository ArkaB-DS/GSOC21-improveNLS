library(Deriv)
rm(x) # ensures x is undefined
Deriv(~ x, "x")  # returns [1] 1 -- clearly a bug!
Deriv(~ x^2, "x")   # returns 2 * x
x <- quote(x^2)
Deriv(x, "x") # returns 2 * x
