## ExDerivs.R -- a file to explore different aspects of R in calculating Jacobians
## This is an experimental file for learning. 2021-5-25

rm(list=ls()) # clear the workspace (does not completely remove things, need to restart R sometimes)
require(nlsr)
printsum <- function(xx){ print(summary(xx))}
traceval  <-  TRUE  # traceval set TRUE to debug or give full history

# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
            38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)
eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))
str(eunsc)
ceunsc <- " y ~ b1/(1+b2*exp(-b3*tt))"
str(ceunsc)
print(as.formula(ceunsc)==eunsc)

cat("LOCAL DATA IN DATA FRAMES\n")
weeddata1  <-  data.frame(y=ydat, tt=tdat)
weedenv <- list2env(weeddata1)
## weedenv$b1 <- start1$b1
## Error in start1$b1 : $ operator is invalid for atomic vectors
weedenv$b1 <- start1[[1]]
weedenv$b2 <- start1[[2]]
weedenv$b3 <- start1[[3]]
## b1<-1
## b2<-1
## b3<-1
## y<-ydat
## tt<-tdat
## all.vars(eunsc)
rexpr<-call("-",eunsc[[3]], eunsc[[2]])
eval(rexpr, weedenv)
ls.str(weedenv)



names(start1)
start1
str(weeddata1)
ldata<-list2env(as.list(start1),parent=weeddata1)
?list2env
ldata<-list2env(as.list(start1),envir=weeddata1)
weedenv<-list2env(weeddata1)
ldata<-list2env(as.list(start1),envir=weedenv)
ldata
ls.str(ldata)
eval(rexpr,envir=ldata)


mfeunsc<-model.frame(eunsc)
str(mfeunsc)
## 2021-5-25 try derivatives
library(numDeriv)
print(eunsc)
funsc <- model2rjfun(eunsc, start1, data=weeddata1) # from nlsr: creates a function
print(funsc(start1))
print(environment(funsc))
print(ls.str(environment(funsc)))
print(ls(environment(funsc)$data))
eval(eunsc, environment(funsc))
# following gives an error "theta should be of type character"
nDfunsc<-try(numericDeriv(eunsc, start1, weeddata1))
print(nDfunsc)
print(start1)

theta <- c("b1", "b2", "b3")
# Following gives error
## Error in numericDeriv(eunsc, theta, weedenv) : 
## 'language' object cannot be coerced to type 'double'
nDfunsc<-try(numericDeriv(eunsc, theta, weedenv))
print(nDfunsc)
str(theta)
str(weedenv)
print(weedenv)
## try other ways
xeunsc<-expression(eunsc)
## Error in numericDeriv(xeunsc, theta, rho = weedenv) : 
## 'list' object cannot be coerced to type 'double'
nDfunsc<-try(numericDeriv(xeunsc, theta, rho=weedenv))
## Error in numericDeriv(funsc, theta, rho = weedenv) : 
##    cannot coerce type 'closure' to vector of type 'double'
nDfunsc<-try(numericDeriv(funsc, theta, rho=weedenv))

y<-weeddata1$y
tt<-weeddata1$tt
# variables now in local environment
## Error in numericDeriv(funsc, theta) : 
##    cannot coerce type 'closure' to vector of type 'double'
nDfunsc<-try(numericDeriv(funsc, theta))

## Error in numericDeriv(xeunsc, theta, rho = weedenv) : 
## 'list' object cannot be coerced to type 'double'
nDfunsc<-try(numericDeriv(xeunsc, theta))

# Error in formula.default(object, env = baseenv()) : invalid formula
nDfunsc<-try(numericDeriv(as.formula(xeunsc), theta))

formunsc<-formula(eunsc)
str(formunsc)
formunsc[1]
formunsc[2]
formunsc[3]

## Error in numericDeriv(formunsc[3], theta) : object 'b1' not found
nDfunsc<-try(numericDeriv(formunsc[3], theta))

b1<-1
b2<-1
b3<-1
## Error in numericDeriv(formunsc[3], theta) : attempt to apply non-function
nDfunsc<-try(numericDeriv(formunsc[3], theta))


