# File: badJlogmod.R
# A problem illustrating poor numeric Jacobian
form<-y ~ 10*a*(8*b-log(0.075*c*x)) # the model formula
# This model uses log near a small argument, which skirts the dangerous
# value of 0. The parameters a, b, c could all be 1 "safely" as a start.
x<-1:20 # define x
a<-1.01
b<-.9
c<-.95
y <- 10*a*(8*b-log(0.075*c*x))+0.2*runif(20) # compute a y
df<-data.frame(x=x, y=y)
# plot(x,y) # for information
st<-c(a=1, b=1,c=1) # set the "default" starting vector
n0<-try(nls(form, start=st, data=df)) # and watch the fun as this fails. 
library(nlsr) # but this will work
n1<-nlxb(form, start=st, data=df)
n1
jmod<-model2rjfun(form, pvec=st,data=data.frame(x=x, y=y)) # extract the model
Jatst<-jmod(st) # compute this at the start from package nlsr
Jatst<-attr(Jatst,"gradient") # and extract the Jacobian
#
# Now try to compute Jacobian produced by nls()
theta <- c("a","b","c") # the parameter names
env<-environment(form) # We need the environment of the formula
eform<-eval(form, envir=env) # and the evaluated expression
jnlsatst<-numericDeriv(form[[3L]], theta=theta, rho=env)
Jnls<-attr(jnlsatst,"gradient")
Jnls # from nls()
Jatst # from nlsr -- analytic derivative
svd(Jnls)$d
svd(Jatst)$d
# what if we try to get closer to the parameters?
n0b<-try(nls(form, start=c(a=1, b=1,c=1.3), data=data.frame(x=x, y=y), trace=TRUE))
# or even start at the solution?
n0c<-try(nls(form, start=coef(n1), data=data.frame(x=x, y=y), trace=TRUE))
