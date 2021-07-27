# File: badJ2.R
# A problem illustrating poor numeric Jacobian
form<-y ~ 10*a*(8+b*log(1-0.049*c*x)) # the model formula
# This model uses log near a small argument, which skirts the dangerous
# value of 0. The parameters a, b, c could all be 1 "safely" as a start.
x<-3*(1:10) # define x
np<-length(x)
a<-1.01
b<-.9
eps<-1e-6
c<-1/(max(x)*.049)-eps
cat("c =",c,"\n")
y <- 10*a*(8+b*log(1-0.049*c*x))+0.2*runif(np) # compute a y
df<-data.frame(x=x, y=y)
plot(x,y) # for information
st<-c(a=1, b=1,c=c) # set the "default" starting vector
n0<-try(nls(form, start=st, data=df)) # and watch the fun as this fails. 
summary(n0)
library(nlsr) # but this will work
n1<-nlxb(form, start=st, data=df)
n1
coef(n1)-coef(n0)
jmod<-model2rjfun(form, pvec=st,data=data.frame(x=x, y=y)) # extract the model
Jatst<-jmod(st) # compute this at the start from package nlsr
Jatst<-attr(Jatst,"gradient") # and extract the Jacobian
#
# Now try to compute Jacobian produced by nls()
env<-environment(form) # We need the environment of the formula
eform<-eval(form, envir=env) # and the evaluated expression
localdata<-list2env(as.list(st), parent=env)
jnlsatst<-numericDeriv(form[[3L]], theta=names(st), rho=localdata)
Jnls<-attr(jnlsatst,"gradient")
Jnls # from nls()
Jatst # from nlsr -- analytic derivative
max(abs(Jnls-Jatst))
svd(Jnls)$d
svd(Jatst)$d
# Even start at the solution?
n0c<-try(nls(form, start=coef(n1), data=data.frame(x=x, y=y), trace=TRUE))
summary(n0c)
## attempts with nlsj 
library(nlsj)
n0jn<-try(nlsj(form, start=st, data=df, trace=TRUE,control=nlsj.control(derivmeth="numericDeriv")))
summary(n0jn)
# tmp<-readline("more.")
n0ja<-try(nlsj(form, start=st, data=df, trace=TRUE,control=nlsj.control(derivmeth="default"))) 
summary(n0ja)
