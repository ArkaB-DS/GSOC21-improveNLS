# Croucher-example1.R -- https://walkingrandomly.com/?p=5254

rm(list=ls())

# A function to evaluate f in environment e
# with_env <- function(f, e=parent.frame()) {
#   stopifnot(is.function(f))
#   environment(f) <- e
#   f
# }
# construct the data vectors using c()
xdata = c(-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9)
ydata = c(0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001)
Cform <- ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata)
Ccall<-call("-",Cform[[3]], Cform[[2]])
Cstart<-list(p1=1,p2=0.2)
Cdata<-data.frame(xdata, ydata)
Ctheta<-c("p1","p2")
mdata<-length(xdata)
Cwts <- rep(0.25, mdata)
Csubset<-1:8

library(nlsj)
# 
# print(nlsj.control())
# tmp <- readline("continue?")
# modj <- nlsjModel(form=Cform, data=Cdata, start=Cstart, wts=Cwts, control=nlsj.control())
# ls(modj)
# tmp <- readline("continue?")

# tj = tmod(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=1,p2=.2), subset=1:8, trace=TRUE)

# do the fit 
fitj = nlsjx(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=1,p2=.2), subset=1:8, trace=TRUE)

# summarise
summary(fitj)

fit = nls(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=1,p2=.2), subset=1:8, trace=TRUE)
summary(fit)

# str(modj)
# tmp <- readline("continue?")

# cat("test the model functions:\n")
# modj$resid()
# tmp <- readline("continue?")
# modj$resfun(Cstart)
# tmp <- readline("continue?")
# modj$rjfun(Cstart)
# tmp <- readline("continue?")
# 
# modj$fitted()
# tmp <- readline("continue?")
# modj$formula()
# tmp <- readline("continue?")
# modj$deviance()
# tmp <- readline("continue?")
# modj$lhs()
# tmp <- readline("continue?")
# # modj$jacobian()
# # tmp <- readline("continue?")
# modj$conv()
# tmp <- readline("continue?")
# # modj$incr()
# # tmp <- readline("continue?")
# modj$getPars()
# tmp <- readline("continue?")
# ls(modj$getEnv())
# tmp <- readline("continue?")
# teq <- modj$parset(newPars=c(p1=1, p2=.2))
# teq
# tmp <- readline("continue?")
# neq <- modj$parset(newPars=c(p1=4, p2=3))
# neq
# tmp <- readline("continue?")
# modj$predict(newdata=c(xdata=1, ydata=2))
# tmp <- readline("continue?")
# 


# look at data
# plot(xdata,ydata)

# do the fit (in default R:: nls()
# fit = nls(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=1,p2=.2), trace=TRUE)

# summarise
# summary(fit)
# do the fit
# tmp<-readline("after fit")

# Note: following works OK
# library(nlsr)
# fitn = nlxb(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=1,p2=.2), trace=TRUE)
# fitn
# tmp<-readline("after fitn")
# 
# # summarise
# summary(fitn)

## Try numericDeriv() function
# ndorig0<-numericDeriv(Ccall, Ctheta)
# ndorig0

# tmod <- nlsjModel(form=Cform, data=Cdata, start=Cstart, wts=Cwts)
# tmp<-readline("cont.")
# test deriv
# Ccalld <- deriv(Ccall, Ctheta)
# str(Ccall)
# str(Ccalld)
# ronly <- eval(Ccall)
# ronly
# rj <- eval(Ccalld)
# rj
# 
# library(nlsalt)
# ndalt0<-numericDeriv(Ccall, Ctheta)
# print(all.equal(ndorig0, ndalt0))

## retest nls
# fit = nls(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), data=Cdata, weights=Cwts, trace=TRUE)

## summarise

 # summary(fit)

# library(nlsralt)
# rj <- model2rjfun(Cform, Cstart)
# modelexpr(rj)
# 
# jexpr <- deriv(Cform, names(Cstart))
# jexpr
# 

# library(nlsr)
# library(minpack.lm)
# print(nlxb(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), trace=TRUE))
# print(nlsLM(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), trace=TRUE))

# fitj <- nlsj(formula=Cform, start=Cstart, data=Cdata, trace=TRUE)
# fitj

# library(nlspkg)
# 
# fitsub = nls(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), data=Cdata, subset=1:8, 
#                trace=TRUE)
# 
# # summarise
# 
# summary(fitsub)
# summary(fit)
# 
# callm <- function(subset, form, data, start, wts){ nlsModel(form, data, start, wts)}
# m2 <- callm(subset=1:8, Cform, Cdata, Cstart, Cwts)
# 
# library(nlsralt)
# fitsubx = nlxbx(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), data=Cdata, subset=1:8, 
#              trace=TRUE)
# summary(fitsubx)
# fitsubx$resid # shows subset now registered
# # NOT quite same -- may be analytic derivs etc. -- seems that fitsub is BETTER. WHY?
# fitsubx$coefficients
# fitsub$m$getPars()
# ls(fitsubx)
# fitsubx$jacobian
# fitsub$m$resid()
