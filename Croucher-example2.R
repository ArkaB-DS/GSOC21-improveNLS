# Croucher-example2.R -- https://walkingrandomly.com/?p=5254

rm(list=ls())

# A function to evaluate f in environment e
# with_env <- function(f, e=parent.frame()) {
#   stopifnot(is.function(f))
#   environment(f) <- e
#   f
# }
# construct the data vectors using c()
source("../Tests/Croucher-base.R")
library(nlsj)
#
# Neither of these work -- good that they match up
#fit0 <- nls(formula=Cform)
#fit0j<-nlsj(formula=Cform)

fitj1 = nlsj(formula=NLSformula, data=NLStestdata, start=NLSstart, control=nlsj.control(derivmeth="default"))
summary(fitj1)
fitj1n = nlsj(formula=NLSformula, data=NLStestdata, start=NLSstart, control=nlsj.control(derivmeth="numericDeriv"))
summary(fitj1n)
fitx1 = nlsr::nlxb(formula=NLSformula, data=NLStestdata, start=NLSstart)
fitx1 # ?? note that nlsr needs fix for summary() and print()
fitn1 <- nls(formula=NLSformula, data=NLStestdata, start=NLSstart)
summary(fitn1)
rj1<-fitj1$m$resid() # ?? wrong sign in resid part
print(rj1)
rn1<-fitn1$m$resid()
print(rn1)
options(digits=16) # to compare coefficiencts
fitj1$m$getPars()
fitj1n$m$getPars()
fitn1$m$getPars()
fitx1$coefficients


tmp <- readline("more?")
fit = nls(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=1,p2=.2), subset=1:8, trace=TRUE)
summary(fit)

# NOTE: nlsj() zeros subsetted residual items, but includes them in the outut of resid()
rj1<-fitj$m$resid() # ?? wrong sign in resid part
print(rj1)
rn1<-fit$m$resid()
print(rn1)

subres <- function(robj, subset) {## Subset a resid() object using subset argument
    J <- attr(robj, "gradient")  
    sres <- robj[subset]
    sJ <- J[subset, ]
    attr(sres, "gradient") <- sJ
    sres  
}

rescomp <- function(r1, r2){ # Compare two residual objects for nls() class outputs
   rdiff <- r1 -r2
   jdiff <- attr(r1,"gradient") - attr(r2,"gradient")
   mr <- max(abs(rdiff))
   mj <- max(abs(jdiff))
   cat("Maximum absolute differences in residual =",mr," in jacobian=",mj,"\n")
   rdiff
}

r1<-subres(r1,Csubset)
rescomp(r1,r2)
as.numeric(r1)
as.numeric(r2)
j1 <-attr(r1,"gradient")
j2 <-attr(r2,"gradient")
j1
j2

fitjn = nlsj(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=1,p2=.2), subset=1:8, control=nlsj.control(derivmeth="numericDeriv"), trace=TRUE)

# summarise
summary(fitjn)
r1n<-fitjn$m$resid() # ?? wrong sign in resid part

r1n<-subres(r1n,Csubset)
rescomp(r1n,r2)
as.numeric(r1n)
as.numeric(r2)
j1n <-attr(r1n,"gradient")
j2 <-attr(r2,"gradient")
j1n
j2



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
