# NLSProbName: Sacch2_1.R
# NLSProbDescription: {The Sacch2 data frame has 10 rows and 2 columns from an experiment on the pharmacokinetics of
# saccharin.
# The two columns are:This data frame contains the following columns:
# `time`: a numeric vector giving the time since drug administration (min).
# `conc`: a numeric vector giving the observed concentration of saccharin.
# }

# Use the Sacch2 data from NRAIA package

## DATA
time=c(0,   5,  15,  30,  45,  60,  75,  90, 105, 120)
conc = c(0.0, 184.3, 102.0,  50.5,  24.9,  14.1,   8.0,   5.7,   4.0,   2.9)
	
NLSdata <- data.frame(time,conc)

## STARTING VALUE
Dose=1
lKa=13
lKe=17
lCl=7
NLSstart <- c(Dose=Dose,lKa=lKa,lKe=lKe,lCl=lCl) # a starting vector (named!)

## MODEL
NLSformula <-conc ~ Dose*exp(lKe+lKa-lCl)*(exp(-exp(lKe)*time) - exp(-exp(lKa)*time))/(exp(lKa) - exp(lKe))
NLSlower<- c(-Inf,-Inf,-Inf,-Inf)
NLSupper<- c(Inf,Inf,Inf,Inf)
NLSweights <- rep(1,length(time))
NLSsubset <- 1:length(time)
rm(Dose,lKa,lKe,lCl,time,conc)
 # library(nlsr)
 # rjsac<-model2rjfun(NLSformula, data=Sacch2, pvec=NLSstart)
 # tnlxb <- nlxb(NLSformula, data=NLSdata, start=NLSstart, trace=TRUE)
 # ormod<-rjsac(NLSstart)
 # ormod
 # NLSstart
# 
# print(tnlxb)
# tnls<- nls(NLSformula, data=NLSdata, start=NLSstart)
# print(tnls)
# coef(tnlxb)
fm1 <- nls(conc ~ SSfol(1.0, time, lKe, lKa, lCl), data = NLSdata, trace=TRUE)
summary(fm1)
# ## seems like not enough data to fit model
# 
# xyplot(conc ~ time, NLSdata, type = c("g", "b"),
#        xlab = "Time since drug administration (min)",
#        ylab = "Saccharin concentration", aspect = "xy")
# xyplot(conc ~ time, data = NLSdata, type = c("g", "b"),
#        scales = list(y = list(log = 2)), aspect = 'xy',
#        xlab = "Time since drug administration (min)",
#        ylab = "Saccharin concentration")
# ## Not run: These fail because SSfol calls nls() to try to find initial parameters
# ## fm1 <- nls(conc ~ SSfol(1.0, time, lKe, lKa, lCl), data = NLSdata, trace=TRUE)
# ## summary(fm1)
# library(minpack.lm)
# fm1m <- nlsLM(conc ~ SSfol(1.0, time, lKe, lKa, lCl), data = NLSdata, trace=TRUE)
# summary(fm1m)
# xpred <- seq(0, 140, len = 51)
# ypred <- predict(fm1, list(time = xpred, Dose = rep(1.0, length(xpred))))
# lines(xpred, ypred)
#  ss <- function(prm){ 
#       as.numeric(crossprod(NLSdata$conc-rjsac(prm))) 
#        }
#  ss(NLSstart)
#  tnlxb <- nlxb(NLSformula, data=NLSdata, start=NLSstart, 
#                  lower=c(1,5, 5, 7), upper=c(1, 20, 20, 7), trace=TRUE)
# tnlxb 
# ssnm<-function(pnm){ ss(c(1, pnm[1], pnm[2], 7))} 
# tnm<-optim(par=c(13, 17),ssnm, control=list(trace=TRUE))
# tnm
