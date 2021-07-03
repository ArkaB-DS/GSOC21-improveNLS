rm(list=ls())
# Croucher-example1.R -- https://walkingrandomly.com/?p=5254
# construct the data vectors using c()
xdata = c(-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9)
ydata = c(0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001)

# look at it
# plot(xdata,ydata)

# some starting values
p1 = 1
p2 = 0.2

# do the fit
fit = nls(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), trace=TRUE)

# summarise
summary(fit)

## Try numericDeriv() function
Cform <- ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata)
Ccall<-call("-",Cform[[3]], Cform[[2]])
Cstart=list(p1=p1,p2=p2)
Cdata<-data.frame(xdata, ydata)
Ctheta<-c("p1","p2")
ndorig0<-numericDeriv(Ccall, Ctheta)
ndorig0

# library(nlsalt)
# ndalt0<-numericDeriv(Ccall, Ctheta)
# print(all.equal(ndorig0, ndalt0))

# retest nls
fit = nls(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), data=Cdata, trace=TRUE)

# summarise

summary(fit)


library(nlsj)
fitj <- nlsj(formula=Cform, start=Cstart, data=Cdata, trace=TRUE)
fitj


# try new all-R form
fit = nlsx(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), trace=TRUE)

# library(nlsr)
# library(minpack.lm)
# print(nlxb(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), trace=TRUE))
# print(nlsLM(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), start=list(p1=p1,p2=p2), trace=TRUE))
