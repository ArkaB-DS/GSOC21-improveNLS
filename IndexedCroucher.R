# Croucher-base.R -- https://walkingrandomly.com/?p=5254
# construct the data vectors using c()
# NLSProbName: Croucher-base
# NLSProbDescription: {This is a fairly simple 2-parameter problem. This file is the "base"
#     from which others in the Croucher family are built.}
xdata <- c(-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9)
ydata <- c(0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001)
y2 <- 1.45*ydata
x2 <- 1.1*xdata
mm <- length(x2)
idx<-c(rep("ONE",mm), rep("TWO",mm))
xx <- c(xdata, x2)
yy <- c(ydata, y2)
p1<- 1
p2<-0.2
NLSformula <- ydata ~ p1[idx]*cos(p2[idx]*xdata) + p2[idx]*sin(p1[idx]*xdata)
# Ccall<-call("-",Cform[[3]], Cform[[2]])
NLSstart<-list(p1=p1,p2=p2) # This is the default start given in the reference.
NLSdata<-data.frame(xdata=xx, ydata=yy,idx=idx)
## need right form of the start: e.g., start = list(a = rep(b[2], 21), b = rep(b[3], 21), th = b[1]))
itry1<-nls(formula=NLSformula, data=NLSdata, start=c(1, 0.2)) ## fails
