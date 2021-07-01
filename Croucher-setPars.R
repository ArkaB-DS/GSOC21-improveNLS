rm(list=ls())
# Croucher-setPars.R -- https://walkingrandomly.com/?p=5254
# construct the data vectors using c()
xdata = c(-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9)
ydata = c(0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001)
# some starting values
p1 = 1
p2 = 0.2

Cform <- ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata)
Ccall<-call("-",Cform[[3]], Cform[[2]])
Cstart=list(p1=p1,p2=p2)
Cdata<-data.frame(xdata, ydata)
Cwts <- rep(1, dim(Cdata)[1])
Ctheta<-c("p1","p2")

# ndorig0<-numericDeriv(Ccall, Ctheta)
# ndorig0

library(nlsalt)

# source("nlsModel-minpack.R")

m <- nlsModelx(ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata), data=Cdata, start=Cstart, wts=Cwts)
# setPars <- function(newPars) { # Original -- can we replace??
#   setPars(newPars)
#   resid <<- .swts * (lhs - (rhs <<- getRHS())) # envir = thisEnv {2 x}
#   dev   <<- sum(resid^2) # envir = thisEnv
#   if(length(gr <- attr(rhs, "gradient")) == 1L) gr <- c(gr)
#   QR <<- qr(.swts * gr) # envir = thisEnv
#   (QR$rank < min(dim(QR$qr))) # to catch the singular gradient matrix
# }
# Error: C stack usage  7973524 is too close to the limit ?? 210630
# ?? where are the elements QR, resid, rhs, gr
# <bytecode: 0x55745560c398>
#   <environment: 0x55745792ffd8>

#?? where is useParams, which is supposed to be global.
cat("looking for useParams\n")
if ("useParams" %in% ls(.GlobalEnv)) cat("in .GlobalEnv\n") else cat("not in .GlobalEnv\n")

ls(m)

wts<-Cwts
rhs <- m$getRHS()
lhs <- m$lhs

sP2 <- function(newPars) {
#   sP2(newPars)
#   rhs <- m$getRHS()
#   lhs <- m$lhs
   resid <- wts * (lhs() - rhs())
   dev   <- sum(resid^2) # envir = thisEnv
   gr <- NULL # to define in case no attr.
#   if(length(gr <- attr(rhs, "gradient")) == 1L) gr <- c(gr)
   if(length(gr <- attr(resid, "gradient")) == 1L) gr <- c(gr)
   QR <- qr(wts * gr) # envir = thisEnv
   singular <- (QR$rank < min(dim(QR$qr))) # to catch the singular gradient matrix
   spobj <- list(singular=singular, QR=QR, gr=gr, dev=dev, resid=resid)
}

np1<-c(p1=0.5, p2=0.6)
sp1<-m$setPars(np1)
sp1a<-sP2(np1)

# ?? Note the incredible amount of stuff!!
sv<-m$setVarying()
ls.str(sv)