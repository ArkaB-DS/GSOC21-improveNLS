nlsderivchk <- function(residexpr, pnames){
# residexpr is the expression we are going to differentiate for jacobian
# pnames is the vector of parameter names
  nnames <- length(pnames)
  dOK <- rep(FALSE, nnames) # start with "bad" derivatives
  for (i in 1:nnames){
     onename <- pnames[i]
     aderiv <- try(deriv(residexpr, onename))
     if (! inherits(aderiv, "try-error") ) dOK[i] <- TRUE
  }
  dOK
}

### A test of the nlsderivchk function


form2expr <- function(modelformula){
    if (length(modelformula) == 2) {
        residexpr <- modelformula[[2]]
    } else if (length(modelformula) == 3) {
        residexpr <- call("-", modelformula[[3]], modelformula[[2]])
    } else stop("Unrecognized formula")
    residexpr
}

# cat("Croucher example\n")
# Cform <- ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata)
# Ctheta<-c("p1","p2")
# Cexpr <- form2expr(Cform)
# tcrouch <- nlsderivchk(Cexpr, Ctheta)
# tcrouch
# 
# cat("Tanh example\n")
# Tform <- ydata ~ a1 * tanh(a2*tdata - a3)
# Ttheta <- c("a1", "a2", "a3")
# Texpr <- form2expr(Tform)
# Ttanh <- nlsderivchk(Texpr, Ttheta)
# Ttanh
# 
# # This one fails. And deriv w.r.t. p1 and p4 should be possible analytically.
# # ?? need to see why we can't do it.
# cat("TDist example\n")
# Dform<- ydata ~ p1*xx + p4 * dt((xx - p2), df=10, ncp=p3)
# Dtheta <- c("p1","p2", "p3", "p4")
# Dexpr <- form2expr(Dform)
# TD <- nlsderivchk(Dexpr, Dtheta)
# TD
# 
