## Trun1.R -- computes jacobian of Hobbs at 1,1,1 using numericDeriv()
# rm(list=ls()) # clear the workspace (does not completely remove things, need to restart R sometimes)

printsum <- function(xx){ print(summary(xx))}
traceval  <-  TRUE  # traceval set TRUE to debug or give full history

# Data for Hobbs problem
ydat  <-  c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443, 
            38.558, 50.156, 62.948, 75.995, 91.972) # for testing
tdat  <-  seq_along(ydat) # for testing

# A simple starting vector -- must have named parameters for nlxb, nls, wrapnlsr.
start1  <-  c(b1=1, b2=1, b3=1)
eunsc  <-   y ~ b1/(1+b2*exp(-b3*tt))
weeddata1  <-  data.frame(y=ydat, tt=tdat)
weedenv <- list2env(weeddata1)
weedenv$b1 <- start1[[1]]
weedenv$b2 <- start1[[2]]
weedenv$b3 <- start1[[3]]
rexpr<-call("-",eunsc[[3]], eunsc[[2]])
r0<-eval(rexpr, weedenv)
cat("Sumsquares at 1,1,1 is ",sum(r0^2),"\n")

theta <- c("b1", "b2", "b3")
ndeunsc<-numericDeriv(rexpr, theta, rho=weedenv)
print(ndeunsc)
print(sum(ndeunsc^2))
