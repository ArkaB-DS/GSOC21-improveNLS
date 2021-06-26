Treated <- Puromycin[Puromycin$state == "treated", ]
weighted.MM <- function(resp, conc, Vm, K)
{
  ## Purpose: exactly as white book p. 451 -- RHS for nls()
  ##  Weighted version of Michaelis-Menten model
  ## ----------------------------------------------------------
  ## Arguments: 'y', 'x' and the two parameters (see book)
  ## ----------------------------------------------------------
  ## Author: Martin Maechler, Date: 23 Mar 2001
  
  pred <- (Vm * conc)/(K + conc)
  (resp - pred) / sqrt(pred)
}
## Here is the estimation using predicted value weights
Pur.wtpred <- nls( ~ weighted.MM(rate, conc, Vm, K), data = Treated,
                   start = list(Vm = 200, K = 0.1))
summary(Pur.wtpred)

## nlxb cannog use this form
Pur.wtxpred <- try(nlxb( ~ weighted.MM(rate, conc, Vm, K), data = Treated,
                         start = list(Vm = 200, K = 0.1)))
## and the structure is wrong for nlfb. See wres.MM below.
Pur.wtfpred <- try(nlfb(resfn=weighted.MM(rate, conc, Vm, K), data = Treated,
                        start = list(Vm = 200, K = 0.1)))



## Now use fixed weights and the weights argument in nls()
Pur.wtw <- nls( rate ~ (Vm * conc)/(K + conc), data = Treated,
                start = list(Vm = 200, K = 0.1), weights=1/rate)
summary(Pur.wtw)

## nlsr::nlxb should give essentially the same result
require(nlsr)
Pur.wtx <- nlxb( rate ~ (Vm * conc)/(K + conc), data = Treated,
                 start = list(Vm = 200, K = 0.1), weights=1/Treated$rate)
summary(Pur.wtx)
## check difference
coef(Pur.wtw) - coef(Pur.wtx)

## Try using a function, that is nlsr::nlfb
stw <- c(Vm = 200, K = 0.1)
wres.MM <- function(prm, rate, conc) {
  Vm <- prm[1]
  K <- prm[2]
  pred <- (Vm * conc)/(K + conc)
  (rate - pred) / sqrt(rate) # NOTE: NOT pred here
}

# test function first to see it works
print(wres.MM(stw, rate=Treated$rate, conc=Treated$conc))

Pur.wtf <- nlfb(start=stw,  resfn=wres.MM, rate = Treated$rate, conc=Treated$conc, 
                trace=FALSE) # Note: weights are inside function already here
summary(Pur.wtf)
print(Pur.wtf)
## check estimates with nls() result
coef(Pur.wtw) - coef(Pur.wtf)
## and with nlxb result
coef(Pur.wtx) - coef(Pur.wtf)

wres0.MM <- function(prm, rate, conc) {
  Vm <- prm[1]
  K <- prm[2]
  pred <- (Vm * conc)/(K + conc)
  (rate - pred) # NOTE: NO weights
}
Pur.wtf0 <- nlfb(start=stw,  resfn=wres0.MM, rate = Treated$rate, conc=Treated$conc, 
                 weights=1/Treated$rate, trace=FALSE) # Note: explicit weights
summary(Pur.wtf0)
print(Pur.wtf0)
## check estimates with nls() result
coef(Pur.wtw) - coef(Pur.wtf0)
## and with nlxb result
coef(Pur.wtx) - coef(Pur.wtf0)

wss.MM <- function(prm, rate, conc) {
  ss <- as.numeric(crossprod(wres.MM(prm, rate, conc)))
}
library(optimr)
osol <- optimr(stw, fn=wss.MM, gr="grnd", method="Nelder-Mead", control=list(trace=0), 
               rate=Treated$rate, conc=Treated$conc)
print(osol)
## difference from nlxb
osol$par - coef(Pur.wtx)


