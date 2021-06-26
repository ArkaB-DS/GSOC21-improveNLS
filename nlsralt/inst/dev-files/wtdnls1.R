DNase1 <- subset(DNase, Run == 1)
fm3DNase1 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
                 data = DNase1,
                 start = list(Asym = 3, xmid = 0, scal = 1))

FITTED <- unique(fitted(fm3DNase1))
DAT <- sapply(FITTED, function(x) rnorm(20, mean = x, sd = 0.02 * x))
matplot(t(DAT), type = "p", pch = 16, lty = 1, col = 1)
lines(FITTED, col = 2)
CONC <- unique(DNase1$conc)
fitDAT <- data.frame(conc = rep(CONC, each = 20), density = matrix(DAT))

## First we create the unweighted fit:
  
  FIT1 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
              data = fitDAT,
              start = list(Asym = 3, xmid = 0, scal = 1))

## Then we fit the data with weights w = 1/var(y). 
##  IMPORTANT: we need to replicate the weight values by 20 in order to match the data length.

VAR <- tapply(fitDAT$density, fitDAT$conc, var)
VAR <- rep(VAR, each = 20)
FIT2 <- nls(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
            data = fitDAT, weights = 1/VAR,
            start = list(Asym = 3, xmid = 0, scal = 1))

## For calculation of \chi^2_\nu and its corresponding p-value, 
## we use the fitchisq function of my ‘qpcR’ package:
  
library(qpcR)
fitchisq(FIT1)
fitchisq(FIT2)

library(nlsr)
fitDAT$wts <- 1/VAR
cat("nfit1 - no weights\n")
nfit1 <-  nlxb(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
              data = fitDAT, start = list(Asym = 3, xmid = 0, scal = 1))
cat("nfit2 - weights\n")
nfit2 <-  nlxb(density ~ Asym/(1 + exp((xmid - log(conc))/scal)),
               data = fitDAT, start = list(Asym = 3, xmid = 0, scal = 1),
               weights=fitDAT$wts)

print(nfit1)
print(FIT1)
print(nfit2)
print(FIT2)
