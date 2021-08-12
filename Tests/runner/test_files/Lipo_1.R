# NLSProbName: Lipo_1.R
# NLSProbDescription: {The Lipo data frame has 12 rows and 2 columns of lipoprotein concentrations over time.
# The two columns are:This data frame contains the following columns:
# `time`: a numeric vector giving the time of the concentration measurement (hr)
# `conc`:  a numeric vector of concentrations.
# }

## DATA
time=c( 0.5,  1.0,  1.5,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0)
conc = c( 46.10, 25.90, 17.00, 12.10,  7.22,  4.51,  3.19,  2.40,  1.82,  1.41,  1.00,
  		    0.94)
NLSdata <- data.frame(time,conc)

## STARTING VALUE
lrc1=1/4
lrc2=-2
A1=100
A2=150
NLSstart <- list(lrc1=lrc1,lrc2=lrc2,A1=A1,A2=A2) # a starting vector (named!)

## MODEL
NLSformula <-conc ~ A1*exp(-exp(lrc1)*time)+A2*exp(-exp(lrc2)*time)
NLSlower<- c(-Inf,-Inf,-Inf,-Inf)
NLSupper<- c(Inf,Inf,Inf,Inf)
NLSsubset <- 1:length(time)
NLSweights <- rep(1,length(time))
rm(time, conc, lrc1,lrc2,A1,A2) 
## tnls <- nls(NLSformula, data=NLSdata, start=NLSstart, trace=TRUE)
## deviance(tnls)
## [1] 0.59319
## summary(tnls)
## Formula: conc ~ A1 * exp(-exp(lrc1) * time) + A2 * exp(-exp(lrc2) * time)
## Parameters:
##   Estimate Std. Error t value Pr(>|t|)    
## lrc1   0.5498     0.0558    9.86  9.5e-06 ***
##   lrc2  -1.0562     0.0628  -16.81  1.6e-07 ***
##   A1    70.7284     1.6791   42.12  1.1e-10 ***
##   A2    19.3882     1.8326   10.58  5.6e-06 ***
## Number of iterations to convergence: 10 
## Achieved convergence tolerance: 4.69e-06
## library(nlsr)
## tnlx <- nlxb(NLSformula, data=NLSdata, start=NLSstart, trace=TRUE)
## tnlx
## residual sumsquares =  0.59319  on  12 observations
## after  14    Jacobian and  15 function evaluations
## name            coeff          SE       tstat      pval      gradient    JSingval   
## lrc1            0.549752       0.05578      9.856  9.459e-06  -7.954e-14       38.44  
## lrc2             -1.0562       0.06284     -16.81   1.59e-07  -2.331e-15       14.68  
## A1               70.7284         1.679      42.12  1.112e-10   8.435e-16       0.163  
## A2               19.3883         1.833      10.58  5.564e-06   1.589e-15      0.1478  
