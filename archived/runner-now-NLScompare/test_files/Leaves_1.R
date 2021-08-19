# NLSProbName: Leaves_1.R
# NLSProbDescription: { The Leaves data frame has 15 rows and 2 columns of leaf length over time.
# The two columns are:This data frame contains the following columns:
# `time`: e time from initial emergence (days).
# `length`: leaf length (cm).
# }

## DATA
time=c( 0.5,  1.5,  2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5, 11.5, 12.5, 13.5,
		  14.5)
length = c( 1.3,  1.3,  1.9,  3.4,  5.3,  7.1, 10.6, 16.0, 16.4, 18.3, 20.9, 20.5, 21.3, 21.2,
			20.9)
NLSdata <- data.frame(time,length)

## STARTING VALUE
Asym=3
xmid=2
scal=1
NLSstart <-list(Asym=Asym,xmid=xmid,scal=scal) # a starting vector (named!)
tnls2 <- nls(NLSform2, data=NLSdata)
## MODEL
NLSformula <-length ~ Asym/(1+exp((xmid-time)/scal))
NLSform2<-length ~ SSlogis(time, Asym, xmid, scal)
NLSlower<- c(-Inf,-Inf,-Inf)
NLSupper<- c(Inf,Inf,Inf)
NLSsubset <- 1:length(time)
NLSweights <- rep(1,length(time))
rm(time,length,Asym,xmid,scal)
## tnls <- nls(NLSformula, NLSstart, NLSdata, trace=TRUE)
## 2280.1    (NA): par = (0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5 1.3 1.3 1.9 3.4 5.3 7.1 10.6 16 16.4 18.3 20.9 20.5 21.3 21.2 20.9)
## Error in qr.default(.swts * gr) : 
##  NA/NaN/Inf in foreign function call (arg 1)
## tnls2<-nls(NLSform2, data=NLSdata) # selfStart SSlogis
## model: length ~ SSlogis(time, Asym, xmid, scal)
## data: NLSdata
## Asym  xmid  scal 
## 21.51  6.36  1.61 
## residual sum-of-squares: 6.21
## Number of iterations to convergence: 0 
## Achieved convergence tolerance: 4.38e-06
## ----- nlxb 
## tnlx <- nlxb(NLSformula, NLSstart, NLSdata, trace=TRUE)
## Start:lamda: 1e-04  SS= 2280.1  at  Asym = 3  xmid = 2  scal = 1  1 / 0
## <<lamda: 4e-05  SS= 1190.5  at  Asym = 18.16  xmid = 14.889  scal = 11.56  2 / 1
## <<lamda: 2.6844e-05  SS= 6.2102  at  Asym = 21.509  xmid = 6.3604  scal = 1.6072  20 / 14
## print(tnlx)
## residual sumsquares =  6.2102  on  15 observations
## after  14    Jacobian and  20 function evaluations
## name            coeff          SE       tstat      pval      gradient    JSingval   
## Asym             21.5089        0.4154      51.78  1.767e-15   1.576e-11       7.937  
## xmid             6.36035        0.1388      45.82  7.629e-15  -1.204e-11       7.072  
## scal             1.60724        0.1152      13.95  8.888e-09   9.966e-11       1.666
