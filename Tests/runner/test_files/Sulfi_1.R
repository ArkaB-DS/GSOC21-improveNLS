# NLSProbName: Sulfi_1.R
# NLSProbDescription: {The Sulfi data frame has 12 rows and 2 columns from an experiment on the pharmacokinetics of
# sulfisoxazole.
# The two columns are:This data frame contains the following columns:
# `time`:  a numeric vector giving the time since drug administration (min).
# `conc`:  a numeric vector giving the observed concentration of sulfisoxazole (µg/ml).

# }


# Use the Sulfi data from NRAIA package

## DATA
time=c( 0.25,  0.50,  0.75,  1.00,  1.50,  2.00,  3.00,  4.00,  6.00, 12.00, 24.00,
		  48.00)
conc = c( 215.6, 189.2, 176.0, 162.8, 138.6, 121.0, 101.2,  88.0,  61.6,  22.0,   4.4,
		    0.1)
NLSdata <- data.frame(time,conc)

## STARTING VALUE
lrc1=1
lrc2=-2
A1=19
A2=31
NLSstart <-c(lrc1=lrc1,lrc2=lrc2,A1=A1,A2=A2) # a starting vector (named!)

## MODEL
NLSformula <-conc ~ A1*exp(-exp(lrc1)*time)+A2*exp(-exp(lrc2)*time)
NLSlower<- c(-Inf,-Inf,-Inf,-Inf)
NLSupper<- c(Inf,Inf,Inf,Inf)
NLSweights <- rep(1,length(time))
NLSsubset <- 1:length(time)
rm(time,conc,lrc1,lrc2,A1,A2)
