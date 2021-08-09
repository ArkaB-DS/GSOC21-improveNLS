# NLSProbName: Rumford_1.R
# NLSProbDescription: {The Rumford data frame has 13 rows and 2 columns from an experiment by Count Rumford on the
# rate of cooling.
# The two columns are:This data frame contains the following columns:
# `time`:  a numeric vector giving the time since the beginning of the experiment (hr)
# `temp`:  a numeric vector giving the temperature (degrees Fahrenheit) of the cannon.

# }


# Use the Rumford data from NRAIA package
## DATA
time <- c(4.0,  5.0,  7.0, 12.0, 14.0, 16.0, 20.0, 24.0, 28.0, 31.0, 34.0, 37.5, 41.0)
temp <- c(126, 125, 123, 120, 119, 118, 116, 115, 114, 113, 112, 111, 110)
	
NLSdata <- data.frame(time,temp)

## STARTING VALUE
tc = -0.01
NLSstart <- list(tc = tc) # a starting vector (named!)

## MODEL
NLSformula <-temp ~ 60 + 70 * exp(tc * time)
NLSlower<- c(-Inf)
NLSupper<- c(Inf)
NLSweights <- rep(1, length(time))
NLSsubset<-1:length(time)
rm(tc,time, temp)