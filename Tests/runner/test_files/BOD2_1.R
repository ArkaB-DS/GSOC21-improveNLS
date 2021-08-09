# NLSProbName: BOD2_1.R
# NLSProbDescription: { The BOD2 data frame has 8 rows and 2 columns
# giving the biochemical oxygen demand versus time in an evaluation of water quality.
# The two columns are:
# `time`: A numeric vector giving the time of the measurement (days).
# `demand`: A numeric vector giving the biochemical oxygen demand (mg/l).
# }


# Use the BOD2 data from NRAIA package
demand <- c(0.47, 0.74, 1.17, 1.42, 1.60, 1.84, 2.19, 2.17)
time <- c( 1,  2,  3,  4,  5,  7,  9, 11)

## DATA 
NLSdata <- data.frame(demand,time) 

## STARTING VALUE
A = 2.2
lrc = log(0.25)
NLSstart <- list(A = A, lrc = lrc) # a starting vector (named!)

## MODEL
NLSformula <- demand ~ A*(1-exp(-exp(lrc)*time))
NLSlower<- c(-Inf,-Inf)
NLSupper<- c(Inf,Inf)
NLSweights <- rep(1,length(time)) ## ?? Find better weights
NLSsubset <- 1:length(time)
rm(demand,time,A,lrc) 