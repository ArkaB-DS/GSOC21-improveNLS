# NLSProbName: Chloride_1.R
# NLSProbDescription: { The Chloride data frame has 54 rows and 2 columns representing measurements of the chloride
# ion concentration in blood cells suspended in a salt solution.
# The two columns are:
# `conc`: A numeric vector giving the chloride ion concentration (%).
# `time`: A numeric vector giving the time of the concentration measurement (min).
# }


# Use the Chloride data from NRAIA package
conc = c(17.3, 17.6, 17.9, 18.3, 18.5, 18.9, 19.0, 19.3, 19.8,
		  19.9, 20.2, 20.5, 20.6, 21.1, 21.5, 21.9, 22.0, 22.3,
		  22.6, 22.8, 23.0, 23.2, 23.4, 23.7, 24.0, 24.2, 24.5,
		  25.0, 25.4, 25.5, 25.9, 25.9, 26.3, 26.2, 26.5, 26.5,
		  26.6, 27.0, 27.0, 27.0, 27.0, 27.3, 27.8, 28.1, 28.1,
		  28.1, 28.4, 28.6, 29.0, 29.2, 29.3, 29.4, 29.4, 29.4)
time = c(2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15, 3.25,
		   3.35, 3.45, 3.55, 3.65, 3.75, 3.85, 3.95, 4.05, 4.15,
		   4.25, 4.35, 4.45, 4.55, 4.65, 4.75, 4.85, 4.95, 5.05,
		   5.15, 5.25, 5.35, 5.45, 5.55, 5.65, 5.75, 5.85, 5.95, 
		   6.05, 6.15, 6.25, 6.35, 6.45, 6.55, 6.65, 6.75, 6.85,
		   6.95, 7.05, 7.15, 7.25, 7.35, 7.45, 7.55, 7.65, 7.75)
## DATA
NLSdata <- data.frame(conc,time)

## STARTING VALUE
Asym = 50
prop=0.6
lrc = log(0.25)
NLSstart <- list(Asym = Asym, prop=prop, lrc = lrc) # a starting vector (named!)

## MODEL
NLSformula <- conc ~ Asym*(1 - prop*exp(-exp(lrc)*time))
NLSlower<- c(-Inf,-Inf,-Inf)
NLSupper<- c(Inf,Inf,Inf)
NLSweights <- rep(1,length(time))
NLSsubset <- 1:length(time)
rm(conc,time,Asym,prop,lrc)