# NLSProbName: Tetra_1.R
# NLSProbDescription: {The Tetra data frame has 9 rows and 2 columns from an experiment on the pharmacokinetics of
# tetracycline.
# The two columns are:This data frame contains the following columns:
# `time`:  a numeric vector of time since drug administration (hr).
# `conc`:  a numeric vector of tetracycline concentrations.

# }


# Use the Tetra data from NRAIA package

## DATA
time=c( 1,  2,  3,  4,  6 , 8, 10, 12, 16)
conc = c( 0.7, 1.2, 1.4, 1.4, 1.1, 0.8, 0.6, 0.5, 0.3)
NLStestdata <- data.frame(time,conc)

## STARTING VALUE
lrc1=-2
lrc2=1/4
A1=150
A2=50
NLSstart <-c(lrc1=lrc1,lrc2=lrc2,A1=A1,A2=A2) # a starting vector (named!)

## MODEL
NLSformula <- conc ~ A1*exp(-exp(lrc1)*time)+A2*exp(-exp(lrc2)*time)
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- "(formula=NLSformula, data=NLStestdata, start=NLSstart)"

test_that("This is a singular gradient problem",
	expect_error(eval(parse(text=paste("nls",NLSrunline))),
		regex="singular gradient",ignore.case=TRUE)
	)

# this is an improvement over nls
output_nlsj <- eval(parse(text=paste("nlsj::nlsj",NLSrunline))) # nlsj is the new nls
summary(output_nlsj)
