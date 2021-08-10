# NLSProbName: Chloride_1.R
# NLSProbDescription: { The Chloride data frame has 54 rows and 2 columns representing measurements of the chloride
# ion concentration in blood cells suspended in a salt solution.
# The two columns are:
# `conc`: A numeric vector giving the chloride ion concentration (%).
# `time`: A numeric vector giving the time of the concentration measurement (min).
# }


# Use the Chloride data from NRAIA package

## DATA
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
NLStestdata <- data.frame(conc,time)

## STARTING VALUE
Asym = 50
prop=0.6
lrc = log(0.25)
NLSstart <-c(Asym = Asym, prop=prop, lrc = lrc) # a starting vector (named!)

## MODEL
NLSformula <- conc ~ Asym*(1 - prop*exp(-exp(lrc)*time))
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- "(formula=NLSformula, data=NLStestdata, start=NLSstart)"
output_nls <- eval(parse(text=paste("nls",NLSrunline))) # nls is our benchmark case
output_nlsj <- eval(parse(text=paste("nlsj::nlsj",NLSrunline))) # nlsj is the new nls

## Test expectations using testthat
library(testthat) # comment out later!!

#### TESTING nls VS nlsj
# SETTING TOLERANCE
epstol <- sqrt(.Machine$double.eps*100) # Can replace 100 with nls.control()$offset

# NLSout/expout has "m", "convInfo", "data", "call",
# "dataClasses", "control"

## testing m values:  "resid"      "fitted"     "formula"    "deviance"   "lhs"       
# "gradient"   "conv"       "incr"       "setVarying" "setPars"   
# "getPars"    "getAllPars" "getEnv"     "trace"      "Rmat"      
# "predict"

test_that("testing m objects",{ #FAILED
      # residuals
	expect_equal(as.vector(resid(output_nls)),
			 as.vector(resid(output_nlsj)),
		    tolerance=epstol*(max(abs(c(as.vector(resid(output_nls)),
					as.vector(resid(output_nlsj)))
					)) + epstol))

#	# fitted
#	expect_equal(as.vector(fitted(output_nls)),
#			 as.vector(fitted(output_nlsj)))
#	# formula
#	expect_equal(formula(output_nls),
#			 formula(output_nlsj))
	# deviance
	expect_equal(deviance(output_nls),
			 deviance(output_nlsj),
		    tolerance=epstol*(max(abs(c(deviance(output_nls),
					deviance(output_nlsj))
					)) + epstol))
	# gradient
	expect_equal( output_nls$m$gradient(),
			   attr(output_nlsj$m$resid(),"gradient"),
		    tolerance=epstol*(max(abs(c(output_nls$m$gradient(),
					 attr(output_nlsj$m$resid(),"gradient"))
					)) + epstol))
#	# conv
#	expect_equal( output_nls$m$conv(),
#			  output_nlsj$m$conv())
#	# incr
#	expect_equal( output_nls$m$incr(),
#			  output_nlsj$m$incr())
#	# getPars # difference between getAllPars adn getPars?
	expect_equal( output_nls$m$getPars(),
			  output_nlsj$m$getPars())
#	# getEnv
#	expect_equal( output_nls$m$igetEnv(),
#			  output_nlsj$m$getEnv())
#	# trace 
#	##expect_equal( output_nls$m$trace(), ## Not run as it prints(devaince,conv,pars)
#	##		  output_nlsj$m$trace())	
	# Rmat
	expect_equal( output_nls$m$Rmat(),
			  output_nlsj$m$Rmat(),
		    tolerance=epstol*(max(abs(c(output_nls$m$Rmat(),
					output_nls$m$Rmat())
					)) + epstol))

#	# predict
#	expect_equal( output_nls$m$predict(),
#			  output_nlsj$m$predict())
	}
)

## testing control #FAILED
#test_that("testing control list items",{
#		expect_equal(output_nls$control,
#				 output_nlsj$control)
#	}
#)

# testing convInfo # FAILED
test_that("testing conInfo list items",{
		expect_equal(as.numeric(output_nls$convInfo$isConv),
			 as.numeric(output_nlsj$convInfo))
	}
)

#rm(conc,time, Asym, prop, lrc)

#rm("NLSformula","NLSrunline","NLSstart","NLStestdata",
#	"NLSupper","NLSlower","output_nls","output_nlsj","epstol")

print("End of test file 'Chloride_1.R' ")
#-----------------------------------------#