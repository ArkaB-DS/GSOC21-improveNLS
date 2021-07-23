# NLSProbName: Chloride_1.R
# NLSProbDescription: { The Chloride data frame has 54 rows and 2 columns representing measurements of the chloride
ion concentration in blood cells suspended in a salt solution.
# The two columns are:
# `conc`: A numeric vector giving the chloride ion concentration (%).
# `time`: A numeric vector giving the time of the concentration measurement (min).
# }


# Use the Chloride data from NRAIA package

## DATA
NLStestdata <- data.frame(
	conc = c(17.3, 17.6, 17.9, 18.3, 18.5, 18.9, 19.0, 19.3, 19.8,
		  19.9, 20.2, 20.5, 20.6, 21.1, 21.5, 21.9, 22.0, 22.3,
		  22.6, 22.8, 23.0, 23.2, 23.4, 23.7, 24.0, 24.2, 24.5,
		  25.0, 25.4, 25.5, 25.9, 25.9, 26.3, 26.2, 26.5, 26.5,
		  26.6, 27.0, 27.0, 27.0, 27.0, 27.3, 27.8, 28.1, 28.1,
		  28.1, 28.4, 28.6, 29.0, 29.2, 29.3, 29.4, 29.4, 29.4),
	time = c(2.45, 2.55, 2.65, 2.75, 2.85, 2.95, 3.05, 3.15, 3.25,
		   3.35, 3.45, 3.55, 3.65, 3.75, 3.85, 3.95, 4.05, 4.15,
		   4.25, 4.35, 4.45, 4.55, 4.65, 4.75, 4.85, 4.95, 5.05,
		   5.15, 5.25, 5.35, 5.45, 5.55, 5.65, 5.75, 5.85, 5.95, 
		   6.05, 6.15, 6.25, 6.35, 6.45, 6.55, 6.65, 6.75, 6.85,
		   6.95, 7.05, 7.15, 7.25, 7.35, 7.45, 7.55, 7.65, 7.75)
	)

## STARTING VALUE

NLSstart <-c(Asym = 50, prop=0.6, lrc = log(0.25)) # a starting vector (named!)

## MODEL
NLSformula <- conc ~ Asym*(1 - prop*exp(-exp(lrc)*time))
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- "(formula=NLSformula, data=NLStestdata, start=NLSstart)"
output_nls <- eval(parse(text=paste("nls",NLSrunline))) # nls is our benchmark case
output_nlsj <- eval(parse(text=paste("nlsj::nlsj",NLSrunline))) # nlsj is the new nls
output_nlsLM <- eval(parse(text=paste("minpack.lm::nlsLM",NLSrunline))) # nlsLM
output_nlxb <- eval(parse(text=paste("nlsr::nlxb",NLSrunline))) # nlxb

# so, if we include nlsr::nlxb, minpack.lm::nlsLM, we can follow the above
# nomenclature as define output_nlxb, output_nlsLM

## Test expectations using testthat
library(testthat) # not required finally!

#### TESTING nls VS nlsj

#library(nlsj) # not required finally!

# NLSout/expout has "m", "convInfo", "data", "call",
# "dataClasses", "control"

## testing m values:  "resid"      "fitted"     "formula"    "deviance"   "lhs"       
# "gradient"   "conv"       "incr"       "setVarying" "setPars"   
# "getPars"    "getAllPars" "getEnv"     "trace"      "Rmat"      
# "predict"

test_that("testing m objects",{ #FAILED
      # residuals
	expect_equal(as.vector(resid(output_nls)),
			 as.vector(resid(output_nlsj)))
	# fitted
	expect_equal(as.vector(fitted(output_nls)),
			 as.vector(fitted(output_nlsj)))
	# formula
	expect_equal(formula(output_nls),
			 formula(output_nlsj))
	# deviance
	expect_equal(deviance(output_nls),
			 deviance(output_nlsj))
	# gradient
	expect_equal( output_nls$m$gradient(),
			  output_nlsj$m$gradient())
	# conv
	expect_equal( output_nls$m$conv(),
			  output_nlsj$m$conv())
	# incr
	expect_equal( output_nls$m$incr(),
			  output_nlsj$m$incr())
	# getAllPars # difference between getAllPars adn getPars?
	expect_equal( output_nls$m$getAllPars(),
			  output_nlsj$m$getAllPars())
	# getEnv
	expect_equal( output_nls$m$igetEnv(),
			  output_nlsj$m$getEnv())
	# trace 
	##expect_equal( output_nls$m$trace(), ## Not run as it prints(devaince,conv,pars)
	##		  output_nlsj$m$trace())	
	# Rmat
	expect_equal( output_nls$m$Rmat(),
			  output_nlsj$m$Rmat())
	# predict
	expect_equal( output_nls$m$predict(),
			  output_nlsj$m$predict())
	}
)

# testing control #FAILED
test_that("testing control list items",{
		expect_equal(output_nls$control,
				 output_nlsj$control)
	}
)

# testing convInfo # FAILED
test_that("testing conInfo list items",{
		expect_equal(output_nls$convInfo,
				 output_nlsj$convInfo)
	}
)

#### TESTING nls VS minpacl.lm::nlsLM

## testing m values:  # PASSED
test_that("testing m objects",{
      # residuals
	expect_equal(as.vector(resid(output_nls)),
			 as.vector(resid(output_nlsLM)))
	# fitted
	expect_equal(as.vector(fitted(output_nls)),
			 as.vector(fitted(output_nlsLM)))
	# formula
	expect_equal(formula(output_nls),
			 formula(output_nlsLM))
	# deviance
	expect_equal(deviance(output_nls),
			 deviance(output_nlsLM))
	# gradient
	expect_equal( output_nls$m$gradient(),
			  output_nlsLM$m$gradient())
	# conv
	expect_equal( output_nls$m$conv(),
			  output_nlsLM$m$conv())
	# incr
	expect_equal( output_nls$m$incr(),
			  output_nlsLM$m$incr())
	# getAllPars # difference between getAllPars adn getPars?
	expect_equal( output_nls$m$getAllPars(),
			  output_nlsLM$m$getAllPars())
	# getEnv
	expect_equal( output_nls$m$getEnv(),
			  output_nlsLM$m$getEnv())
	# trace 
	##expect_equal( output_nls$m$trace(), ## Not run as it prints(devaince,conv,pars)
	##		  output_nlsLM$m$trace())	
	# Rmat
	expect_equal( output_nls$m$Rmat(),
			  output_nlsLM$m$Rmat())
	# predict
	expect_equal( output_nls$m$predict(),
			  output_nlsLM$m$predict())
	}
)

# testing control #FAILED
test_that("testing control list items",{
		expect_equal(output_nls$control,
				 output_nlsLM$control)
	}
)

# testing convInfo ## FAILED
test_that("testing conInfo list items",{
		expect_equal(output_nls$convInfo,
				 output_nlsLM$convInfo)
	}
)

## TESTING nls VS nlsr::nlxb

# testing coefficients # FAILED----> need to change tolerance - suggestions??
test_that("testing parameter estimates",{
		expect_equal(as.vector(output_nls$m$getAllPars()),
				as.vector(coefficients(output_nlxb)))
	}
)
# testing residuals # FAILED----> signs all opposite??
test_that("testing residuals",{
		expect_equal(as.vector(residuals(output_nls)),
				as.vector(output_nlxb$resid))
	}
)
# testing jacobian # FAILED----> should I change dimnames here? then it should pass
test_that("testing jacobian",{
		expect_equal(as.matrix(output_nls$m$gradient()),
				 output_nlxb$jacobian)
	}
)
# testing formula  #PASSED
test_that("testing formula",{
		expect_equal(formula(output_nls),
			 	formula(output_nlxb))
	}
)
# testing weights  #PASSED
test_that("testing weights",{
		expect_equal(weights(output_nls),
			 	weights(output_nlxb))
	}
)
# testing ssquares  #PASSED
test_that("testing sum of squares",{
		expect_equal(sum(as.vector(residuals(output_nls))^2),
			 	output_nlxb$ssquares)
	}
)

 
# How do I get feval, jeval, lower, upper in nls? Also, what is maskidx?Not documented in the help file of nlxb
				





print("End of test file 'Chloride_1.R' ")
#-----------------------------------------#

## To delete the new variables due to this test file from the workspace 

# ls() ## see what variables are there in thr workspace

## "demand", "NLSformula", "NLSlower", "NLSrunline", 
## "NLSstart", "NLStestdata", "NLSupper", "output_nls" 
## "output_nlsj", "time"     
## note that "demand" and "time" won't be there in other test cases, generally
## SUGGESTION: use NLSdata <- data.frame(time=c(1,2),demand=c(3,4))

## then we can use rm("NLSformula","NLSrunline","NLSstart","NLStestdata",
## "NLSupper","NLSlower","output_nls","output_nlsj","output_nlxb",
## "output_nlsLM")

## What do you think Dr. Nash? Further, the idea of uisng BOD2_x, where x represents a number 
## would be better. Tags would indicate the differences in the x's










