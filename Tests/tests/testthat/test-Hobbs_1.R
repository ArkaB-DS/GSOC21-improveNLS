# NLSProbName: Tetra_1.R
# NLSProbDescription: {The Hobbs weed infestation problem to estimate a 
#    3-parameter logistic S curve in its unscaled form from a reasonably
#    easy starting point of (200, 50, 0.3)
# }


## DATA
NLStestdata <- data.frame(
	y=c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972),
	tt=1:12
	)

## STARTING VALUE

NLSstart <-c(b1=200, b2=50, b3=0.3) # a starting vector (named!)

## MODEL
NLSformula <- y ~ b1/(1+b2*exp(-b3*tt))
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










