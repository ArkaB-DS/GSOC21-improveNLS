# NLSProbName: Isom_1.R
# NLSProbDescription: { The Isom data frame has 24 rows and 4 columns from an isomerization experiment.
# The four columns are:This data frame contains the following columns:
# `hyd`: partial pressure of hydrogen (psia).
# `n.pent`: partial pressure of n-pentane (psia).
# `iso.pen`: partial pressure of isopentane (psia).
# `rate`: reaction rate for isomerization of n-pentane to isopentane (1/hr).
# }


# Use the Isom data from NRAIA package

## DATA
NLStestdata <- data.frame(
	rate=c(3.541,  2.397,  6.694,  4.722,  0.593,  0.268,  2.797,
		 2.451,  3.196,  2.021,  0.896,  5.084,  5.686,  1.193,
		 2.648,  3.303,  3.054,  3.302,  1.271, 11.648,  2.002,
		 9.604,  7.754, 11.590),
	hyd = c(205.8, 404.8, 209.7, 401.6, 224.9, 402.6, 212.7, 406.2, 133.3, 470.9, 300.0,
		  301.6, 297.3, 314.0, 305.7, 300.1, 305.4, 305.2, 300.1, 106.6, 417.2, 251.0,
		  250.3, 145.1),
	iso.pen = c(37.1,  36.3,  49.4,  44.9, 116.3, 128.9, 134.4, 134.9,  87.6,  86.9,  81.7,
			101.7,  10.5, 157.1,  86.0,  90.2,  87.4,  87.0,  66.4,  33.0,  32.9,  41.5,
			14.7,  50.2),
	n.pent = c( 90.9,  92.9, 174.9, 187.2,  92.7, 102.2, 186.9, 192.6, 140.8, 144.2,  68.3,
			214.6, 142.2, 146.7, 142.0, 143.7, 141.1, 141.5,  83.0, 209.6,  83.9, 294.4,
			148.0, 291.0)
	)

## STARTING VALUE

NLSstart <-c(b2 = 0.1, b3 = 0.1, b4 = 0.1) # a starting vector (named!)

## MODEL
NLSformula <- rate ~ b3*(n.pent - iso.pen/1.632)/(1+b2*hyd+b3*n.pent+b4*iso.pen)
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- '(formula=NLSformula, data=NLStestdata, start=NLSstart,algorithm="plinear")'
output_nls <- eval(parse(text=paste("nls",NLSrunline))) # nls is our benchmark case
output_nls2 <- eval(parse(text=paste("nls",NLSrunline)))
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










