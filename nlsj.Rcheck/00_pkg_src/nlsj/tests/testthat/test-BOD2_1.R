# NLSProbName: BOD2_1.R
# NLSProbDescription: { The BOD2 data frame has 8 rows and 2 columns
# giving the biochemical oxygen demand versus time in an evaluation of water quality.
# The two columns are:
# `time`: A numeric vector giving the time of the measurement (days).
# `demand`: A numeric vector giving the biochemical oxygen demand (mg/l).
# }


# Use the BOD2 data from NRAIA package

## DATA
time <- c( 1,  2,  3,  4,  5,  7,  9, 11)
demand <- c(0.47, 0.74, 1.17, 1.42, 1.60, 1.84, 2.19, 2.17)
NLStestdata <- data.frame(demand,time) 

## STARTING VALUE
A = 2.2
lrc = log(0.25)
NLSstart <- c(A=A,lrc=lrc) # a starting vector (named!)

## MODEL
NLSformula <- demand ~ A*(1-exp(-exp(lrc)*time))
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- "(formula=NLSformula, data=NLStestdata, start=NLSstart)"
output_nls <- eval(parse(text=paste("nls",NLSrunline))) # nls is our benchmark case
output_nlsj <- eval(parse(text=paste("nlsj::nlsj",NLSrunline))) # nlsj is the new nls

## Test expectations using testthat
#library(testthat) # comment out later!!

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
#	## incr
#	#expect_equal( output_nls$m$incr(),
#	#		  output_nlsj$m$incr())
	# getPars # difference between getAllPars adn getPars?
	expect_equal( output_nls$m$getPars(),
			  output_nlsj$m$getPars())
#	## getEnv
#	#expect_equal( output_nls$m$igetEnv(),
#	#		  output_nlsj$m$getEnv())
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

#rm(demand,time,A,lrc)

#rm("NLSformula","NLSrunline","NLSstart","NLStestdata",
#	"NLSupper","NLSlower","output_nls","output_nlsj","epstol")

print("End of test file 'BOD2_1.R' ")
#-----------------------------------------#
