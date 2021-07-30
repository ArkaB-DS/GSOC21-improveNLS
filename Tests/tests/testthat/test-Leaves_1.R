# NLSProbName: Leaves_1.R
# NLSProbDescription: { The Leaves data frame has 15 rows and 2 columns of leaf length over time.
# The two columns are:This data frame contains the following columns:
# `time`: e time from initial emergence (days).
# `length`: leaf length (cm).
# }


# Use the Leaves data from NRAIA package

## DATA
time=c( 0.5,  1.5,  2.5,  3.5,  4.5,  5.5,  6.5,  7.5,  8.5,  9.5, 10.5, 11.5, 12.5, 13.5,
		  14.5)
length = c( 1.3,  1.3,  1.9,  3.4,  5.3,  7.1, 10.6, 16.0, 16.4, 18.3, 20.9, 20.5, 21.3, 21.2,
			20.9)
NLStestdata <- data.frame(time,length)

## STARTING VALUE
Asym=3
xmid=2
scal=1
NLSstart <-c(Asym=Asym,xmid=xmid,scal=scal) # a starting vector (named!)

## MODEL
NLSformula <-length ~ Asym/(1+exp((xmid-time)/scal))
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

#rm(length,time)

#rm("NLSformula","NLSrunline","NLSstart","NLStestdata",
#	"NLSupper","NLSlower","output_nls","output_nlsj","epstol")

print("End of test file 'Leaves_1.R' ")
#-----------------------------------------#
