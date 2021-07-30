# NLSProbName: Lipo_1.R
# NLSProbDescription: {The Lipo data frame has 12 rows and 2 columns of lipoprotein concentrations over time.
# The two columns are:This data frame contains the following columns:
# `time`: a numeric vector giving the time of the concentration measurement (hr)
# `conc`:  a numeric vector of concentrations.
# }


# Use the Isom data from NRAIA package

## DATA
time=c( 0.5,  1.0,  1.5,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0,  8.0,  9.0, 10.0)
conc = c( 46.10, 25.90, 17.00, 12.10,  7.22,  4.51,  3.19,  2.40,  1.82,  1.41,  1.00,
  		    0.94)
NLStestdata <- data.frame(time,conc)

## STARTING VALUE
lrc1=1/4
lrc2=-2
A1=100
A2=150
NLSstart <-c(lrc1=lrc1,lrc2=lrc2,A1=A1,A2=A2) # a starting vector (named!)

## MODEL
NLSformula <-conc ~ A1*exp(-exp(lrc1)*time)+A2*exp(-exp(lrc2)*time)
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

#rm(conc,time)

#rm("NLSformula","NLSrunline","NLSstart","NLStestdata",
#	"NLSupper","NLSlower","output_nls","output_nlsj","epstol")

print("End of test file 'Lipo_1.R' ")
#-----------------------------------------#
