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
rate=c(3.541,  2.397,  6.694,  4.722,  0.593,  0.268,  2.797,
		 2.451,  3.196,  2.021,  0.896,  5.084,  5.686,  1.193,
		 2.648,  3.303,  3.054,  3.302,  1.271, 11.648,  2.002,
		 9.604,  7.754, 11.590)
hyd = c(205.8, 404.8, 209.7, 401.6, 224.9, 402.6, 212.7, 406.2, 133.3, 470.9, 300.0,
		  301.6, 297.3, 314.0, 305.7, 300.1, 305.4, 305.2, 300.1, 106.6, 417.2, 251.0,
		  250.3, 145.1)
iso.pen = c(37.1,  36.3,  49.4,  44.9, 116.3, 128.9, 134.4, 134.9,  87.6,  86.9,  81.7,
			101.7,  10.5, 157.1,  86.0,  90.2,  87.4,  87.0,  66.4,  33.0,  32.9,  41.5,
			14.7,  50.2)
n.pent = c( 90.9,  92.9, 174.9, 187.2,  92.7, 102.2, 186.9, 192.6, 140.8, 144.2,  68.3,
			214.6, 142.2, 146.7, 142.0, 143.7, 141.1, 141.5,  83.0, 209.6,  83.9, 294.4,
			148.0, 291.0)
	
NLStestdata <- data.frame(rate,hyd,iso.pen,n.pent)

## STARTING VALUE
b2 = 0.1
b3 = 0.1
b4 = 0.1
NLSstart <-c(b2 = b2, b3 = b3, b4 = b4) # a starting vector (named!)
NLSstart2 <-c(b1=1, b2 = b2, b3 = b3, b4 = b4) # a starting vector (named!)

## MODEL
NLSformula <- rate ~ b3*(n.pent - iso.pen/1.632)/(1+b2*hyd+b3*n.pent+b4*iso.pen)
NLSformula2 <- rate ~ b1 * b3*(n.pent - iso.pen/1.632)/(1+b2*hyd+b3*n.pent+b4*iso.pen)
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- '(formula=NLSformula, data=NLStestdata, start=NLSstart,algorithm="plinear")'
NLSrunline2 <- '(formula=NLSformula2, data=NLStestdata, start=NLSstart2, trace=TRUE)'
NLSrunline3 <- '(formula=NLSformula2, data=NLStestdata, start=NLSstart2, trace=TRUE, algorithm="marquardt")'
output_nls <- eval(parse(text=paste("nls",NLSrunline))) # nls is our benchmark case
output_nlsj <- eval(parse(text=paste("nlsj::nlsj",NLSrunline))) # nlsj is the new nls
output_nlsj2 <- eval(parse(text=paste("nlsj",NLSrunline2)))
output_nlsj3 <- eval(parse(text=paste("nlsj",NLSrunline3)))
tmp <- readline("continue")
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


print("End of test file 'Isom_1.R' file.")
#-----------------------------------------#
