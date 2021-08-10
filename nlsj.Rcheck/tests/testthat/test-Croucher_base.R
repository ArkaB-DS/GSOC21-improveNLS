# Croucher-base.R -- https://walkingrandomly.com/?p=5254
# NLSProbName: Croucher-base
# NLSProbDescription: {This is a fairly simple 2-parameter problem. This file is the "base"
#     from which others in the Croucher family are built.}

xdata <- c(-2,-1.64,-1.33,-0.7,0,0.45,1.2,1.64,2.32,2.9)
ydata <- c(0.699369,0.700462,0.695354,1.03905,1.97389,2.41143,1.91091,0.919576,-0.730975,-1.42001)
p1<- 1
p2<-0.2
NLSformula <- ydata ~ p1*cos(p2*xdata) + p2*sin(p1*xdata)
NLSstart<-list(p1=p1,p2=p2) # This is the default start given in the reference.
NLStestdata<-data.frame(xdata, ydata)
NLSweights <- rep(0.25, length(xdata))
NLSsubset<-1:8
NLSlower<-c(0,0)
NLSupper<-c(1.5, 1.5)
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

print("End of test file 'Croucher_base.R' ")
