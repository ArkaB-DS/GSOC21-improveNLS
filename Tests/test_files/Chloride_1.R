# NLSProbName: Chloride_1.R
# NLSProbDescription: { The Chloride data frame has 54 rows and 2 columns representing measurements of the chloride
# ion concentration in blood cells suspended in a salt solution.
# The two columns are:
# `conc`: A numeric vector giving the chloride ion concentration (%).
# `time`: A numeric vector giving the time of the concentration measurement (min).
# }


# Use the Chloride data from NRAIA package
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
## DATA
NLStestdata <- data.frame(conc,time)

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
NLSmatrix <- as.data.frame(matrix(NA,nrow=3,ncol=6,dimnames=list(c("nlsj","nlsLM","nlxb"),
						c("residual","gradient","deviance","getPars","Rmat","isConv"))
			)
	)	

#### TESTING nls VS nlsj
# SETTING TOLERANCE

epstol <- sqrt(.Machine$double.eps*100) # Can replace 100 with nls.control()$offset

# NLSout/expout has "m", "convInfo", "data", "call",
# "dataClasses", "control"

## testing m values:  "resid"      "fitted"     "formula"    "deviance"   "lhs"       
# "gradient"   "conv"       "incr"       "setVarying" "setPars"   
# "getPars"    "getAllPars" "getEnv"     "trace"      "Rmat"      
# "predict"

      # residuals
NLSmatrix[1,1] <-	all.equal(as.vector(resid(output_nls)),
		    as.vector(resid(output_nlsj)),
		    tolerance=epstol*(max(abs(c(as.vector(resid(output_nls)),
					as.vector(resid(output_nlsj)))
					)) + epstol))
	# fitted
	all.equal(as.vector(fitted(output_nls)),
		    as.vector(fitted(output_nlsj)),
		    tolerance=epstol*(max(abs(c(as.vector(fitted(output_nls)),
					as.vector(fitted(output_nlsj)))
					)) + epstol))
#	# formula
#	all.equal(formula(output_nls),
#			 formula(output_nlsj))
	# deviance
NLSmatrix[1,3] <- all.equal(deviance(output_nls),
			 deviance(output_nlsj))
	# gradient
NLSmatrix[1,2] <- all.equal( output_nls$m$gradient(),
			  attr(output_nlsj$m$resid(),"gradient"))
#	# conv
#	all.equal( output_nls$m$conv(),
#			  output_nlsj$m$conv())
#	# incr
#	all.equal( output_nls$m$incr(),
#			  output_nlsj$m$incr())
	# getPars # difference between getAllPars adn getPars?
NLSmatrix[1,4] <- all.equal( output_nls$m$getPars(),
			  output_nlsj$m$getPars())
#	# getEnv
#	all.equal( output_nls$m$igetEnv(),
#			  output_nlsj$m$getEnv())
	# Rmat
NLSmatrix[1,5] <- all.equal( output_nls$m$Rmat(),
			  output_nlsj$m$Rmat())
#	# predict
#	all.equal( output_nls$m$predict(),
#			  output_nlsj$m$predict())

# testing control #FAILED
#	all.equal(output_nls$control,
#			 output_nlsj$control)

# testing convInfo # FAILED
NLSmatrix[1,6] <- all.equal(as.numeric(output_nls$convInfo$isConv),
			 as.numeric(output_nlsj$convInfo))

#### TESTING nls VS minpack.lm::nlsLM

## testing m values:  # PASSED
      # residuals
NLSmatrix[2,1] <-	all.equal(as.vector(resid(output_nls)),
			 as.vector(resid(output_nlsLM)))
	# fitted
	all.equal(as.vector(fitted(output_nls)),
			 as.vector(fitted(output_nlsLM)))
#	# formula
#	all.equal(formula(output_nls),
#			 formula(output_nlsLM))
	# deviance
NLSmatrix[2,3] <-	all.equal(deviance(output_nls),
			 deviance(output_nlsLM))
	# gradient
NLSmatrix[2,2] <-	all.equal( output_nls$m$gradient(),
			  output_nlsLM$m$gradient())
#	# conv
#	all.equal( output_nls$m$conv(),
#			  output_nlsLM$m$conv())
#	# incr
#	all.equal( output_nls$m$incr(),
#			  output_nlsLM$m$incr())
# getPars # difference between getAllPars adn getPars?
NLSmatrix[2,4] <- all.equal( output_nls$m$getAllPars(),
			  output_nlsLM$m$getAllPars())
#	# getEnv
#	all.equal( output_nls$m$getEnv(),
#			  output_nlsLM$m$getEnv())
	# Rmat
NLSmatrix[2,5] <- all.equal( output_nls$m$Rmat(),
			  output_nlsLM$m$Rmat())
#	# predict
#	all.equal( output_nls$m$predict(),
#			  output_nlsLM$m$predict())

# testing control #FAILED
#	all.equal(output_nls$control,
#			 output_nlsLM$control)

# testing convInfo ## FAILED
NLSmatrix[2,6] <- all.equal(output_nls$convInfo$isConv,
			 output_nlsLM$convInfo$isConv)
	all.equal(output_nls$convInfo$finIter,
			 output_nlsLM$convInfo$finIter)

## TESTING nls VS nlsr::nlxb

# testing coefficients # FAILED----> need to change tolerance - suggestions??
NLSmatrix[3,4] <-	all.equal(as.vector(output_nls$m$getAllPars()),
			as.vector(coefficients(output_nlxb)),
			tolerance=epstol*(max(abs(c(as.vector(output_nls$m$getAllPars()),
					as.vector(coefficients(output_nlxb)))
					)) + epstol))
# testing residuals # FAILED----> signs all opposite??
NLSmatrix[3,1] <-	all.equal(as.vector(residuals(output_nls)),
			-1*as.vector(output_nlxb$resid),
			tolerance=epstol*(max(abs(c(as.vector(residuals(output_nls)),
					-1*as.vector(output_nlxb$resid))
					)) + epstol))

# testing jacobian # FAILED----> should I change dimnames here? then it should pass
NLSmatrix[3,2] <-	all.equal(as.numeric(as.matrix(output_nls$m$gradient())),
			 as.numeric(output_nlxb$jacobian))
# testing formula  #PASSED
#	all.equal(formula(output_nls),
#			 	formula(output_nlxb))
# testing weights  #PASSED
#	all.equal(weights(output_nls),
#			 	weights(output_nlxb))
# testing ssquares  #PASSED
	all.equal(sum(as.vector(residuals(output_nls))^2),
		 	output_nlxb$ssquares)
 
# How do I get feval, jeval, lower, upper in nls? Also, what is maskidx?Not documented in the help file of nlxb
				
rm(conc, time)

rm("NLSformula","NLSrunline","NLSstart","NLStestdata",
"NLSupper","NLSlower","output_nls","output_nlsj","output_nlxb",
"output_nlsLM","epstol","tolerance")

print("End of test file 'Chloride_1.R' ")

#-----------------------------------------#
