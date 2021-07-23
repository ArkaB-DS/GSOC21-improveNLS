# NLSProbName: Tetra_1.R
# NLSProbDescription: {The Hobbs weed infestation problem to estimate a 
#    3-parameter logistic S curve in its unscaled form from a reasonably
#    easy starting point of (200, 50, 0.3)
# }


## DATA
y=c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
          38.558, 50.156, 62.948, 75.995, 91.972)
tt=1:12
NLStestdata <- data.frame(y,tt)

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
				
rm(y, tt)

rm("NLSformula","NLSrunline","NLSstart","NLStestdata",
"NLSupper","NLSlower","output_nls","output_nlsj","output_nlxb",
"output_nlsLM","epstol","tolerance")

print("End of test file 'Hobbs_1.R' ")

#-----------------------------------------#
