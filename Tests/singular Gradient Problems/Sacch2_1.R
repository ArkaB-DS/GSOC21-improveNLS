# NLSProbName: Sacch2_1.R
# NLSProbDescription: {The Sacch2 data frame has 10 rows and 2 columns from an experiment on the pharmacokinetics of
# saccharin.
# The two columns are:This data frame contains the following columns:
# `time`: a numeric vector giving the time since drug administration (min).
# `conc`: a numeric vector giving the observed concentration of saccharin.

# }


# Use the Sacch2 data from NRAIA package

## DATA
time=c(0,   5,  15,  30,  45,  60,  75,  90, 105, 120)
conc = c(0.0, 184.3, 102.0,  50.5,  24.9,  14.1,   8.0,   5.7,   4.0,   2.9)
	
NLStestdata <- data.frame(time,conc)

## STARTING VALUE

NLSstart <- c(Dose=1,lKa=13,lKe=17,lCl=7) # a starting vector (named!)

## MODEL
NLSformula <-conc ~ Dose * exp(lKe+lKa-lCl) * (exp(-exp(lKe)*time) - exp(-exp(lKa)*time))/(exp(lKa) - exp(lKe))
NLSlower <- NULL
NLSupper <- NULL
NLSrunline <- "(formula=NLSformula, data=NLStestdata, start=NLSstart)"
output_nls <- eval(parse(text=paste("nls",NLSrunline))) # nls is our benchmark case
output_nlsj <- eval(parse(text=paste("nlsj::nlsj",NLSrunline))) # nlsj is the new nls
output_nlsLM <- eval(parse(text=paste("minpack.lm::nlsLM",NLSrunline))) # nlsLM
output_nlxb <- eval(parse(text=paste("nlsr::nlxb",NLSrunline))) # nlxb

#### TESTING nls VS nlsj

# NLSout/expout has "m", "convInfo", "data", "call",
# "dataClasses", "control"

## testing m values:  "resid"      "fitted"     "formula"    "deviance"   "lhs"       
# "gradient"   "conv"       "incr"       "setVarying" "setPars"   
# "getPars"    "getAllPars" "getEnv"     "trace"      "Rmat"      
# "predict"

      # residuals
	all.equal(as.vector(resid(output_nls)),
			 as.vector(resid(output_nlsj)))
	# fitted
	all.equal(as.vector(fitted(output_nls)),
			 as.vector(fitted(output_nlsj)))
#	# formula
#	all.equal(formula(output_nls),
#			 formula(output_nlsj))
	# deviance
	all.equal(deviance(output_nls),
			 deviance(output_nlsj))
	# gradient
	all.equal( output_nls$m$gradient(),
			  attr(output_nlsj$m$resid(),"gradient"))
#	# conv
#	all.equal( output_nls$m$conv(),
#			  output_nlsj$m$conv())
#	# incr
#	all.equal( output_nls$m$incr(),
#			  output_nlsj$m$incr())
	# getPars # difference between getAllPars adn getPars?
	all.equal( output_nls$m$getPars(),
			  output_nlsj$m$getPars())
#	# getEnv
#	all.equal( output_nls$m$igetEnv(),
#			  output_nlsj$m$getEnv())
	# Rmat
	all.equal( output_nls$m$Rmat(),
			  output_nlsj$m$Rmat())
#	# predict
#	all.equal( output_nls$m$predict(),
#			  output_nlsj$m$predict())

# testing control #FAILED
#	all.equal(output_nls$control,
#			 output_nlsj$control)

# testing convInfo # FAILED
	all.equal(as.numeric(output_nls$convInfo$isConv),
			 as.numeric(output_nlsj$convInfo))

#### TESTING nls VS minpacl.lm::nlsLM

## testing m values:  # PASSED
      # residuals
	all.equal(as.vector(resid(output_nls)),
			 as.vector(resid(output_nlsLM)))
	# fitted
	all.equal(as.vector(fitted(output_nls)),
			 as.vector(fitted(output_nlsLM)))
#	# formula
#	all.equal(formula(output_nls),
#			 formula(output_nlsLM))
	# deviance
	all.equal(deviance(output_nls),
			 deviance(output_nlsLM))
	# gradient
	all.equal( output_nls$m$gradient(),
			  output_nlsLM$m$gradient())
#	# conv
#	all.equal( output_nls$m$conv(),
#			  output_nlsLM$m$conv())
#	# incr
#	all.equal( output_nls$m$incr(),
#			  output_nlsLM$m$incr())
	# getPars # difference between getAllPars adn getPars?
	all.equal( output_nls$m$getAllPars(),
			  output_nlsLM$m$getAllPars())
#	# getEnv
#	all.equal( output_nls$m$getEnv(),
#			  output_nlsLM$m$getEnv())
	# Rmat
	all.equal( output_nls$m$Rmat(),
			  output_nlsLM$m$Rmat())
#	# predict
#	all.equal( output_nls$m$predict(),
#			  output_nlsLM$m$predict())

# testing control #FAILED
#	all.equal(output_nls$control,
#			 output_nlsLM$control)

# testing convInfo ## FAILED
	all.equal(output_nls$convInfo$isConv,
			 output_nlsLM$convInfo$isConv)
	all.equal(output_nls$convInfo$finIter,
			 output_nlsLM$convInfo$finIter)

## TESTING nls VS nlsr::nlxb

# testing coefficients # FAILED----> need to change tolerance - suggestions??
	all.equal(as.vector(output_nls$m$getAllPars()),
			as.vector(coefficients(output_nlxb)))
# testing residuals # FAILED----> signs all opposite??
	all.equal(as.vector(residuals(output_nls)),
			-1*as.vector(output_nlxb$resid))
# testing jacobian # FAILED----> should I change dimnames here? then it should pass
	all.equal(as.numeric(as.matrix(output_nls$m$gradient())),
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
"output_nlsLM")

print("End of test file 'Sacch2_1.R' ")